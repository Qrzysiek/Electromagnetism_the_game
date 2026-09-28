// Time-dependent field views on the GPU (visual only, f32): antennas, plane waves and the
// retarded (Liénard–Wiechert) fields of moving charges, plus for the total view a static
// field computed on the CPU (f64, exact, including metal) and stored on a grid. Every
// pixel is evaluated at the animation time; all times are relative to it (τ = t − t_now),
// so f32 keeps full precision over long flights. Mirrors the CPU code in radiation.rs
// (`sample`, `charges`), which is still used for the colour scales and the E arrows.

#import bevy_sprite::mesh2d_vertex_output::VertexOutput

struct Params {
    // min.x, min.y, size.x, size.y of the arena.
    area: vec4<f32>,
    // c (0: infinite), B saturation, E saturation, dynamic range (10^decades) of the
    // logarithmic scale or gain of the linear one.
    scales: vec4<f32>,
    // Static grid width, height, 1 if present; number of moving charges.
    grid: vec4<u32>,
    // Number of antennas, of waves; flags: 1 radiation part only, 2 only what the
    // quasi-static beam interaction leaves out, 4 colour |E| (else B_z), 8 linear
    // colour scale; unused.
    counts: vec4<u32>,
};

@group(#{MATERIAL_BIND_GROUP}) @binding(0) var<uniform> params: Params;
// World-line samples, two per sample: (τ, x, y, vx), (vy, ax, ay, 0).
@group(#{MATERIAL_BIND_GROUP}) @binding(1) var<storage, read> samples: array<vec4<f32>>;
// Items: charges, three vec4 each: (offset, count, q, has_end), (τ_end, x_end, y_end,
// 1 if its charge stays where it was absorbed, 0 if drained), (magnetic moment, 0, 0, 0);
// then antennas, two vec4 each: (x, y, p0x, p0y), (ω, phase now, radius, 0); then waves,
// two vec4 each: (k̂x, k̂y, êx, êy), (E0, ω, phase now at x = 0, 0).
@group(#{MATERIAL_BIND_GROUP}) @binding(2) var<storage, read> items: array<vec4<f32>>;
// Static field on a grid of texel centres: (Ex, Ey, Bz, inside a body).
@group(#{MATERIAL_BIND_GROUP}) @binding(3) var<storage, read> statics: array<vec4<f32>>;

struct State {
    x: vec2<f32>,
    v: vec2<f32>,
    a: vec2<f32>,
};

fn sample_t(o: u32, i: u32) -> f32 {
    return samples[2u * (o + i)].x;
}

fn sample_state(o: u32, i: u32) -> State {
    let p = samples[2u * (o + i)];
    let q = samples[2u * (o + i) + 1u];
    return State(p.yz, vec2<f32>(p.w, q.x), q.yz);
}

// State of the world line at τ: cubic Hermite position, linear velocity and acceleration
// between samples; uniform motion before the first and after the last (as the CPU's
// SampledWorldline). `i` is the interval index (τ_i ≤ τ ≤ τ_{i+1}).
fn state_in(o: u32, i: u32, tau: f32) -> State {
    let s0 = sample_state(o, i);
    let s1 = sample_state(o, i + 1u);
    let t0 = sample_t(o, i);
    let h = sample_t(o, i + 1u) - t0;
    let s = clamp((tau - t0) / h, 0.0, 1.0);
    let h00 = (1.0 + 2.0 * s) * (1.0 - s) * (1.0 - s);
    let h10 = s * (1.0 - s) * (1.0 - s);
    let h01 = s * s * (3.0 - 2.0 * s);
    let h11 = s * s * (s - 1.0);
    let x = s0.x * h00 + s0.v * (h10 * h) + s1.x * h01 + s1.v * (h11 * h);
    return State(x, mix(s0.v, s1.v, s), mix(s0.a, s1.a, s));
}

fn uniform_from(s: State, t0: f32, tau: f32) -> State {
    return State(s.x + s.v * (tau - t0), s.v, vec2<f32>(0.0));
}

// Retarded time of charge `k` for the point `p` (now, τ = 0): the root of
// g(τ) = −c τ − |p − x(τ)|, which decreases strictly. Returns the state there.
fn retarded(o: u32, n: u32, p: vec2<f32>, c: f32) -> State {
    let first = sample_state(o, 0u);
    let t_first = sample_t(o, 0u);
    let g_first = -c * t_first - length(p - first.x);
    if (n == 1u || g_first < 0.0) {
        // Before the first sample: uniform motion, Newton from there.
        var tau = min(t_first, 0.0);
        for (var k = 0; k < 12; k = k + 1) {
            let s = uniform_from(first, t_first, tau);
            let d = p - s.x;
            let dist = max(length(d), 1e-6);
            let g = -c * tau - dist;
            let slope = -c + dot(d, s.v) / dist;
            tau = tau - g / slope;
        }
        return uniform_from(first, t_first, tau);
    }
    let last = sample_state(o, n - 1u);
    let t_last = sample_t(o, n - 1u);
    let g_last = -c * t_last - length(p - last.x);
    if (g_last >= 0.0) {
        // After the last sample: uniform motion.
        var tau = t_last;
        for (var k = 0; k < 12; k = k + 1) {
            let s = uniform_from(last, t_last, tau);
            let d = p - s.x;
            let dist = max(length(d), 1e-6);
            let g = -c * tau - dist;
            let slope = -c + dot(d, s.v) / dist;
            tau = tau - g / slope;
        }
        return uniform_from(last, t_last, tau);
    }
    // Binary search over the samples: g(τ_lo) ≥ 0 > g(τ_hi).
    var lo = 0u;
    var hi = n - 1u;
    loop {
        if (hi - lo <= 1u) {
            break;
        }
        let mid = (lo + hi) / 2u;
        let s = sample_state(o, mid);
        if (-c * sample_t(o, mid) - length(p - s.x) >= 0.0) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    // Inside the interval: bisection, then Newton.
    var a = sample_t(o, lo);
    var b = sample_t(o, hi);
    for (var k = 0; k < 14; k = k + 1) {
        let m = 0.5 * (a + b);
        let s = state_in(o, lo, m);
        if (-c * m - length(p - s.x) >= 0.0) {
            a = m;
        } else {
            b = m;
        }
    }
    var tau = 0.5 * (a + b);
    for (var k = 0; k < 2; k = k + 1) {
        let s = state_in(o, lo, tau);
        let d = p - s.x;
        let dist = max(length(d), 1e-6);
        let g = -c * tau - dist;
        let slope = -c + dot(d, s.v) / dist;
        tau = clamp(tau - g / slope, sample_t(o, lo), sample_t(o, hi));
    }
    return state_in(o, lo, tau);
}

// State now (τ = 0).
fn present(o: u32, n: u32) -> State {
    let t_first = sample_t(o, 0u);
    if (n == 1u || 0.0 <= t_first) {
        return uniform_from(sample_state(o, 0u), t_first, 0.0);
    }
    let t_last = sample_t(o, n - 1u);
    if (0.0 >= t_last) {
        return uniform_from(sample_state(o, n - 1u), t_last, 0.0);
    }
    var lo = 0u;
    var hi = n - 1u;
    loop {
        if (hi - lo <= 1u) {
            break;
        }
        let mid = (lo + hi) / 2u;
        if (sample_t(o, mid) <= 0.0) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return state_in(o, lo, 0.0);
}

struct Field {
    e: vec2<f32>,
    bz: f32,
};

fn cross_z(a: vec2<f32>, b: vec2<f32>) -> f32 {
    return a.x * b.y - a.y * b.x;
}

// Liénard–Wiechert fields at `p` of charge q whose retarded state is `s` (Jackson §14.1):
// E = q (n − β)(1 − β²)/(κ³R²) + (q/c) n × ((n − β) × β̇)/(κ³R), B = n × E / c.
fn lienard(q: f32, c: f32, p: vec2<f32>, s: State, radiation_only: bool) -> Field {
    let d = p - s.x;
    let dist = max(length(d), 1e-6);
    let n = d / dist;
    let beta = s.v / c;
    let beta_dot = s.a / c;
    let kappa = 1.0 - dot(n, beta);
    let k3 = kappa * kappa * kappa;
    let nb = n - beta;
    // n × (nb × β̇) in the plane: nb × β̇ is along z with value w; n × (w ẑ) = w (n.y, −n.x).
    let w = cross_z(nb, beta_dot);
    let e_rad = vec2<f32>(n.y * w, -n.x * w) * (q / (c * k3 * dist));
    var e = e_rad;
    if (!radiation_only) {
        e = e + nb * (q * (1.0 - dot(beta, beta)) / (k3 * dist * dist));
    }
    return Field(e, cross_z(n, e) / c);
}

// The quasi-static beam interaction's field of a charge now in state `s`: its past taken
// as constant acceleration while |a τ| ≤ 0.1 c, uniform before (beam::accelerated_fields).
fn accelerated(q: f32, c: f32, p: vec2<f32>, s: State) -> Field {
    let tau_lim = 0.1 * c / max(length(s.a), 1e-20);
    var tau = -length(p - s.x) / c;
    for (var k = 0; k < 20; k = k + 1) {
        let tc = max(tau, -tau_lim);
        let vc = s.v + s.a * tc;
        let xc = s.x + s.v * tc + s.a * (0.5 * tc * tc) + vc * (tau - tc);
        let d = p - xc;
        let dist = max(length(d), 1e-6);
        let g = -c * tau - dist;
        let slope = -c + dot(d, vc) / dist;
        tau = tau - g / slope;
    }
    let tc = max(tau, -tau_lim);
    let vc = s.v + s.a * tc;
    let xc = s.x + s.v * tc + s.a * (0.5 * tc * tc) + vc * (tau - tc);
    var ac = s.a;
    if (tau < -tau_lim) {
        ac = vec2<f32>(0.0);
    }
    return lienard(q, c, p, State(xc, vc, ac), false);
}

// Field of the moving charges at `p`; `ok` false right at one.
struct Charges {
    f: Field,
    ok: bool,
};

fn charges(p: vec2<f32>) -> Charges {
    let c = params.scales.x;
    let radiation_only = (params.counts.z & 1u) != 0u;
    let neglected_only = (params.counts.z & 2u) != 0u;
    var f = Field(vec2<f32>(0.0), 0.0);
    for (var k = 0u; k < params.grid.w; k = k + 1u) {
        let head = items[3u * k];
        let end = items[3u * k + 1u];
        let moment = items[3u * k + 2u].x;
        let o = u32(head.x);
        let n = u32(head.y);
        let q = head.z;
        let has_end = head.w > 0.5;
        if (n == 0u) {
            continue;
        }
        let stays = end.w > 0.5;
        // A magnetic moment's dipole field, B_z = -m/r^3 in the plane (radiation.rs
        // `charges`): from the position now for c = inf, the retarded one otherwise.
        if (moment != 0.0 && !radiation_only) {
            var at = vec2<f32>(0.0);
            var there = true;
            if (c == 0.0) {
                if (has_end && end.x <= 0.0) {
                    there = stays;
                    at = end.yz;
                } else {
                    at = present(o, n).x;
                }
            } else {
                if (has_end && c * (-end.x) >= length(p - end.yz)) {
                    there = stays;
                    at = end.yz;
                } else {
                    at = retarded(o, n, p, c).x;
                }
            }
            if (there) {
                let r = length(p - at);
                if (r <= 0.15) {
                    return Charges(f, false);
                }
                f.bz = f.bz - moment / (r * r * r);
            }
        }
        if (q == 0.0) {
            continue;
        }
        if (c == 0.0) {
            // c = ∞: the Coulomb field of the present position; once absorbed, of where
            // it stopped (its charge stays) or none (drained).
            if (has_end && end.x <= 0.0) {
                if (stays) {
                    let d = p - end.yz;
                    let r = length(d);
                    if (r <= 0.15) {
                        return Charges(f, false);
                    }
                    f.e = f.e + d * (q / (r * r * r));
                }
                continue;
            }
            let s = present(o, n);
            let d = p - s.x;
            let r = length(d);
            if (r <= 0.15) {
                return Charges(f, false);
            }
            f.e = f.e + d * (q / (r * r * r));
            continue;
        }
        // Once the light cone of its absorption has passed: at rest where it stopped
        // (its charge stays), or gone (drained).
        let absorbed = has_end && c * (-end.x) >= length(p - end.yz);
        if (absorbed) {
            if (stays && !radiation_only) {
                let d = p - end.yz;
                let r = length(d);
                if (r <= 0.15) {
                    return Charges(f, false);
                }
                f.e = f.e + d * (q / (r * r * r));
            }
        } else {
            let s = retarded(o, n, p, c);
            if (length(p - s.x) <= 0.15) {
                return Charges(f, false);
            }
            let l = lienard(q, c, p, s, radiation_only);
            f.e = f.e + l.e;
            f.bz = f.bz + l.bz;
        }
        if (neglected_only && has_end && end.x <= 0.0) {
            // The dynamics has the absorbed charge at rest at once (or drained).
            if (stays) {
                let d = p - end.yz;
                let r = max(length(d), 1e-6);
                f.e = f.e - d * (q / (r * r * r));
            }
        } else if (neglected_only) {
            let s = present(o, n);
            if (length(p - s.x) <= 0.15) {
                return Charges(f, false);
            }
            let a = accelerated(q, c, p, s);
            f.e = f.e - a.e;
            f.bz = f.bz - a.bz;
        }
    }
    return Charges(f, true);
}

// Antennas: oscillating dipoles p(t) = p0 cos(ωt + φ), exact retarded fields
// (antenna.rs); quasi-static for c = ∞. `ok` false inside an antenna body.
fn antennas(p: vec2<f32>) -> Charges {
    let c = params.scales.x;
    let base = 3u * params.grid.w;
    var f = Field(vec2<f32>(0.0), 0.0);
    for (var k = 0u; k < params.counts.x; k = k + 1u) {
        let a0 = items[base + 2u * k];
        let a1 = items[base + 2u * k + 1u];
        let d = p - a0.xy;
        let r = length(d);
        if (r < a1.z) {
            return Charges(f, false);
        }
        let n = d / r;
        let w = a1.x;
        var ph = a1.y;
        if (c > 0.0) {
            ph = ph - w * r / c;
        }
        let p0 = a0.zw;
        let pm = p0 * cos(ph);
        let pd = p0 * (-w * sin(ph));
        let pdd = p0 * (-w * w * cos(ph));
        var e = (n * (3.0 * dot(n, pm)) - pm) / (r * r * r);
        var bz = 0.0;
        if (c > 0.0) {
            e = e + (n * (3.0 * dot(n, pd)) - pd) / (c * r * r);
            // n × (n × p̈) = n (n·p̈) − p̈
            e = e + (n * dot(n, pdd) - pdd) / (c * c * r);
            bz = (cross_z(pd, n) / (r * r) + cross_z(pdd, n) / (c * r)) / (c * c);
        }
        f.e = f.e + e;
        f.bz = f.bz + bz;
    }
    return Charges(f, true);
}

// Plane waves E = E0 ê cos(ω (t − k̂·x/c) + φ), B = k̂ × E / c.
fn waves(p: vec2<f32>) -> Field {
    let c = params.scales.x;
    let base = 3u * params.grid.w + 2u * params.counts.x;
    var f = Field(vec2<f32>(0.0), 0.0);
    for (var k = 0u; k < params.counts.y; k = k + 1u) {
        let w0 = items[base + 2u * k];
        let w1 = items[base + 2u * k + 1u];
        var ph = w1.z;
        if (c > 0.0) {
            ph = ph - w1.y * dot(w0.xy, p) / c;
        }
        let e = w0.zw * (w1.x * cos(ph));
        f.e = f.e + e;
        if (c > 0.0) {
            f.bz = f.bz + cross_z(w0.xy, e) / c;
        }
    }
    return f;
}

// Static field, bilinear between grid points; `ok` false inside a body.
fn static_field(p: vec2<f32>) -> Charges {
    let w = params.grid.x;
    let h = params.grid.y;
    let u = (p.x - params.area.x) / params.area.z * f32(w) - 0.5;
    let v = (p.y - params.area.y) / params.area.w * f32(h) - 0.5;
    let i0 = u32(clamp(floor(u), 0.0, f32(w - 1u)));
    let j0 = u32(clamp(floor(v), 0.0, f32(h - 1u)));
    let i1 = min(i0 + 1u, w - 1u);
    let j1 = min(j0 + 1u, h - 1u);
    let fu = clamp(u - f32(i0), 0.0, 1.0);
    let fv = clamp(v - f32(j0), 0.0, 1.0);
    let s00 = statics[j0 * w + i0];
    let s10 = statics[j0 * w + i1];
    let s01 = statics[j1 * w + i0];
    let s11 = statics[j1 * w + i1];
    let s = mix(mix(s00, s10, fu), mix(s01, s11, fu), fv);
    return Charges(Field(s.xy, s.z), s.w < 0.5);
}

fn srgb_to_linear(c: vec3<f32>) -> vec3<f32> {
    return pow(c, vec3<f32>(2.2));
}

// Signed compression into [−1, 1] on an asinh scale (radiation.rs `compress`).
fn compress(v: f32, sat: f32, range: f32) -> f32 {
    if (sat <= 0.0) {
        // The quantity vanishes everywhere: nothing to colour.
        return 0.0;
    }
    if ((params.counts.z & 8u) != 0u) {
        // Linear, with a gain.
        return clamp(v * range / sat, -1.0, 1.0);
    }
    let r = sat / range;
    return clamp(asinh(v / r) / asinh(sat / r), -1.0, 1.0);
}

@fragment
fn fragment(in: VertexOutput) -> @location(0) vec4<f32> {
    let p = in.world_position.xy;
    var f = Field(vec2<f32>(0.0), 0.0);
    var ok = true;
    if (params.grid.z == 1u) {
        let s = static_field(p);
        f = s.f;
        ok = ok && s.ok;
    }
    let a = antennas(p);
    f.e = f.e + a.f.e;
    f.bz = f.bz + a.f.bz;
    ok = ok && a.ok;
    let wv = waves(p);
    f.e = f.e + wv.e;
    f.bz = f.bz + wv.bz;
    let ch = charges(p);
    f.e = f.e + ch.f.e;
    f.bz = f.bz + ch.f.bz;
    ok = ok && ch.ok;
    if (!ok) {
        return vec4<f32>(0.0);
    }
    let range = params.scales.w;
    var col: vec3<f32>;
    var s: f32;
    if ((params.counts.z & 4u) != 0u) {
        s = compress(length(f.e), params.scales.z, range);
        col = vec3<f32>(1.0, 0.92, 0.55);
    } else {
        s = compress(f.bz, params.scales.y, range);
        if (s >= 0.0) {
            col = vec3<f32>(1.0, 0.59, 0.16);
        } else {
            col = vec3<f32>(0.16, 0.75, 1.0);
        }
    }
    let alpha = pow(abs(s), 0.8) * 0.9;
    return vec4<f32>(srgb_to_linear(col), alpha);
}
