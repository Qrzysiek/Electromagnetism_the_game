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
    // colour scale, 16 colour |S| (the energy flow), 32 colour of the charges' field alone
    // (the particle-field view); energy flow: part (bits 0-1: total,
    // the charges' own, exchange, the rest's own), 4 averaged over a period, number of
    // frequency groups (bits 8-15).
    counts: vec4<u32>,
    // Energy flow: saturation of |E x B_z| (the flow without its constant c^2/4 pi), 0, 0, 0.
    flow: vec4<f32>,
};

@group(#{MATERIAL_BIND_GROUP}) @binding(0) var<uniform> params: Params;
// World-line samples, two per sample: (τ, x, y, vx), (vy, ax, ay, 0).
@group(#{MATERIAL_BIND_GROUP}) @binding(1) var<storage, read> samples: array<vec4<f32>>;
// Items: charges, seven vec4 each: (offset, count, q, has_end), (τ_end, x_end, y_end,
// 1 if its charge stays where it was absorbed, 0 if drained), (magnetic moment, fade rate
// in a screening cup (0: none), 0, 0), and the quasi-static interaction's continued past
// (its motion in the fields it feels now, field_motion.rs): (U₀, ω²), (Λ U₀, 1 if valid),
// (Λ² U₀, 0), and the jerk of their change along its path (Δȧ, T, 0);
// then antennas, two vec4 each: (x, y, p0x, p0y), (ω, phase now, radius, frequency group
// (255: static)); then waves, two vec4 each: (k̂x, k̂y, êx, êy), (E0, ω, phase now at x = 0,
// frequency group).
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

// The past of a charge now in state `s` continued back by τ (beam::tapered): constant
// acceleration while |a τ| ≤ 0.7 Δv, then fading as sech², so that the velocity changes
// by at most Δv = min(0.1 c, 0.9 (c − |v|)) (beam::taper_reach: never reaching c).
fn tapered(s: State, tau: f32, c: f32) -> State {
    let reach = max(min(0.1 * c, 0.9 * max(c - length(s.v), 0.0)), 1e-20);
    let t_scale = reach / max(length(s.a), 1e-20);
    let u = tau / t_scale;
    if (abs(u) <= 0.7) {
        return State(s.x + s.v * tau + s.a * (0.5 * tau * tau), s.v + s.a * tau, s.a);
    }
    let w = 0.3;
    let z = (abs(u) - 0.7) / w;
    let ln_cosh = z + log(1.0 + exp(-2.0 * z)) - 0.6931472;
    let th = tanh(z);
    return State(
        s.x + s.v * tau + s.a * (t_scale * t_scale * (0.245 + 0.7 * w * z + w * w * ln_cosh)),
        s.v + s.a * (t_scale * sign(u) * (0.7 + w * th)),
        s.a * (1.0 - th * th),
    );
}

// The quasi-static beam interaction's field of a charge now in state `s`: its past
// continued back with `tapered` (beam::accelerated_fields).
fn accelerated(q: f32, c: f32, p: vec2<f32>, s: State) -> Field {
    // Newton on the strictly decreasing g(τ) = −c τ − |p − x(τ)| (g(0) ≤ 0), kept inside
    // the bracket of the signs seen (bisection where it would leave it).
    var tau = -length(p - s.x) / c;
    var lo = -1e30;
    var hi = 0.0;
    for (var k = 0; k < 24; k = k + 1) {
        let st = tapered(s, tau, c);
        let d = p - st.x;
        let dist = max(length(d), 1e-6);
        let g = -c * tau - dist;
        if (g > 0.0) {
            lo = max(lo, tau);
        } else {
            hi = min(hi, tau);
        }
        let slope = -c + dot(d, st.v) / dist;
        var next = tau - g / slope;
        if (!(next < hi && next > lo)) {
            if (lo > -1e29) {
                next = 0.5 * (lo + hi);
            } else {
                next = hi - 2.0 * max(hi - tau, dist / c);
            }
        }
        tau = next;
    }
    return lienard(q, c, p, tapered(s, tau, c), false);
}

// (C₀, S₁, C₁, S₂) of the motion in uniform fields at proper time s (field_motion.rs): the
// series where ω²s² is small (here below 0.5: f32 would lose the closed forms' small
// differences), the closed forms otherwise.
fn fm_functions(w2: f32, s: f32) -> vec4<f32> {
    let z = w2 * s * s;
    if (abs(z) < 0.5) {
        let s1 = s * (1.0 - z / 6.0 * (1.0 - z / 20.0 * (1.0 - z / 42.0 * (1.0 - z / 72.0))));
        let c1 = s * s * 0.5 * (1.0 - z / 12.0 * (1.0 - z / 30.0 * (1.0 - z / 56.0 * (1.0 - z / 90.0))));
        let s2 = s * s * s / 6.0 * (1.0 - z / 20.0 * (1.0 - z / 42.0 * (1.0 - z / 72.0 * (1.0 - z / 110.0))));
        return vec4<f32>(1.0 - w2 * c1, s1, c1, s2);
    }
    if (w2 > 0.0) {
        let w = sqrt(w2);
        let h = sin(0.5 * w * s);
        return vec4<f32>(cos(w * s), sin(w * s) / w, 2.0 * h * h / w2, (s - sin(w * s) / w) / w2);
    }
    let k = sqrt(-w2);
    let h = sinh(0.5 * k * s);
    return vec4<f32>(cosh(k * s), sinh(k * s) / k, 2.0 * h * h / (-w2), (sinh(k * s) / k - s) / (-w2));
}

struct FmPoint {
    st: State,
    // c t(s) relative to now, and the excursion |U(s) − U₀|/c.
    ct: f32,
    excursion: f32,
    slope_u: vec3<f32>,
};

// The jerk correction of charge k's continued past at the coordinate time tau from now
// (field_motion.rs JerkCorrection): δa = Δȧ T g(τ/T), g(u) = u (1 − u²)³ inside the taper
// and 0 beyond, with δv and δx its integrals (written in τ inside, so that f32 keeps the
// small times).
struct Corr {
    dx: vec2<f32>,
    dv: vec2<f32>,
    da: vec2<f32>,
};

fn jerk_corr(k: u32, tau: f32) -> Corr {
    let p3 = items[7u * k + 6u];
    let jerk = p3.xy;
    let big_t = max(p3.z, 1e-30);
    let u = tau / big_t;
    if (abs(u) < 1.0) {
        let u2 = u * u;
        let w = 1.0 - u2;
        return Corr(
            jerk * (tau * tau * tau * (1.0 / 6.0 - u2 * (0.15 - u2 * (1.0 / 14.0 - u2 / 72.0)))),
            jerk * (tau * tau * (0.5 - u2 * (0.75 - u2 * (0.5 - u2 / 8.0)))),
            jerk * (tau * w * w * w),
        );
    }
    let su = sign(u);
    return Corr(
        jerk * (su * big_t * big_t * big_t * (187.0 / 2520.0) + (tau - su * big_t) * big_t * big_t / 8.0),
        jerk * (big_t * big_t / 8.0),
        vec2<f32>(0.0),
    );
}

// The continued past of charge k (now in state `now`) at proper time s, with its jerk
// correction; slope_u: (U⁰, dX/ds).
fn fm_at(k: u32, now: State, s: f32, c: f32) -> FmPoint {
    let p0 = items[7u * k + 3u];
    let p1 = items[7u * k + 4u];
    let p2 = items[7u * k + 5u];
    let f = fm_functions(p0.w, s);
    let u = p0.xyz + f.y * p1.xyz + f.z * p2.xyz;
    let du = f.x * p1.xyz + f.y * p2.xyz;
    let dx = s * p0.xyz + f.z * p1.xyz + f.w * p2.xyz;
    let v = u.yz * (c / u.x);
    let a = (du.yz * u.x - u.yz * du.x) * (c * c / (u.x * u.x * u.x));
    let corr = jerk_corr(k, dx.x / c);
    let slope = vec3<f32>(u.x, u.yz + corr.dv * (u.x / c));
    return FmPoint(
        State(now.x + dx.yz + corr.dx, v + corr.dv, a + corr.da),
        dx.x,
        length(u - p0.xyz) / c,
        slope,
    );
}

// The quasi-static interaction's field at p of charge k, now in state `now`: its past
// continued along its motion in the fields it feels now with the jerk of their change,
// blending into `accelerated` where that swings far (excursion 0.5 to 1, quintic
// smoothstep; beam::Continuation::fields), or `accelerated` alone where no motion was
// given or its retarded point is not found.
fn continued(q: f32, c: f32, p: vec2<f32>, now: State, k: u32) -> Field {
    if (items[7u * k + 4u].w < 0.5) {
        return accelerated(q, c, p, now);
    }
    // Newton on G(s) = −c t(s) − |p − x(s)| from the uniform-motion guess.
    let u0 = items[7u * k + 3u].xyz;
    let gamma0 = u0.x / c;
    let v0 = u0.yz / gamma0;
    let r = p - now.x;
    let rv = dot(r, v0);
    let v2 = dot(v0, v0);
    var s = -(rv + sqrt(rv * rv + (c * c - v2) * dot(r, r))) / (c * c - v2) / gamma0;
    var pt = fm_at(k, now, s, c);
    var converged = false;
    for (var i = 0; i < 12; i = i + 1) {
        let d = p - pt.st.x;
        let dist = max(length(d), 1e-6);
        let g = -pt.ct - dist;
        let slope = -pt.slope_u.x + dot(d, pt.slope_u.yz) / dist;
        let next = s - g / slope;
        if (abs(next - s) <= 1e-6 * (abs(s) + 1e-3)) {
            converged = true;
        }
        s = next;
        pt = fm_at(k, now, s, c);
        if (converged) {
            break;
        }
    }
    var w = 0.0;
    if (converged) {
        let u = clamp((pt.excursion - 0.5) / 0.5, 0.0, 1.0);
        w = 1.0 - u * u * u * (10.0 - 15.0 * u + 6.0 * u * u);
    }
    var f = Field(vec2<f32>(0.0), 0.0);
    if (w > 0.0) {
        let l = lienard(q, c, p, pt.st, false);
        f = Field(l.e * w, l.bz * w);
    }
    if (w < 1.0) {
        let a = accelerated(q, c, p, now);
        f = Field(f.e + a.e * (1.0 - w), f.bz + a.bz * (1.0 - w));
    }
    return f;
}

// The quasi-static interaction's field of a charge fading in a screening cup, now in the
// uniform state `s` (beam::Fade::quasi_static_fields): the field of its uniform motion
// (Heaviside, Jackson §11.10; B = v × E / c²) with the charge at the retarded time of
// that motion, faded since it entered the cup at τ_end.
fn fading_quasi_static(q: f32, rate: f32, tau_end: f32, c: f32, p: vec2<f32>, s: State) -> Field {
    let d = p - s.x;
    let dv = dot(d, s.v);
    let v2 = dot(s.v, s.v);
    let tau = (dv + sqrt(dv * dv + (c * c - v2) * dot(d, d))) / (c * c - v2);
    var qe = q;
    if (-tau > tau_end) {
        qe = q * exp(-rate * (-tau - tau_end));
    }
    let dist = max(length(d), 1e-6);
    let beta = s.v / c;
    let b2 = dot(beta, beta);
    let bn = dot(beta, d / dist);
    let e = d * (qe * (1.0 - b2) / (dist * dist * dist * pow(1.0 - b2 + bn * bn, 1.5)));
    return Field(e, cross_z(s.v, e) / (c * c));
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
        let head = items[7u * k];
        let end = items[7u * k + 1u];
        let moment = items[7u * k + 2u].x;
        let fade = items[7u * k + 2u].y;
        let o = u32(head.x);
        let n = u32(head.y);
        let q = head.z;
        let has_end = head.w > 0.5;
        if (n == 0u) {
            continue;
        }
        let stays = end.w > 0.5;
        // Flying on into a screening cup, its charge fading (radiation.rs `charges`).
        let fades = has_end && !stays && fade > 0.0;
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
            // it stopped (its charge stays), of its fading charge flying on into a
            // screening cup (factor e^{−rate (t − t_end)}, τ_end = end.x ≤ 0), or none
            // (drained).
            if (has_end && end.x <= 0.0) {
                if (stays) {
                    let d = p - end.yz;
                    let r = length(d);
                    if (r <= 0.15) {
                        return Charges(f, false);
                    }
                    f.e = f.e + d * (q / (r * r * r));
                } else if (fades) {
                    // Screened inside the cup: not masked like a free charge.
                    let d = p - present(o, n).x;
                    let r = max(length(d), 1e-6);
                    f.e = f.e + d * (q * exp(fade * end.x) / (r * r * r));
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
        // (its charge stays), or gone (drained). Flying on into a screening cup it keeps
        // its world line, its charge fading at the retarded time.
        let absorbed = has_end && !fades && c * (-end.x) >= length(p - end.yz);
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
            var qk = q;
            // Seen inside the cup (screened): not masked like a free charge.
            var in_cup = false;
            if (fades) {
                // The retarded time, relative to now: |p − x(τ)| = −c τ.
                let tau_r = -length(p - s.x) / c;
                if (tau_r > end.x) {
                    qk = q * exp(-fade * (tau_r - end.x));
                    in_cup = true;
                }
            }
            if (length(p - s.x) <= 0.15 && !in_cup) {
                return Charges(f, false);
            }
            let l = lienard(qk, c, p, s, radiation_only);
            f.e = f.e + l.e;
            f.bz = f.bz + l.bz;
        }
        if (neglected_only && has_end && end.x <= 0.0) {
            // The dynamics has the absorbed charge at rest at once, fading in its cup as it
            // moves on uniformly, or drained.
            if (stays) {
                let d = p - end.yz;
                let r = max(length(d), 1e-6);
                f.e = f.e - d * (q / (r * r * r));
            } else if (fades) {
                let a = fading_quasi_static(q, fade, end.x, c, p, present(o, n));
                f.e = f.e - a.e;
                f.bz = f.bz - a.bz;
            }
        } else if (neglected_only) {
            let s = present(o, n);
            if (length(p - s.x) <= 0.15) {
                return Charges(f, false);
            }
            let a = continued(q, c, p, s, k);
            f.e = f.e - a.e;
            f.bz = f.bz - a.bz;
        }
    }
    return Charges(f, true);
}

// Antennas: oscillating dipoles p(t) = p0 cos(ωt + φ), exact retarded fields
// (antenna.rs); quasi-static for c = ∞. `ok` false inside an antenna body.
fn antennas(p: vec2<f32>) -> Charges {
    return antennas_at(p, 0.0, -1);
}

// The antennas of frequency group `group` (all for −1), their phases advanced by `dph`
// (π/2: a quarter period later, for the average over a period).
fn antennas_at(p: vec2<f32>, dph: f32, group: i32) -> Charges {
    let c = params.scales.x;
    let base = 7u * params.grid.w;
    var f = Field(vec2<f32>(0.0), 0.0);
    for (var k = 0u; k < params.counts.x; k = k + 1u) {
        let a0 = items[base + 2u * k];
        let a1 = items[base + 2u * k + 1u];
        let d = p - a0.xy;
        let r = length(d);
        if (r < a1.z) {
            return Charges(f, false);
        }
        if (group >= 0 && i32(a1.w + 0.5) != group) {
            continue;
        }
        let n = d / r;
        let w = a1.x;
        var ph = a1.y + dph;
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
    return waves_at(p, 0.0, -1);
}

// The waves of frequency group `group` (all for −1), their phases advanced by `dph`.
fn waves_at(p: vec2<f32>, dph: f32, group: i32) -> Field {
    let c = params.scales.x;
    let base = 7u * params.grid.w + 2u * params.counts.x;
    var f = Field(vec2<f32>(0.0), 0.0);
    for (var k = 0u; k < params.counts.y; k = k + 1u) {
        let w0 = items[base + 2u * k];
        let w1 = items[base + 2u * k + 1u];
        if (group >= 0 && i32(w1.w + 0.5) != group) {
            continue;
        }
        var ph = w1.z + dph;
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

// E × (b ẑ) in the plane: the energy flow without its constant c²/4π.
fn flow_of(e: vec2<f32>, b: f32) -> vec2<f32> {
    return vec2<f32>(e.y * b, -e.x * b);
}

// The energy flow of the antennas and waves averaged over a period (without c²/4π): per
// frequency group ½ (S + S⊥), S⊥ a quarter period later (exact for one frequency; the cross
// terms of different frequencies average out); static antennas (group 255) as they are.
fn averaged_oscillators(p: vec2<f32>) -> vec2<f32> {
    let n_groups = (params.counts.w >> 8u) & 255u;
    var s = vec2<f32>(0.0);
    for (var g = 0u; g < n_groups; g = g + 1u) {
        let a0 = antennas_at(p, 0.0, i32(g)).f;
        let w0 = waves_at(p, 0.0, i32(g));
        let a1 = antennas_at(p, 1.5707963, i32(g)).f;
        let w1 = waves_at(p, 1.5707963, i32(g));
        s = s + 0.5 * (flow_of(a0.e + w0.e, a0.bz + w0.bz) + flow_of(a1.e + w1.e, a1.bz + w1.bz));
    }
    return s;
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
    var st = Field(vec2<f32>(0.0), 0.0);
    if (params.grid.z == 1u) {
        let s = static_field(p);
        st = s.f;
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
    // What the view shows besides the moving charges.
    let rest = f;
    let ch = charges(p);
    f.e = f.e + ch.f.e;
    f.bz = f.bz + ch.f.bz;
    ok = ok && ch.ok;
    if (!ok) {
        return vec4<f32>(0.0);
    }
    if ((params.counts.z & 32u) != 0u) {
        // The particle-field view: the charges' field alone (the rest enters only the
        // energy flow's exchange and external parts below).
        f = ch.f;
    }
    let range = params.scales.w;
    var col: vec3<f32>;
    var s: f32;
    if ((params.counts.z & 16u) != 0u) {
        // The energy flow of the part chosen (radiation.rs `flows_at`), without c²/4π:
        // the charges' own terms, the exchange terms, the rest's own terms, or all.
        let part = params.counts.w & 3u;
        var s_rest = flow_of(rest.e, rest.bz);
        var s_x = flow_of(rest.e, ch.f.bz) + flow_of(ch.f.e, rest.bz);
        if ((params.counts.w & 4u) != 0u) {
            // Averaged over a period, the charges held still: the rest's static part plus
            // each frequency's average; the charges exchange with the static part only.
            var still = st;
            for (var k = 0u; k < params.counts.x; k = k + 1u) {
                if (items[7u * params.grid.w + 2u * k + 1u].w > 254.5) {
                    let sa = antennas_at(p, 0.0, 255).f;
                    still.e = still.e + sa.e;
                    still.bz = still.bz + sa.bz;
                    break;
                }
            }
            s_rest = flow_of(still.e, still.bz) + averaged_oscillators(p);
            s_x = flow_of(still.e, ch.f.bz) + flow_of(ch.f.e, still.bz);
        }
        let s_own = flow_of(ch.f.e, ch.f.bz);
        var flow = s_own + s_x + s_rest;
        if (part == 1u) {
            flow = s_own;
        } else if (part == 2u) {
            flow = s_x;
        } else if (part == 3u) {
            flow = s_rest;
        }
        s = compress(length(flow), params.flow.x, range);
        col = vec3<f32>(0.55, 1.0, 0.62);
    } else if ((params.counts.z & 4u) != 0u) {
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
