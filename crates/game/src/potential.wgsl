// Field maps, computed per pixel (visual only, f32).
//
// Mode 0, potential: colour shows the particle's potential energy relative to the launch
// point, U / T0 = [q (phi - phi_A) - m (B_z - B_z,A)] / T0 (m: its magnetic moment), in
// the static part of the flight's field: red uphill, blue downhill, contours every T0/4.
// Dark: forbidden by energy conservation for every particle shown (each group, a shot or
// a beam, forbids U > E with E its highest total energy); its boundary, the turning line,
// is bright.
//
// Mode 1, magnetic field: B_z in units of b_ref: for a charged particle the field in which
// it circles with a 5-cell gyroradius (so the value is 5 / r_gyro); for a neutral one with
// a moment m, the field with |m B_z| = T0. Orange: B out of the plane, teal: into the
// plane; contours every 0.25.
//
// Lines are antialiased with screen-space derivatives, so the maps are sharp at any zoom.

#import bevy_sprite::mesh2d_vertex_output::VertexOutput

struct Params {
    // Energy limits, one per group of particles: x: q / T, y: -m b_ref / T, z: E / T.
    // Forbidden for the group where x Phi + y B_z / b_ref > z.
    limits: array<vec4<f32>, 64>,
    // Numbers of solid charges, charge clouds, induced charges, magnets.
    counts: vec4<u32>,
    // Numbers of circular coils, straight segments, static antennas, electrode panels.
    counts2: vec4<u32>,
    // Uniform stray field: E_x, E_y (potential -E.(x - origin)), B_z / b_ref.
    uniform_field: vec4<f32>,
    // xy: origin of the uniform field's potential, z: antenna body radius, w: coil wire
    // radius.
    origin: vec4<f32>,
    // Number of energy limits (0: no dark region).
    limit_count: u32,
    // q / T0 of the shot the colours show: weight of Phi in U / T0.
    phi_weight: f32,
    // U / T0 at the launch point.
    u_a: f32,
    // 0: potential, 1: magnetic field.
    mode: u32,
    // -m b_ref / T0: weight of B_z / b_ref in U / T0 (0 without a magnetic moment).
    moment_weight: f32,
    // Tube levels (z-invariant): numbers of charged segments (prism cross-sections) and
    // of line charges, after the electrode panels.
    line_segments: u32,
    line_charges: u32,
};

@group(#{MATERIAL_BIND_GROUP}) @binding(0) var<uniform> params: Params;
// The sources, one category after the other (potential.rs, `Items`): solid charges
// (x, y, Q, radius), clouds (x, y, Q, R), induced charges (x, y, Q, z), magnets
// (x, y, mu / b_ref, radius), coils (x, y, radius, kappa / b_ref), segments (two each:
// (ax, ay, bx, by), (kappa / b_ref, 0, 0, 0)), static antennas (x, y, px, py), electrode
// panels (three each: (a, sigma), (b, size), (c, area)), charged segments of a
// z-invariant level (two each: (ax, ay, bx, by), (sigma, 0, 0, 0)), line charges
// (x, y, lambda, epsilon^2: softening).
@group(#{MATERIAL_BIND_GROUP}) @binding(1) var<storage, read> items: array<vec4<f32>>;

fn first_cloud() -> u32 {
    return params.counts.x;
}
fn first_induced() -> u32 {
    return first_cloud() + params.counts.y;
}
fn first_magnet() -> u32 {
    return first_induced() + params.counts.z;
}
fn first_loop() -> u32 {
    return first_magnet() + params.counts.w;
}
fn first_segment() -> u32 {
    return first_loop() + params.counts2.x;
}
fn first_antenna() -> u32 {
    return first_segment() + 2u * params.counts2.y;
}
fn first_panel() -> u32 {
    return first_antenna() + params.counts2.z;
}
fn first_line_segment() -> u32 {
    return first_panel() + 3u * params.counts2.w;
}
fn first_line_charge() -> u32 {
    return first_line_segment() + 2u * params.line_segments;
}

// Potential of a segment of unit surface density in the z-invariant world (as
// physics::zinv::segment_integrals): -[G(u) - G(u - L)], G(w) = w ln(w^2 + v^2) - 2w
// + 2v atan(w / v).
fn line_segment_potential(p: vec2<f32>, a: vec2<f32>, b: vec2<f32>) -> f32 {
    let d = b - a;
    let l = length(d);
    let t = d / l;
    let r = p - a;
    let u = dot(r, t);
    let v = r.x * -t.y + r.y * t.x;
    let w0 = u;
    let w1 = u - l;
    let s0 = max(w0 * w0 + v * v, 1e-12);
    let s1 = max(w1 * w1 + v * v, 1e-12);
    var g0 = w0 * log(s0) - 2.0 * w0;
    var g1 = w1 * log(s1) - 2.0 * w1;
    if (abs(v) > 1e-9) {
        g0 = g0 + 2.0 * v * atan(w0 / v);
        g1 = g1 + 2.0 * v * atan(w1 / v);
    }
    return -(g0 - g1);
}

fn srgb_to_linear(c: vec3<f32>) -> vec3<f32> {
    return pow(c, vec3<f32>(2.2));
}

// Coverage of a line of `width_px` pixels at distance `d` (in the units of `fw`).
fn line_cover(d: f32, fw: f32, width_px: f32) -> f32 {
    return 1.0 - smoothstep(0.5 * width_px * fw, (0.5 * width_px + 1.0) * fw, d);
}

// Complete elliptic integrals K(m), E(m) by the arithmetic-geometric mean.
fn elliptic_ke(m: f32, m1: f32) -> vec2<f32> {
    var a = 1.0;
    var b = sqrt(m1);
    var c = sqrt(m);
    var sum = 0.5 * c * c;
    var p = 0.5;
    for (var i = 0; i < 12; i = i + 1) {
        let an = 0.5 * (a + b);
        let bn = sqrt(a * b);
        c = 0.5 * (a - b);
        a = an;
        b = bn;
        p = p * 2.0;
        sum = sum + p * c * c;
    }
    let k = 3.14159265 / (2.0 * a);
    return vec2<f32>(k, k * (1.0 - sum));
}

// Signed distance to the nearest solid (charge, magnet, static antenna, coil wire).
fn solid_distance(p: vec2<f32>) -> f32 {
    var d = 1e9;
    for (var i = 0u; i < params.counts.x; i = i + 1u) {
        let c = items[i];
        d = min(d, distance(p, c.xy) - c.w);
    }
    for (var i = 0u; i < params.counts.w; i = i + 1u) {
        let m = items[first_magnet() + i];
        d = min(d, distance(p, m.xy) - m.w);
    }
    for (var i = 0u; i < params.counts2.z; i = i + 1u) {
        let a = items[first_antenna() + i];
        d = min(d, distance(p, a.xy) - params.origin.z);
    }
    for (var i = 0u; i < params.counts2.x; i = i + 1u) {
        let l = items[first_loop() + i];
        d = min(d, abs(distance(p, l.xy) - l.z) - params.origin.w);
    }
    for (var i = 0u; i < params.counts2.y; i = i + 1u) {
        let s = items[first_segment() + 2u * i];
        let ab = s.zw - s.xy;
        let t = clamp(dot(p - s.xy, ab) / dot(ab, ab), 0.0, 1.0);
        d = min(d, distance(p, s.xy + ab * t) - params.origin.w);
    }
    return d;
}

// Potential of a uniformly charged triangle (unit density) at r (Wilton et al. 1984,
// as physics::panel): sum over the edges of P0 f - |h| beta, with f written without
// cancellation for f32 (from the identity (R+ + l+)(R+ - l+) = R0^2 = (R- + l-)(R- - l-)).
fn panel_potential(r: vec3<f32>, a: vec3<f32>, b: vec3<f32>, c: vec3<f32>) -> f32 {
    let n = normalize(cross(b - a, c - a));
    let h = dot(n, r - a);
    let abs_h = abs(h);
    let rho = r - n * h;
    var pot = 0.0;
    var beta_sum = 0.0;
    for (var e = 0; e < 3; e = e + 1) {
        var p = a;
        var q = b;
        if (e == 1) {
            p = b;
            q = c;
        } else if (e == 2) {
            p = c;
            q = a;
        }
        let l = normalize(q - p);
        let m = cross(l, n);
        let p0 = dot(p - rho, m);
        let l_minus = dot(p - rho, l);
        let l_plus = dot(q - rho, l);
        let r_minus = length(r - p);
        let r_plus = length(r - q);
        let r0_sq = p0 * p0 + h * h;
        // On the edge's line (R0 = 0) the term P0 f vanishes.
        if (r0_sq > 1e-12 * (r_minus * r_minus + r_plus * r_plus)) {
            var f = 0.0;
            if (l_minus >= 0.0) {
                f = log((r_plus + l_plus) / (r_minus + l_minus));
            } else if (l_plus <= 0.0) {
                f = log((r_minus - l_minus) / (r_plus - l_plus));
            } else {
                f = log((r_plus + l_plus) * (r_minus - l_minus) / r0_sq);
            }
            pot = pot + p0 * f;
            beta_sum = beta_sum + atan(p0 * l_plus / (r0_sq + abs_h * r_plus))
                - atan(p0 * l_minus / (r0_sq + abs_h * r_minus));
        }
    }
    return pot - abs_h * beta_sum;
}

// The electric potential of the static sources at p (z = 0).
fn potential(p: vec2<f32>) -> f32 {
    var phi = 0.0;
    // Solid charges (in the plane), clouds, induced charges (off the plane).
    for (var i = 0u; i < params.counts.x; i = i + 1u) {
        let c = items[i];
        phi = phi + c.z / max(distance(p, c.xy), 1e-4);
    }
    for (var i = 0u; i < params.counts.y; i = i + 1u) {
        let c = items[first_cloud() + i];
        let r = distance(p, c.xy);
        if (r < c.w) {
            // Inside the uniform sphere.
            phi = phi + c.z * (3.0 * c.w * c.w - r * r) / (2.0 * c.w * c.w * c.w);
        } else {
            phi = phi + c.z / r;
        }
    }
    for (var i = 0u; i < params.counts.z; i = i + 1u) {
        let c = items[first_induced() + i];
        let d = p - c.xy;
        phi = phi + c.z / max(sqrt(dot(d, d) + c.w * c.w), 1e-4);
    }
    // Static antennas: n.p / r^2.
    for (var i = 0u; i < params.counts2.z; i = i + 1u) {
        let a = items[first_antenna() + i];
        let d = p - a.xy;
        let r = max(length(d), 1e-4);
        phi = phi + dot(d, a.zw) / (r * r * r);
    }
    // Electrode panels and their mirror images (equal in the plane): exactly within two
    // panel sizes, by the three-point rule beyond.
    let r3 = vec3<f32>(p, 0.0);
    for (var i = 0u; i < params.counts2.w; i = i + 1u) {
        let k = first_panel() + 3u * i;
        let ta = items[k];
        let tb = items[k + 1u];
        let tc = items[k + 2u];
        let centroid = (ta.xyz + tb.xyz + tc.xyz) / 3.0;
        if (distance(r3, centroid) < 2.0 * tb.w) {
            phi = phi + 2.0 * ta.w * panel_potential(r3, ta.xyz, tb.xyz, tc.xyz);
        } else {
            let y1 = (4.0 * ta.xyz + tb.xyz + tc.xyz) / 6.0;
            let y2 = (ta.xyz + 4.0 * tb.xyz + tc.xyz) / 6.0;
            let y3 = (ta.xyz + tb.xyz + 4.0 * tc.xyz) / 6.0;
            let s = 1.0 / distance(r3, y1) + 1.0 / distance(r3, y2) + 1.0 / distance(r3, y3);
            phi = phi + 2.0 * ta.w * tc.w * s / 3.0;
        }
    }
    // Tube levels: charged segments and line charges (potential -2 lambda ln r).
    for (var i = 0u; i < params.line_segments; i = i + 1u) {
        let k = first_line_segment() + 2u * i;
        let s = items[k];
        phi = phi + items[k + 1u].x * line_segment_potential(p, s.xy, s.zw);
    }
    for (var i = 0u; i < params.line_charges; i = i + 1u) {
        let c = items[first_line_charge() + i];
        let d = p - c.xy;
        // Softened by w = epsilon^2 (the map shows the space charge, not each
        // macroparticle).
        phi = phi - c.z * log(dot(d, d) + c.w);
    }
    // The uniform stray field.
    return phi - dot(params.uniform_field.xy, p - params.origin.xy);
}

fn potential_colour(p: vec2<f32>) -> vec3<f32> {
    let phi = potential(p);
    var b = 0.0;
    if (params.moment_weight != 0.0 || params.limit_count > 0u) {
        b = magnetic_field(p);
    }
    // (Zero weights are skipped: 0 times a huge field would be NaN.)
    var u = params.phi_weight * phi - params.u_a;
    if (params.moment_weight != 0.0) {
        u = u + params.moment_weight * b;
    }
    // Excess of U over the allowed energy, the least over the groups: > 0 forbidden.
    var g = 0.0;
    if (params.limit_count > 0u) {
        g = 1e30;
        for (var i = 0u; i < params.limit_count; i = i + 1u) {
            let l = params.limits[i];
            var excess = l.x * phi - l.z;
            if (l.y != 0.0) {
                excess = excess + l.y * b;
            }
            g = min(g, excess);
        }
    } else {
        g = -1e30;
    }

    let s = tanh(u / 1.5);
    let base = vec3<f32>(0.10, 0.11, 0.14);
    var col: vec3<f32>;
    if (s >= 0.0) {
        col = base + s * vec3<f32>(0.55, 0.10, 0.02);
    } else {
        col = base - s * vec3<f32>(0.02, 0.18, 0.55);
    }
    // Forbidden region, antialiased edge.
    let fw_g = max(fwidth(g), 1e-6);
    let forbidden = smoothstep(-0.5 * fw_g, 0.5 * fw_g, g);
    col = mix(col, col * 0.3, forbidden);
    // Contours every T0/4, faded out smoothly where they crowd closer than a few pixels
    // (a hard cut-off would leave a visible edge, and dense lines alias into moiré).
    let k = u * 4.0;
    let fw_k = fwidth(k);
    let fade = 1.0 - smoothstep(0.06, 0.25, fw_k);
    if (fade > 0.0) {
        let f = fract(k);
        let d = min(f, 1.0 - f);
        let cover = line_cover(d, fw_k, 1.0) * (1.0 - forbidden) * fade;
        col = mix(col, col + vec3<f32>(0.10), cover);
    }
    // Turning line (g = 0).
    let turning = line_cover(abs(g), fw_g, 2.0);
    return mix(col, vec3<f32>(0.95, 0.95, 0.85), turning);
}

fn magnetic_field(p: vec2<f32>) -> f32 {
    // The uniform stray field.
    var b = params.uniform_field.z;
    // Dipoles with moment along z: B_z = -mu / r^3 in the plane.
    for (var i = 0u; i < params.counts.w; i = i + 1u) {
        let m = items[first_magnet() + i];
        let r = max(distance(p, m.xy), 1e-3);
        b = b - m.z / (r * r * r);
    }
    // Circular coils in the plane (z = 0): B_z = 2k / (alpha^2 beta) [(a^2 - rho^2) E + alpha^2 K].
    for (var i = 0u; i < params.counts2.x; i = i + 1u) {
        let l = items[first_loop() + i];
        let a = l.z;
        let rho = distance(p, l.xy);
        let alpha2 = max((a - rho) * (a - rho), 1e-8);
        let beta2 = (a + rho) * (a + rho);
        let ke = elliptic_ke(4.0 * a * rho / beta2, alpha2 / beta2);
        b = b + 2.0 * l.w / (alpha2 * sqrt(beta2)) * ((a * a - rho * rho) * ke.y + alpha2 * ke.x);
    }
    // Straight segments.
    for (var i = 0u; i < params.counts2.y; i = i + 1u) {
        let s = items[first_segment() + 2u * i];
        let kappa = items[first_segment() + 2u * i + 1u].x;
        let ra = s.xy - p;
        let rb = s.zw - p;
        let la = length(ra);
        let lb = length(rb);
        let cross_z = ra.x * rb.y - ra.y * rb.x;
        let denom = max(la * lb * (la * lb + dot(ra, rb)), 1e-12);
        b = b + kappa * cross_z * (la + lb) / denom;
    }
    return b;
}

// Colour of the magnetic map for field value b; `fw_k` is the screen-space derivative
// of the contour variable 4b, computed by the caller (no derivatives in here, so it can
// be used for supersampling in non-uniform control flow).
fn magnetic_colour_of(b: f32, fw_k: f32) -> vec3<f32> {
    let s = tanh(b / 1.5);
    let base = vec3<f32>(0.10, 0.11, 0.14);
    var col: vec3<f32>;
    if (s >= 0.0) {
        col = base + s * vec3<f32>(0.60, 0.35, 0.02);
    } else {
        col = base - s * vec3<f32>(0.02, 0.45, 0.42);
    }
    let k = b * 4.0;
    let fade = 1.0 - smoothstep(0.06, 0.25, fw_k);
    if (fade > 0.0) {
        let f = fract(k);
        let d = min(f, 1.0 - f);
        col = mix(col, col + vec3<f32>(0.10), line_cover(d, fw_k, 1.0) * fade);
    }
    return col;
}

// Magnetic map with adaptive supersampling: where the colour changes quickly between
// neighbouring pixels (sign changes next to wires and magnets), the field is evaluated
// at 4 × 4 points inside the pixel and the colours are averaged.
fn magnetic_colour(p: vec2<f32>) -> vec3<f32> {
    let b = magnetic_field(p);
    let fw_k = fwidth(b * 4.0);
    let fw_s = fwidth(tanh(b / 1.5));
    let dx = dpdx(p);
    let dy = dpdy(p);
    if (fw_s < 0.05) {
        return magnetic_colour_of(b, fw_k);
    }
    var acc = vec3<f32>(0.0);
    for (var i = 0; i < 4; i = i + 1) {
        for (var j = 0; j < 4; j = j + 1) {
            let o = (f32(i) + 0.5) / 4.0 - 0.5;
            let q = (f32(j) + 0.5) / 4.0 - 0.5;
            acc = acc + magnetic_colour_of(magnetic_field(p + o * dx + q * dy), fw_k);
        }
    }
    return acc / 16.0;
}

@fragment
fn fragment(in: VertexOutput) -> @location(0) vec4<f32> {
    let p = in.world_position.xy;
    var col: vec3<f32>;
    if (params.mode == 1u) {
        col = magnetic_colour(p);
    } else {
        col = potential_colour(p);
    }
    // Inside a solid.
    let inside = solid_distance(p);
    let fw_p = max(fwidth(inside), 1e-6);
    let solid = 1.0 - smoothstep(-0.5 * fw_p, 0.5 * fw_p, inside);
    col = mix(col, vec3<f32>(0.08, 0.08, 0.09), solid);
    return vec4<f32>(srgb_to_linear(clamp(col, vec3<f32>(0.0), vec3<f32>(1.0))), 1.0);
}
