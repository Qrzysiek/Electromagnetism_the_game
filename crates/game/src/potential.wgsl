// Field maps, computed per pixel (visual only, f32).
//
// Mode 0, potential: colour shows the particle's potential energy relative to the launch
// point, U / T0 = [q (phi - phi_A) - m (B_z - B_z,A)] / T0 (m: its magnetic moment): red
// uphill, blue downhill, contours every T0/4. U > T0 is forbidden by energy conservation
// (dark); its boundary, the turning line U = T0, is bright.
//
// Mode 1, magnetic field: B_z in units of b_ref: for a charged particle the field in which
// it circles with a 5-cell gyroradius (so the value is 5 / r_gyro); for a neutral one with
// a moment m, the field with |m B_z| = T0. Orange: B out of the plane, teal: into the
// plane; contours every 0.25.
//
// Lines are antialiased with screen-space derivatives, so the maps are sharp at any zoom.

#import bevy_sprite::mesh2d_vertex_output::VertexOutput

struct Params {
    // xy: position (cells), z: w = q Q / T0, w: sphere radius (cells) for the first
    // `solid` charges (fixed charges, in the plane), else the charge's z (induced charges
    // of metal: images and equivalent charges off the plane).
    charges: array<vec4<f32>, 1024>,
    // xy: position, z: mu / b_ref, w: sphere radius.
    magnets: array<vec4<f32>, 64>,
    // xy: centre, z: radius, w: kappa / b_ref.
    loops: array<vec4<f32>, 16>,
    // Straight coil segments: (ax, ay, bx, by), current from a to b.
    segments: array<vec4<f32>, 64>,
    // x: kappa / b_ref of the segment.
    segment_kappa: array<vec4<f32>, 64>,
    // Numbers of charges, magnets, loops, segments.
    counts: vec4<u32>,
    // U / T0 at the launch point.
    u_a: f32,
    // Number of leading charges that are solid spheres in the plane.
    solid: u32,
    // 0: potential, 1: magnetic field.
    mode: u32,
    // Coil wire radius (cells).
    wire: f32,
    // -m b_ref / T0: weight of B_z / b_ref in U / T0 (0 without a magnetic moment).
    moment_weight: f32,
};

@group(#{MATERIAL_BIND_GROUP}) @binding(0) var<uniform> params: Params;

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

// Signed distance to the nearest solid (charge, magnet, coil wire).
fn solid_distance(p: vec2<f32>) -> f32 {
    var d = 1e9;
    for (var i = 0u; i < params.solid; i = i + 1u) {
        let c = params.charges[i];
        d = min(d, distance(p, c.xy) - c.w);
    }
    for (var i = 0u; i < params.counts.y; i = i + 1u) {
        let m = params.magnets[i];
        d = min(d, distance(p, m.xy) - m.w);
    }
    for (var i = 0u; i < params.counts.z; i = i + 1u) {
        let l = params.loops[i];
        d = min(d, abs(distance(p, l.xy) - l.z) - params.wire);
    }
    for (var i = 0u; i < params.counts.w; i = i + 1u) {
        let s = params.segments[i];
        let ab = s.zw - s.xy;
        let t = clamp(dot(p - s.xy, ab) / dot(ab, ab), 0.0, 1.0);
        d = min(d, distance(p, s.xy + ab * t) - params.wire);
    }
    return d;
}

fn potential_colour(p: vec2<f32>) -> vec3<f32> {
    var u = 0.0;
    for (var i = 0u; i < params.counts.x; i = i + 1u) {
        let c = params.charges[i];
        var z = 0.0;
        if (i >= params.solid) {
            z = c.w;
        }
        let d = p - c.xy;
        u = u + c.z / max(sqrt(dot(d, d) + z * z), 1e-4);
    }
    if (params.moment_weight != 0.0) {
        u = u + params.moment_weight * magnetic_field(p);
    }
    u = u - params.u_a;

    let s = tanh(u / 1.5);
    let base = vec3<f32>(0.10, 0.11, 0.14);
    var col: vec3<f32>;
    if (s >= 0.0) {
        col = base + s * vec3<f32>(0.55, 0.10, 0.02);
    } else {
        col = base - s * vec3<f32>(0.02, 0.18, 0.55);
    }
    // Forbidden region (U > T0), antialiased edge.
    let fw_u = max(fwidth(u), 1e-6);
    let forbidden = smoothstep(-0.5 * fw_u, 0.5 * fw_u, u - 1.0);
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
    // Turning line U = T0.
    let turning = line_cover(abs(u - 1.0), fw_u, 2.0);
    return mix(col, vec3<f32>(0.95, 0.95, 0.85), turning);
}

fn magnetic_field(p: vec2<f32>) -> f32 {
    var b = 0.0;
    // Dipoles with moment along z: B_z = -mu / r^3 in the plane.
    for (var i = 0u; i < params.counts.y; i = i + 1u) {
        let m = params.magnets[i];
        let r = max(distance(p, m.xy), 1e-3);
        b = b - m.z / (r * r * r);
    }
    // Circular coils in the plane (z = 0): B_z = 2k / (alpha^2 beta) [(a^2 - rho^2) E + alpha^2 K].
    for (var i = 0u; i < params.counts.z; i = i + 1u) {
        let l = params.loops[i];
        let a = l.z;
        let rho = distance(p, l.xy);
        let alpha2 = max((a - rho) * (a - rho), 1e-8);
        let beta2 = (a + rho) * (a + rho);
        let ke = elliptic_ke(4.0 * a * rho / beta2, alpha2 / beta2);
        b = b + 2.0 * l.w / (alpha2 * sqrt(beta2)) * ((a * a - rho * rho) * ke.y + alpha2 * ke.x);
    }
    // Straight segments.
    for (var i = 0u; i < params.counts.w; i = i + 1u) {
        let s = params.segments[i];
        let ra = s.xy - p;
        let rb = s.zw - p;
        let la = length(ra);
        let lb = length(rb);
        let cross_z = ra.x * rb.y - ra.y * rb.x;
        let denom = max(la * lb * (la * lb + dot(ra, rb)), 1e-12);
        b = b + params.segment_kappa[i].x * cross_z * (la + lb) / denom;
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
