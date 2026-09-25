// Potential map, computed per pixel (visual only, f32).
//
// Colour shows U / T0 = q (phi - phi_A) / T0 from the particle's point of view: red uphill,
// blue downhill, contours every T0/4. U > T0 is forbidden by energy conservation (dark);
// its boundary, the turning line U = T0, is bright. Lines are antialiased with screen-space
// derivatives, so the map is sharp at any zoom.

#import bevy_sprite::mesh2d_vertex_output::VertexOutput

struct Params {
    // xy: position (cells), z: w = q Q / T0, w: sphere radius (cells).
    charges: array<vec4<f32>, 256>,
    count: u32,
    // U / T0 at the launch point.
    u_a: f32,
};

@group(#{MATERIAL_BIND_GROUP}) @binding(0) var<uniform> params: Params;

fn srgb_to_linear(c: vec3<f32>) -> vec3<f32> {
    return pow(c, vec3<f32>(2.2));
}

// Coverage of a line of `width_px` pixels at distance `d` (in the units of `fw`).
fn line_cover(d: f32, fw: f32, width_px: f32) -> f32 {
    return 1.0 - smoothstep(0.5 * width_px * fw, (0.5 * width_px + 1.0) * fw, d);
}

@fragment
fn fragment(in: VertexOutput) -> @location(0) vec4<f32> {
    let p = in.world_position.xy;
    var u = 0.0;
    var inside = 1e9; // signed distance to the nearest sphere surface
    for (var i = 0u; i < params.count; i = i + 1u) {
        let c = params.charges[i];
        let d = distance(p, c.xy);
        inside = min(inside, d - c.w);
        u = u + c.z / max(d, 1e-4);
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

    // Contours every T0/4; skipped where they would be denser than about 3 pixels.
    let k = u * 4.0;
    let fw_k = fwidth(k);
    if (fw_k < 0.3) {
        let f = fract(k);
        let d = min(f, 1.0 - f);
        let cover = line_cover(d, fw_k, 1.0) * (1.0 - forbidden);
        col = mix(col, col + vec3<f32>(0.10), cover);
    }

    // Turning line U = T0.
    let turning = line_cover(abs(u - 1.0), fw_u, 2.0);
    col = mix(col, vec3<f32>(0.95, 0.95, 0.85), turning);

    // Inside a charge sphere.
    let fw_p = max(fwidth(inside), 1e-6);
    let solid = 1.0 - smoothstep(-0.5 * fw_p, 0.5 * fw_p, inside);
    col = mix(col, vec3<f32>(0.08, 0.08, 0.09), solid);

    return vec4<f32>(srgb_to_linear(clamp(col, vec3<f32>(0.0), vec3<f32>(1.0))), 1.0);
}
