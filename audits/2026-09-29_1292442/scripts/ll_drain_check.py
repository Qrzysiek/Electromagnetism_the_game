"""Central difference of landau_lifshitz across an instant-drain light cone.

Transcribes the stencil in crates/physics/src/dynamics.rs (h = 1e-5 / max(|v|, 1),
samples at x ± v h, t ± h) and the cutoff in crates/physics/src/beam.rs fields_of:
c * (t - t_off) >= |x - x_off| returns a zero field. This does not call the crate.
"""

from __future__ import annotations

H_NUM = 1e-5


def stencil_h(speed: float) -> float:
    return H_NUM / max(abs(speed), 1.0)


def drain_e(x: float, t: float, c: float, t_off: float, x_off: float, e0: float) -> float:
    if c * (t - t_off) >= abs(x - x_off):
        return 0.0
    return e0


def central(x: float, t: float, v: float, c: float, t_off: float, x_off: float, e0: float):
    h = stencil_h(v)
    ep = drain_e(x + v * h, t + h, c, t_off, x_off, e0)
    em = drain_e(x - v * h, t - h, c, t_off, x_off, e0)
    return (ep - em) / (2.0 * h), h, ep, em


def main() -> None:
    c, t_off, x_off, e0 = 1.0, 0.0, 0.0, 1.0
    for v in (0.0, 0.5):
        de, h, ep, em = central(1.0, 1.0, v, c, t_off, x_off, e0)
        expect = -e0 / (2.0 * h)
        print(
            f"on cone v={v}: de {de:.17g} expect {expect:.17g} "
            f"residual {abs(de - expect):.3e} samples {ep} {em} h {h:.3e}"
        )
    de_before, _, ep_b, em_b = central(1.0, 0.5, 0.0, c, t_off, x_off, e0)
    de_after, _, ep_a, em_a = central(1.0, 2.0, 0.0, c, t_off, x_off, e0)
    print(f"before cone: de {de_before:.3e} samples {ep_b} {em_b}")
    print(f"after cone: de {de_after:.3e} samples {ep_a} {em_a}")
    h = stencil_h(0.0)
    de_lin = ((1.0 + h) - (1.0 - h)) / (2.0 * h)
    print(f"linear field residual {abs(de_lin - 1.0):.3e}")


if __name__ == "__main__":
    main()
