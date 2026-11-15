"""A floating contact's potential trace (`generator tube-scan --trace`): its spread and its
correlation time (the autocorrelation's integral to its first zero, and an exponential's
1/e time fitted to its first part)."""
import math
import sys


def load(path, t0):
    t, v = [], []
    for line in open(path, encoding="utf-8"):
        if line.startswith("trace"):
            _, a, b = line.split()
            if float(a) >= t0:
                t.append(float(a))
                v.append(float(b))
    return t, v


def analyse(path, t0=10.0, label=None):
    t, v = load(path, t0)
    n = len(v)
    dt = (t[-1] - t[0]) / (n - 1)
    mean = sum(v) / n
    d = [x - mean for x in v]
    var = sum(x * x for x in d) / n
    acf = []
    for k in range(n // 3):
        acf.append(sum(d[i] * d[i + k] for i in range(n - k)) / ((n - k) * var))
    # The integral to the first zero crossing.
    integral = 0.0
    for k, a in enumerate(acf):
        if a <= 0:
            break
        integral += a * dt if k else 0.5 * a * dt
    # The 1/e time.
    e_time = next((k * dt for k, a in enumerate(acf) if a < math.exp(-1)), float("nan"))
    span = t[-1] - t[0]
    print(f"{label or path}: n {n}, dt {dt:.4f}, mean {mean:+.4f}, rms {math.sqrt(var):.4f}, "
          f"tau(integral) {integral:.3f}, tau(1/e) {e_time:.3f}, window {span:.0f} "
          f"= {span / max(integral, 1e-9):.0f} tau; acf at 0.1,0.5,1,2,5: "
          + " ".join(f"{acf[min(len(acf) - 1, int(round(s / dt)))]:.2f}" for s in (0.1, 0.5, 1, 2, 5)))
    return math.sqrt(var), integral, e_time


if __name__ == "__main__":
    t0 = float(sys.argv[2]) if len(sys.argv) > 2 else 10.0
    analyse(sys.argv[1], t0)
