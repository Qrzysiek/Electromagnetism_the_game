"""A square-wave-kicked contact's trace (`--trace`): folded at the period, the mean response,
its 1/e time after each edge (up and down edges averaged with their signs), and the noise
left when the mean response is subtracted: its spread and its correlation time."""
import math, sys
path, period, t0, t1 = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
pts = [tuple(map(float, l.split()[1:])) for l in open(path) if l.startswith("trace")]
pts = [(t, v) for t, v in pts if t0 <= t < t1]
dt = (pts[-1][0] - pts[0][0]) / (len(pts) - 1)
nb = int(round(period / dt))
sums, counts = [0.0] * nb, [0] * nb
for t, v in pts:
    b = int(round((t % period) / dt)) % nb
    sums[b] += v; counts[b] += 1
mean = [s / max(c, 1) for s, c in zip(sums, counts)]
# A half period's response, the down edge's with its sign flipped: average of the two halves.
half = nb // 2
resp = [(mean[i] - mean[(i + half) % nb]) / 2 for i in range(nb)]
p = max(range(nb), key=lambda i: abs(resp[i]))
peak = resp[p]; level = abs(peak) / math.e
tau = None
for k in range(1, half):
    a, b = abs(resp[(p + k - 1) % nb]), abs(resp[(p + k) % nb])
    if b <= level:
        tau = (k - 1 + (a - level) / (a - b)) * dt
        break
resid = []
for t, v in pts:
    b = int(round((t % period) / dt)) % nb
    resid.append(v - mean[b])
m = sum(resid) / len(resid)
d = [x - m for x in resid]
c0 = sum(x * x for x in d) / len(d)
pace = None; prev = 1.0
for k in range(1, len(d) // 4):
    rho = sum(d[i] * d[i + k] for i in range(len(d) - k)) / ((len(d) - k) * c0)
    if rho < 1 / math.e:
        pace = (k - 1 + (prev - 1 / math.e) / (prev - rho)) * dt
        break
    prev = rho
print(f"{path}: peak {peak:.3f}, decay(1/e) {tau}, residual spread {math.sqrt(c0):.4f}, residual pace {pace}")
