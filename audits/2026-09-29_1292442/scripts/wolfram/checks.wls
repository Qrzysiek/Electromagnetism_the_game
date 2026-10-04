(* Independent checks for the audit in audit/FINDINGS.md.
   Dipole potential, uniform-sphere field, launch beta, the tapered velocity limit,
   and the coil brackets g(m) and h(m). Regenerated with:
   wolframscript -file audit/scripts/wolfram/checks.wls
   Wolfram Engine 14.2. *)

$fails = 0;
fail[msg_] := ($fails++; Print["FAILED: ", msg]);
show[label_, expr_] := Print[label, ": ", ToString[expr, InputForm]];

(* 1. Static dipole. phi = n.p/r^2, E = (3 n (n.p) - p)/r^3. *)
phiDipole = (px x + py y)/(x^2 + y^2)^(3/2);
eJackson = With[{r2 = x^2 + y^2, nd = px x + py y},
  { (3 x nd/r2 - px)/r2^(3/2), (3 y nd/r2 - py)/r2^(3/2) }];
gradDiff = Simplify[-Grad[phiDipole, {x, y}] - eJackson, x^2 + y^2 > 0];
If[gradDiff =!= {0, 0}, fail["dipole gradient"]];
show["dipole -grad(n.p/r^2) - E", gradDiff];

(* 2. Uniform sphere, k = 1. Continuous potential, E = -grad phi. *)
phiIn = Q (3 R^2 - x^2 - y^2)/(2 R^3);
phiOut = Q/Sqrt[x^2 + y^2];
eIn = Simplify[-Grad[phiIn, {x, y}]];
eOut = Simplify[-Grad[phiOut, {x, y}], x^2 + y^2 > 0];
onSurface = {x -> R Cos[th], y -> R Sin[th]};
surface = Simplify[(phiIn - phiOut) /. onSurface, R > 0];
eSurface = Simplify[(eIn - eOut) /. {x -> R, y -> 0}, R > 0];
If[surface =!= 0, fail["sphere potential"]];
If[eSurface =!= {0, 0}, fail["sphere field"]];
show["sphere phi_in - phi_out at r=R", surface];
show["sphere E_in", eIn];
show["sphere E_out", eOut];
show["sphere E_in - E_out at (R,0)", eSurface];

(* 3. beta from kinetic energy, and the taper thresholds. *)
beta = Sqrt[eps (eps + 2)]/(1 + eps);
gamma = 1 + eps;
betaFromGamma = FullSimplify[beta^2 - (1 - 1/gamma^2), eps > 0];
If[betaFromGamma =!= 0, fail["beta"]];
show["beta^2 - (1 - 1/gamma^2)", betaFromGamma];
(* Pure brake: |v_inf| = |v| + delta c. Equals c at beta = 1 - delta. *)
show["brake reaches c after a 0.1c shift", 1 - 1/10];
show["brake reaches c after a 0.07c shift", 1 - 7/100];

(* 4. Taper. T = 0.1 c/A, u = tau/T, u1 = 7/10, w = 3/10.
   accel = A for |u|<=u1, and A Sech[(|u|-u1)/w]^2 beyond.
   Integrate the positive side and compare with beam::tapered. *)
u1 = 7/10; w = 3/10;
vClosed[s_] := v0 + A T (u1 + w Tanh[s]);
vIntegrated = v0 + A (u1 T) + Integrate[A Sech[sp]^2 * (T w), {sp, 0, s}];
velDiff = Simplify[vClosed[s] - vIntegrated, s > 0];
If[velDiff =!= 0, fail["taper velocity"]];
show["taper velocity minus integral", velDiff];
show["taper future v as s->Infinity", Limit[vClosed[s], s -> Infinity]];
(* Past: s increases as tau decreases, d(tau) = -T w ds. *)
vPast[s_] := v0 - A T (u1 + w Tanh[s]);
vPastInt = v0 - A (u1 T) + Integrate[A Sech[sp]^2 * (-T w), {sp, 0, s}];
pastDiff = Simplify[vPast[s] - vPastInt, s > 0];
If[pastDiff =!= 0, fail["taper past"]];
show["taper past velocity minus integral", pastDiff];
show["taper past v as s->Infinity", Limit[vPast[s], s -> Infinity]];
(* Position from tau = 0, matching the code's disp for u > u1.
   x = x0 + v0 tau + A T^2 (u1^2/2 + u1 w s + w^2 Log[Cosh[s]]). *)
tauOf[s_] := T (u1 + w s);
xCode[s_] := x0 + v0 tauOf[s] + A T^2 (u1^2/2 + u1 w s + w^2 Log[Cosh[s]]);
xInt = x0 + v0 (u1 T) + A (u1 T)^2/2 + Integrate[vClosed[sp] * (T w), {sp, 0, s}];
posDiff = Simplify[xCode[s] - xInt, s > 0];
If[posDiff =!= 0, fail["taper position"]];
show["taper position minus integral", posDiff];

(* 5. Coil brackets. Parameter m = k^2, Mathematica's EllipticK/EllipticE. *)
g[m_] := (1 - m/2) EllipticE[m] - (1 - m) EllipticK[m];
h[m_] := (1 - m/2) EllipticK[m] - EllipticE[m];
gSeries = Series[g[m], {m, 0, 4}];
hSeries = Series[h[m], {m, 0, 4}];
show["g(m) series", gSeries];
show["h(m) series", hSeries];
show["g leading / (3 Pi m^2/32)", Limit[g[m]/(3 Pi m^2/32), m -> 0]];
show["h leading / (Pi m^2/32)", Limit[h[m]/(Pi m^2/32), m -> 0]];
If[SeriesCoefficient[gSeries, 0] =!= 0 || SeriesCoefficient[gSeries, 1] =!= 0, fail["g low orders"]];
If[SeriesCoefficient[hSeries, 0] =!= 0 || SeriesCoefficient[hSeries, 1] =!= 0, fail["h low orders"]];

(* 6. Cup illustration, the same ODE as audit/scripts/audit_checks.py.
   Survivor starts at rest at the origin. Absorbed charge starts at (2,0),
   speed 1 along +x, factor Exp[-Pi t/2]. *)
sol = NDSolve[{
  x'[t] == vx[t], y'[t] == vy[t],
  vx'[t] == (x[t] - (2 + t)) Exp[-Pi t/2] / ((x[t] - (2 + t))^2 + y[t]^2)^(3/2),
  vy'[t] == (y[t] - 0) Exp[-Pi t/2] / ((x[t] - (2 + t))^2 + y[t]^2)^(3/2),
  x[0] == 0, y[0] == 0, vx[0] == 0, vy[0] == 0
  }, {x, y, vx, vy}, {t, 0, 8},
  Method -> "StiffnessSwitching", PrecisionGoal -> 10, AccuracyGoal -> 10];
If[!MatchQ[sol, {{_Rule ..}}], fail["ndsolve"]];
ke = (vx[8]^2 + vy[8]^2)/2 /. sol[[1]];
If[!NumericQ[ke], fail["cup kinetic"]];
show["cup survivor kinetic at t=8", N[ke, 17]];
If[$fails > 0, Print["wolfram checks failed: ", $fails]; Exit[1]];
Print["wolfram checks finished"];
