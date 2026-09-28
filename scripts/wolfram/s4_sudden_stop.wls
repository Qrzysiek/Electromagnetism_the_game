#!/usr/bin/env wolframscript
(* Test S4 (PHYSICS.md §3.4): the spectrum of a charge stopped smoothly, velocity
   v0 s(t) with s = 1/(1 + e^(t/tau)), against the sudden stop at t = 0.

   Along a fixed line the radiation vector n x (n x beta)/kappa has a fixed direction, so
   Jackson's amplitude (14.65, the integrand is the derivative of F = n x (n x beta)/kappa)
   is the scalar  A = Integral f'(t) e^(i w phi(t)) dt,  f = beta sin(th)/(1 - beta cos(th)),
   phi = t - n.x/c. Divided by the sudden stop's  -f0 e^(i w phi(0)):

     R(w) = Integral W(t) e^(i w psi(t)) dt,  W = f'/(-f0) (Integral W = 1),
     psi(t) = phi(t) - phi(0) = Integral_0^t (1 - beta(t') cos th) dt',

   so |R|^2 = 1 - w^2 Var_W(psi) + O(w^4). With u = t/tau, Var_W(psi) = tau^2 V(beta0, th):
   the relative difference of the smooth stop's spectrum is -(w tau)^2 V to leading order.
   Here: V for beta0 = 0.8 at the test's angles, and the exact |R|^2 at the test's w tau
   (by quadrature) to show the size of the higher orders. *)

s[u_] := 1/(1 + Exp[u]);
b0 = 8/10;
f[u_, th_] := b0 s[u] Sin[th]/(1 - b0 s[u] Cos[th]);
f0[th_] := b0 Sin[th]/(1 - b0 Cos[th]);
w[u_, th_] := -D[f[x, th], x]/f0[th] /. x -> u;
(* psi/tau = u - b0 cos(th) (u - Log[1 + E^u] + Log[2]): the integral of 1 - b0 s cos th. *)
psi[u_, th_] := u - b0 Cos[th] (u - Log[1 + Exp[u]] + Log[2]);

moments[th_] := Module[{m0, m1, m2},
  m0 = NIntegrate[w[u, th], {u, -80, 0, 80}, WorkingPrecision -> 40, PrecisionGoal -> 25];
  m1 = NIntegrate[w[u, th] psi[u, th], {u, -80, 0, 80}, WorkingPrecision -> 40, PrecisionGoal -> 25];
  m2 = NIntegrate[w[u, th] psi[u, th]^2, {u, -80, 0, 80}, WorkingPrecision -> 40, PrecisionGoal -> 25];
  {m0, m2/m0 - (m1/m0)^2}];

exact[th_, wt_] := Module[{re, im},
  re = NIntegrate[w[u, th] Cos[wt psi[u, th]], {u, -80, 0, 80}, WorkingPrecision -> 40, PrecisionGoal -> 25];
  im = NIntegrate[w[u, th] Sin[wt psi[u, th]], {u, -80, 0, 80}, WorkingPrecision -> 40, PrecisionGoal -> 25];
  re^2 + im^2 - 1];

(* The tails beyond |u| = 80 are below e^-80. Output: the table for the test (Rust). *)
Print["// theta (deg), V = Var(psi/tau), then |R|^2 - 1 at w tau = 1e-3, 4e-3, 1.6e-2"];
Do[
  th = deg Degree;
  {m0, v} = moments[th];
  vals = Table[exact[th, wt], {wt, {1/1000, 4/1000, 16/1000}}];
  Print["(", ToString[N[deg, 6], CForm], ", ", ToString[N[v, 17], CForm], ", [",
    StringRiffle[ToString[N[#, 17], CForm] & /@ vals, ", "], "]),"],
  {deg, {10, 3687/100, 60, 120}}];
