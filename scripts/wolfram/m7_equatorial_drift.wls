#!/usr/bin/env wolframscript
(* Test M7 (PHYSICS.md §2.2, Jackson Pr. 12.9b): a charge (q = m = 1, c = infinity)
   gyrating in the equatorial plane of a dipole, B_z = mu/r^3 (mu = 1e5; the test's moment
   is -1e5 z, which makes B_z positive in the plane). Launched at r0 = R - a = 9.9 on the
   x axis with velocity (0, v), v = a omega, omega = mu/R^3 = 100 (a = 0.1, R = 10).
   Exact: the speed is constant and, with A_phi = -mu/r^2 (B_z = (1/r) d(r A_phi)/dr),
   the canonical angular momentum L = r v_phi + q r A_phi = r v_phi - mu/r is conserved:
     v_phi(r) = (L + mu/r)/r,  v_r^2 = v^2 - v_phi^2,
   so the radial motion runs between the roots r1 < r2 of v_phi^2 = v^2, with period
   T = 2 Int dr/|v_r| and azimuth gain Dphi = 2 Int (v_phi/r) dr/|v_r|. The longitude
   drifts at Dphi/T exactly. Jackson's leading term: (3/2)(a/R)^2 omega. *)

prec = 40;
mu = 10^5; bigR = 10; a = 1/10; omega = mu/bigR^3; v = a omega;
r0 = bigR - a;
L = r0 v - mu/r0;
vphi[r_] := (L + mu/r)/r;
(* turning points: vphi = +v or -v; the orbit through r0 (where vphi = v, a turning point) *)
roots = r /. NSolve[{vphi[r]^2 == v^2, r > 0}, r, WorkingPrecision -> prec];
r1 = Max[Select[roots, # <= r0 (1 + 10^-20) &]];
r2 = Min[Select[roots, # > r0 (1 + 10^-20) &]];
(* r = rm + dr sin(th), th from -pi/2 to pi/2: dr/|v_r| is then regular *)
rm = (r1 + r2)/2; dr = (r2 - r1)/2;
jac[th_] := Module[{r = rm + dr Sin[th]}, dr Cos[th]/Sqrt[v^2 - vphi[r]^2]];
period = 2 NIntegrate[jac[th], {th, -Pi/2, Pi/2}, WorkingPrecision -> prec, PrecisionGoal -> 25];
dphi = 2 NIntegrate[vphi[rm + dr Sin[th]]/(rm + dr Sin[th]) jac[th], {th, -Pi/2, Pi/2},
  WorkingPrecision -> prec, PrecisionGoal -> 25];
drift = dphi/period;
jackson = 3/2 (a/bigR)^2 omega;
Print["// r1, r2 (turning points), exact |drift| dphi/dt, Jackson (3/2)(a/R)^2 omega, ratio"];
Print[ToString[N[r1, 17], CForm], ", ", ToString[N[r2, 17], CForm], ", ",
  ToString[N[Abs[drift], 17], CForm], ", ", ToString[N[jackson, 17], CForm], ", ",
  ToString[N[Abs[drift]/jackson, 12], CForm]];
