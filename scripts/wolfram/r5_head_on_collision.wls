#!/usr/bin/env wolframscript
(* Test R5 (PHYSICS.md §3.1): a charge q (mass m) fired head-on at a repulsive fixed charge
   Q from x = -400 with kinetic energy T0 = m v0^2/2 (v0 = 0.8, c = 20), relativistic.
   Along a line a = F/(gamma^3 m), so Lienard's power (2 q^2/3 c^3) gamma^6 a^2 is
   (2 q^2/(3 m^2 c^3)) (q Q/r^2)^2, independent of gamma, and the energy radiated along the
   exact (radiation-free) trajectory is
     W = (2 q^4 Q^2/(3 m^2 c^3)) [Int_{rmin}^{400} + Int_{rmin}^{400.5}] dr/(r^4 v(r)),
   the way in from the launch and out to the arena's edge (the cube of half-size 400.5),
   with gamma(r) m c^2 = E - q Q/r and rmin = q Q/(T0 + q Q/400). The test flies it with
   radiation reaction: the trajectory then differs by O(W/T0) ~ 2e-5, which bounds the
   difference. Also printed: Jackson's nonrelativistic (8/45)(q/Q) m v_inf^5/c^3. *)

prec = 40;
q = 1; bigQ = 1; m = 1; c = 20; v0 = 8/10;
t0 = m v0^2/2;
(* The test launches with the relativistic momentum of kinetic energy T0. *)
e = m c^2 + t0 + q bigQ/400;
gam[r_] := (e - q bigQ/r)/(m c^2);
vel[r_] := c Sqrt[1 - 1/gam[r]^2];
rmin = q bigQ/(t0 + q bigQ/400);
pre = 2 q^4 bigQ^2/(3 m^2 c^3);
leg[rmax_] := NIntegrate[1/(r^4 vel[r]), {r, rmin, rmin (1 + 10^-3), rmin 2, 10, rmax},
  WorkingPrecision -> prec, PrecisionGoal -> 20, MaxRecursion -> 40];
w = pre (leg[400] + leg[801/2]);
vInf = Sqrt[v0^2 + 2 q bigQ/(m 400)];
jackson = 8/45 (q/bigQ) m vInf^5/c^3;
Print["// exact radiated energy (relativistic, radiation-free trajectory), Jackson Pr. 14.5b"];
Print[ToString[N[w, 17], CForm], ", ", ToString[N[jackson, 17], CForm],
  "  // ratio ", ToString[N[w/jackson, 10], CForm]];
