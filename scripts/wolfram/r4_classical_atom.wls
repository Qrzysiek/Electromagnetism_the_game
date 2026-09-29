#!/usr/bin/env wolframscript
(* Test R4 (PHYSICS.md §3.1, Jackson Pr. 16.2): a charge -1 (m = 1) around a fixed +Z
   (Z = 1, Coulomb constant 1), c = 8, radiating (Landau-Lifshitz). Adiabatic inspiral
   through relativistic circular orbits:
     circular orbit: gamma m v^2/r = Z/r^2, i.e. u v = Z/r with u = gamma v,
       so u^4 - (Z/r)^2 u^2/c^2 - (Z/r)^2 = 0;
     energy E(r) = gamma c^2 - Z/r;  power P(r) = (2/(3 c^3)) gamma^4 (v^2/r)^2
       (relativistic Larmor for circular motion; LL gives it up to O(tau^2), and the Schott
       energy of a circular orbit is constant);
     dr/dt = -P(r)/E'(r).
   Output: the launch (tangential momentum gamma v and radial velocity dr/dt at r0 = 5,
   so that the flight starts on the inspiral and no eccentricity is excited), and r(T) at
   T = 1500, i.e. the exact d(r^3) over the flight; Jackson's nonrelativistic
   r^3 = r0^3 - 6 Z tau t (tau = 2/(3 c^3)) for comparison. *)

prec = 40;
z = 1; c = 8; r0 = 5; tEnd = 1500;
u2[r_] := ((z/r)^2/c^2 + Sqrt[(z/r)^4/c^4 + 4 (z/r)^2])/2;
gam[r_] := Sqrt[1 + u2[r]/c^2];
vel[r_] := Sqrt[u2[r]]/gam[r];
energy[r_] := gam[r] c^2 - z/r;
power[r_] := 2/(3 c^3) gam[r]^4 (vel[r]^2/r)^2;
drdt[r_] := -power[r]/(D[energy[x], x] /. x -> r);
sol = NDSolve[{rr'[t] == drdt[rr[t]], rr[0] == r0}, rr, {t, 0, tEnd},
  WorkingPrecision -> prec, PrecisionGoal -> 25, AccuracyGoal -> 30, MaxSteps -> 10^6];
rEnd = rr[tEnd] /. First[sol];
tau = 2/(3 c^3);
Print["// launch at r0 = 5: tangential momentum gamma v, radial velocity dr/dt; r(T=1500); exact d(r^3)/dt mean; Jackson -6 Z tau"];
Print[ToString[N[Sqrt[u2[r0]], 17], CForm], ", ", ToString[N[drdt[r0], 17], CForm], ", ",
  ToString[N[rEnd, 17], CForm], ", ", ToString[N[(rEnd^3 - r0^3)/tEnd, 17], CForm], ", ",
  ToString[N[-6 z tau, 17], CForm]];
