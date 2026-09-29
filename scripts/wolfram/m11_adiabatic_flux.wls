#!/usr/bin/env wolframscript
(* Test M11 (PHYSICS.md §2.2, Jackson §12.5): a charge gyrating in the plane of a loop
   (radius a = 20, current kappa rising 4x) with its guiding centre at rho0 ~ 2.
   The canonical angular momentum about the axis, L = rho p_phi + q rho A_phi, is conserved
   exactly (axial symmetry, test M10). Averaged over a gyration it is
   (q/2 pi) Phi(rho_gc) - (q/2) B r_g^2 (+ O(r_g^2 grad^2 B)), and B r_g^2 = 2 m mu/q^2 is the
   adiabatic invariant, so kappa(t) Phi(rho_gc) is conserved: the guiding centre keeps the
   loop field's flux through its circle. Phi(rho) = 2 pi rho A_phi(rho) with the loop's
   exact vector potential (units of the game: mu0/4 pi = 1, so mu0 I = 4 pi kappa):
     A_phi = (4 kappa/k) sqrt(a/rho) [(1 - k^2/2) K(k^2) - E(k^2)],  k^2 = 4 a rho/(a + rho)^2.
   Output: the guiding centre's start rho0 (as the test defines it), and the end radius
   rho_end solving 4 Phi(rho_end) = Phi(rho0): the adiabatic limit of the test's result. *)

prec = 30;
a = 20; kap0 = a/(2 Pi); q = 1; m = 1;
aPhi[rho_, kap_] := Module[{k2 = 4 a rho/(a + rho)^2},
  4 kap/Sqrt[k2] Sqrt[a/rho] ((1 - k2/2) EllipticK[k2] - EllipticE[k2])];
flux[rho_, kap_] := 2 Pi rho aPhi[rho, kap];
(* B_z in the plane from the flux: B = (1/(2 pi rho)) dPhi/drho *)
bz[rho_, kap_] := D[flux[r, kap], r]/(2 Pi r) /. r -> rho;

(* The test: b0 = B_z at (2, 0); p0 = 0.5 q b0 along +x, launched at (2, 0.5); its
   guiding centre (x + p x z/(q B)) with B at the launch point, rho = sqrt(2^2 + 0.5^2). *)
b0 = N[bz[2, kap0], prec];
p0 = b0/2;
bLaunch = N[bz[Sqrt[4 + 1/4], kap0], prec];
rho0 = Sqrt[4 + (1/2 - p0/(q bLaunch))^2];
rhoEnd = rho /. FindRoot[4 flux[rho, kap0] == flux[rho0, kap0], {rho, 1},
  WorkingPrecision -> prec];
Print["// guiding centre: start rho0, adiabatic end rho_end (flux kept), uniform-field 2/sqrt(4) = 1"];
Print[ToString[N[rho0, 17], CForm], ", ", ToString[N[rhoEnd, 17], CForm]];
Print["// B non-uniformity: B(rho0)/B(0) - 1 = ", ToString[N[bz[rho0, kap0]/bz[10^-8, kap0] - 1, 6], CForm]];
