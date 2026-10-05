#!/usr/bin/env wolframscript
(* Tests Z2 and Z3 (PHYSICS.md §2.2; Jackson §5.17, Problems 5.28, 5.32, 5.34): references
   for the inductances of thin coplanar circular loops, in units of mu0/4pi (the game's
   1/c^2), i.e. the Neumann integral N = (closed loop)(closed loop) dl1.dl2/|x1 - x2|.
   Independent of the game's elliptic integrals: the double integrals are taken directly.
   (1) Concentric coplanar loops, radii a = 3, b = 5: the Neumann double integral, and
       Jackson's (5.28) closed form with d = 0, 4 pi Sqrt[a b] ((2/k - k) K - (2/k) E),
       k^2 = 4ab/(a + b)^2 (Mathematica's EllipticK takes m = k^2).
   (2) Identical coplanar loops (radius a = 1) with centres R = 2.5, 5, 10 apart: the
       Neumann double integral, and Problem 5.34(c)'s series
       -pi^2 a (a/R)^3 (1 + (9/4)(a/R)^2 + (375/64)(a/R)^4), i.e. -(mu0 pi a/4)(...).
   (3) A loop's own flux with a uniform surface current on a wire of radius b (Problem
       5.32): the flux 2 pi rho A_phi of the filament through the circle at (rho, z) =
       (a + b cos psi, b sin psi), averaged over psi; its expansion in b/a, against
       4 pi a (Log[8a/b] - 2). *)

prec = 30;
neumann[a_, b_, dd_] := NIntegrate[
  a b Cos[t1 - t2]/Sqrt[a^2 + b^2 + dd^2 - 2 a dd Cos[t1] + 2 b dd Cos[t2] - 2 a b Cos[t1 - t2]],
  {t1, 0, 2 Pi}, {t2, 0, 2 Pi},
  Method -> {"TrapezoidalRule", "Points" -> 400}, WorkingPrecision -> prec, PrecisionGoal -> 20];

(* (1) *)
{a1, b1} = {3, 5};
k2 = 4 a1 b1/(a1 + b1)^2;
closed = 4 Pi Sqrt[a1 b1] ((2/Sqrt[k2] - Sqrt[k2]) EllipticK[k2] - (2/Sqrt[k2]) EllipticE[k2]);
Print["(1) concentric a = 3, b = 5: Neumann ", N[neumann[a1, b1, 0], 20],
  "; Jackson (5.28) ", N[closed, 20]];

(* (2) *)
Do[Module[{n = neumann[1, 1, rr], s},
    s = -Pi^2 (1/rr)^3 (1 + 9/4 (1/rr)^2 + 375/64 (1/rr)^4);
    Print["(2) identical, R = ", rr, ": Neumann ", N[n, 20], "; Problem 5.34(c) series ",
      N[s, 20], " (relative ", N[s/n - 1, 3], ")"]],
  {rr, {5/2, 5, 10}}];

(* (3) The filament's A_phi per mu0 I/4pi, Jackson (5.37):
   A = (4/k) Sqrt[a/rho] ((1 - k^2/2) K - E), k^2 = 4 a rho/((a + rho)^2 + z^2). *)
aphi[a_, rho_, z_] := Module[{m = 4 a rho/((a + rho)^2 + z^2)},
   4/Sqrt[m] Sqrt[a/rho] ((1 - m/2) EllipticK[m] - EllipticE[m])];
flux[a_, b_] := NIntegrate[2 Pi (a + b Cos[p]) aphi[a, a + b Cos[p], b Sin[p]], {p, 0, 2 Pi},
    Method -> {"TrapezoidalRule", "Points" -> 400}, WorkingPrecision -> prec, PrecisionGoal -> 20]/(2 Pi);
(* The next order, found from these values and checked here: the difference is
   4 pi a (b/a)^2 ((1/8) ln(8a/b) + 1/16) + O((b/a)^4 ln): the last column tends to a
   constant (or grows like a logarithm), so the (b/a)^2 coefficient is exact. *)
Do[Module[{f = flux[1, bb], lead = 4 Pi (Log[8/bb] - 2)},
    Print["(3) a = 1, b = ", N[bb], ": flux ", N[f, 20], "; 4 pi a (ln(8a/b) - 2) ", N[lead, 20],
      "; difference / (4 pi a (b/a)^2) ", N[(f - lead)/(4 Pi bb^2), 12],
      "; minus ((1/8) ln(8a/b) + 1/16), over (b/a)^2: ",
      N[((f - lead)/(4 Pi bb^2) - (Log[8/bb]/8 + 1/16))/bb^2, 8]]],
  {bb, {1/10, 1/30, 1/100, 1/300, 1/1000}}];
