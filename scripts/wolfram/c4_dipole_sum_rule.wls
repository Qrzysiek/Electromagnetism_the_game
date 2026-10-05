#!/usr/bin/env wolframscript
(* Tests C3 and C4 (PHYSICS.md §2.1; Jackson §16.8 and Problem 16.13): a charge q, mass m,
   bound harmonically (omega0), with radiation reaction, driven by E0 e^(-i w t).
   Jackson's equation of motion (16.73) with no other damping (Gamma' = 0) is the
   Landau-Lifshitz reduction (his (16.10)): x'' + Gamma x' + omega0^2 x
   = (q E0/m)(1 - i w tau) e^(-i w t), Gamma = omega0^2 tau, tau = 2 q^2/(3 m c^3).
   Its polarizability alpha = q x/E0; Problem 16.13(a): sigma_t = (4 pi w/c) Im alpha.
   Checks: (1) sigma_t is Jackson's (16.78) with Gamma' = 0, = sigma_T w^4/D;
   (2) sigma_t -> sigma_T at high frequency (Thomson's plateau), so the integral of
   sigma_t diverges: the premise of Problem 16.13(b), alpha -> -q^2/(m w^2), fails by
   the factor (1 - i w tau); (3) the integral of sigma_t - sigma_T over (0, oo), exactly;
   (4) its value for the tests' parameters; (5) aside: Abraham-Lorentz's steady state
   (x''' -> i w^3 x), whose plateau ends at w ~ 1/tau, integrated. *)

ClearAll["Global`*"];
assume = q != 0 && m > 0 && c > 0 && w0 > 0 && tau > 0 && w > 0;
gam = w0^2 tau;
alpha = (q^2/m) (1 - I w tau)/(w0^2 - w^2 - I w gam);
sigmaT = (8 Pi/3) (q^2/(m c^2))^2;
sigmaLL = FullSimplify[(4 Pi w/c) ComplexExpand[Im[alpha]] /. tau -> 2 q^2/(3 m c^3), assume];
d = (w0^2 - w^2)^2 + w^2 gam^2 /. tau -> 2 q^2/(3 m c^3);
Print["(1) sigma_t - sigma_T w^4/D = ",
  FullSimplify[sigmaLL - sigmaT w^4/d, assume]];
Print["(2) high-frequency limit of sigma_t/sigma_T: ",
  Limit[sigmaLL/sigmaT, w -> Infinity, Assumptions -> assume]];
(* (3) In terms of g = Gamma: sigma_t - sigma_T = sigma_T (w^4/D - 1). *)
(* With w = omega0 x it is omega0 times a function of gg = Gamma/omega0; underdamped
   (0 < gg < 2; the tests have gg = 1.3e-3), gg = 2 sin(ph) lets FullSimplify close it. *)
int1 = Integrate[x^4/((1 - x^2)^2 + gg^2 x^2) - 1, {x, 0, Infinity}, Assumptions -> 0 < gg < 2];
Print["(3) Integral of (w^4/D - 1) dw over (0, oo), minus pi (omega0^2 - Gamma^2)/(2 Gamma), ",
  "for 0 < Gamma < 2 omega0: omega0 times ",
  FullSimplify[(int1 - Pi (1 - gg^2)/(2 gg)) /. gg -> 2 Sin[ph], 0 < ph < Pi/2]];
Do[Module[{num},
    num = NIntegrate[w^4/((pw0^2 - w^2)^2 + w^2 pg^2) - 1, {w, 0, pw0 - 30 pg, pw0, pw0 + 30 pg, 4 pw0, Infinity},
      WorkingPrecision -> 40, PrecisionGoal -> 20, MaxRecursion -> 50];
    Print["    numerically at omega0 = ", pw0, ", Gamma = ", pg, ": ", N[num, 20],
      ", minus the formula ", N[num - Pi (pw0^2 - pg^2)/(2 pg), 3]]],
  {pw0, {1, 3/2}}, {pg, {1/768, 1/10, 1}}];
sumRule = 2 Pi^2 q^2/(m c);
val = sigmaT Pi (w0^2 - g^2)/(2 g) /. g -> w0^2 tau /. tau -> 2 q^2/(3 m c^3);
Print["    sigma_T times it, over 2 pi^2 q^2/(m c), in terms of tau: ",
  FullSimplify[(val/sumRule) /. q -> Sqrt[3 m c^3 tau/2], assume]];
(* (4) The tests' parameters: q = -1, m = 1, c = 8, omega0 = 1. *)
pars = {q -> -1, m -> 1, c -> 8, w0 -> 1};
Print["(4) q = -1, m = 1, c = 8, omega0 = 1: integral = ", N[val /. pars, 20],
  "; 2 pi^2 q^2/(m c) = ", N[sumRule /. pars, 20],
  "; ratio - 1 = ", N[(val/sumRule - 1) /. pars, 20]];
(* (5) Abraham-Lorentz: sigma = sigma_T w^4/((w0^2 - w^2)^2 + tau^2 w^6). *)
tauV = 2/(3 8^3);
alInt = NIntegrate[w^4/((1 - w^2)^2 + tauV^2 w^6), {w, 0, 1 - 50 tauV, 1, 1 + 50 tauV, 10, 1/tauV, Infinity},
  WorkingPrecision -> 30, PrecisionGoal -> 12, MaxRecursion -> 40];
Print["(5) Abraham-Lorentz, same parameters: integral of sigma over (0, oo) / (2 pi^2 q^2/(m c)) = ",
  N[(8 Pi/3) (1/64)^2 alInt/(2 Pi^2/8), 12]];
