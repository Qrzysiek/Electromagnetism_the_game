#!/usr/bin/env wolframscript
(* Test C5 (PHYSICS.md §2.1; Jackson Problems 13.2-13.3; Brau §5.2.3): a charge z = 16 passes
   at speed v (c = 100, v = 80: gamma = 5/3) and impact parameter b a charge q = -1 (m = 1)
   bound harmonically (omega0 = 1/8) and at rest. Problem 13.2: the energy transferred is
   (q^2/2m) |Int E(t) e^(i omega0 t) dt|^2 (an undamped oscillator from rest, in the
   dipole approximation). The field of the uniformly moving charge at the oscillator:
   E_perp = z gamma b/(b^2 + gamma^2 v^2 t^2)^(3/2), E_par = -z gamma v t/(...)^(3/2).
   Problem 13.3: with xi = omega0 b/(gamma v),
     Delta E = (2 z^2 q^2/(m b^2 v^2)) (xi^2 K1(xi)^2 + xi^2 K0(xi)^2/gamma^2).
   For xi = 0.3, 1, 2, 4 (b = xi gamma v/omega0): the closed form, and the Fourier integral
   of the field taken numerically (independent of the Bessel-function algebra). The scale
   (c = 100) keeps the electron slow (v_e/c ~ 1e-5) at swings of 1e-2: the projectile's
   magnetic force on it changes the transfer at first order in v_e v/c^2. *)

prec = 40;
{z, q, m, c, v, w0} = {16, -1, 1, 100, 80, 1/8};
gam = 1/Sqrt[1 - v^2/c^2];
Do[Module[{b = xi gam v/w0, closed, fperp, fpar, numeric},
    closed = 2 z^2 q^2/(m b^2 v^2) (xi^2 BesselK[1, xi]^2 + xi^2 BesselK[0, xi]^2/gam^2);
    (* E_perp is even in t, E_par odd: Int E e^(i w t) = 2 Int_0^oo E_perp cos, 2 i Int_0^oo E_par sin. *)
    fperp = 2 NIntegrate[z gam b/(b^2 + gam^2 v^2 t^2)^(3/2) Cos[w0 t], {t, 0, Infinity},
        Method -> "LevinRule", WorkingPrecision -> prec, PrecisionGoal -> 25, MaxRecursion -> 60];
    fpar = 2 NIntegrate[z gam v t/(b^2 + gam^2 v^2 t^2)^(3/2) Sin[w0 t], {t, 0, Infinity},
        Method -> "LevinRule", WorkingPrecision -> prec, PrecisionGoal -> 25, MaxRecursion -> 60];
    numeric = q^2/(2 m) (fperp^2 + fpar^2);
    Print["xi = ", N[xi], ", b = ", N[b, 20], ": Delta E = ", N[closed, 20],
      "; the Fourier integral ", N[numeric, 20], " (relative ", N[numeric/closed - 1, 3],
      "); the impulse limit 2 z^2 q^2/(m b^2 v^2) ", N[2 z^2 q^2/(m b^2 v^2), 6]]],
  {xi, {3/10, 1, 2, 4}}];
