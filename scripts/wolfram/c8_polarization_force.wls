#!/usr/bin/env wolframscript
(* Test C8 (PHYSICS.md §2.1; Brau §3.1.3 and Ex. 3.6): a charge z passes a neutral
   Thomson atom (cloud Q = 1, R = 4, electron q = -1, m = 1 at rest at its centre:
   omega0 = 1/8, polarizability alpha = q^2/(m omega0^2) = 64) at speed v along x, impact
   parameter b along y, from x = -L to x = +L (L = 20 b): the transverse impulse it
   receives from the dipole it induces.
   Linear response, the projectile on its straight line (infinitely heavy), with t0 = -L/v:
     E(t)  = field of z at the atom's centre = -z R(t)/|R|^3, R(t) = (v t, b),
     x(t)  = (q/(m omega0)) Int_{t0}^{t} sin(omega0 (t - s)) E(s) ds   (from rest at t0),
     p(t)  = q x(t), the induced dipole,
     F(t)  = z E_dip(R(t)) = z (3 (p.n) n - p)/|R|^3,
     Delta p_y = Int_{t0}^{-t0} F_y dt.
   The adiabatic limit (slow passage, b/v >> 1/omega0): the potential -alpha z^2/(2 r^4),
   Delta p_y = -3 pi alpha z^2/(4 v b^4). z = 1/800 keeps the swing at 1e-5 of b (its
   quadrupole changes the force at first order in swing/b). *)

prec = 30;
{z, q, m, w0, b} = {1/800, -1, 1, 1/8, 20};
alpha = q^2/(m w0^2);
Do[Module[{v = w0 b/xi, t0, ex, ey, xs, ys, fy, dp, adiabatic},
    t0 = -20 b/v;
    ex[s_] := -z v s/((v s)^2 + b^2)^(3/2);
    ey[s_] := -z b/((v s)^2 + b^2)^(3/2);
    (* x(t) by the convolution, written as Im[e^(i w0 t) Int e^(-i w0 s) E(s) ds] / (m w0) q:
       solve the oscillator ODE instead (equivalent, and NDSolve at 30 digits). *)
    sol = NDSolve[{xx''[t] == -w0^2 xx[t] + (q/m) ex[t], yy''[t] == -w0^2 yy[t] + (q/m) ey[t],
        xx[t0] == 0, yy[t0] == 0, xx'[t0] == 0, yy'[t0] == 0,
        pp'[t] == Module[{rx = v t, ry = b, r2, dx, dy, pdn},
          r2 = rx^2 + ry^2; dx = q xx[t]; dy = q yy[t];
          pdn = (dx rx + dy ry)/r2;
          z (3 pdn ry - dy)/r2^(3/2)],
        pp[t0] == 0},
      {xx, yy, pp}, {t, t0, -t0}, WorkingPrecision -> prec, PrecisionGoal -> 20,
      AccuracyGoal -> 40, MaxSteps -> Infinity];
    dp = pp[-t0] /. First[sol];
    adiabatic = -3 Pi alpha z^2/(4 v b^4);
    Print["xi = omega0 b/v = ", xi, ", v = ", N[v, 6], ": Delta p_y = ", N[dp, 20],
      "; adiabatic ", N[adiabatic, 20], " (relative ", N[dp/adiabatic - 1, 4], ")"]],
  {xi, {10, 25, 50}}];
