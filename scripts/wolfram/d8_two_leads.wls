#!/usr/bin/env wolframscript
(* Test `two_leads_make_the_loop_inductance` (PHYSICS.md §2.12, DEC-27): two straight
   antiparallel leads in the z-invariant world (current sheets along z, current I and -I
   along y), length l, d apart, in the infinite plane, their field in the Darwin
   (Coulomb-gauge) limit. B_z = mu0 psi, psi the stream function of the currents'
   transverse part: -Laplacian psi = (curl J)_z, a dipole layer along each lead. The
   energy W = (mu0/2) Integral psi^2 = (mu0/2) Integral Integral w G w' with G = r^2 ln r/(8 Pi)
   the Green's function of the squared Laplacian (the constant parts cancel between the two
   leads: their currents sum to zero), so with u = y - y':
     L = 2W/I^2 = (mu0/(2 Pi)) Integral_0^l Integral_0^l
                    [ (1/2) ln(1 + d^2/u^2) + d^2/(u^2 + d^2) ] dy dy'
       = (mu0/Pi) Integral_0^l (l - u) f(u) du,
   which tends to mu0 d l for long leads (Integral_R f = 2 Pi d). Checked independently by
   summing psi^2 on a fine grid from the angles the leads subtend, psi = (Omega1 - Omega2)/(2 Pi)
   (target/scratch/two_leads_reference.py: 0.13 % below, its tail cut at 40).
   Usage: wolframscript -file d8_two_leads.wls [d] [l] [c]  (defaults 4, 24, 3). Output: the
   closed form and its value to 15 digits. *)

args = Rest[$ScriptCommandLine];
d = If[Length[args] >= 1, ToExpression[args[[1]]], 4];
l = If[Length[args] >= 2, ToExpression[args[[2]]], 24];
c = If[Length[args] >= 3, ToExpression[args[[3]]], 3];
mu0 = 4 Pi/c^2;
f[u_] := Log[1 + dd^2/u^2]/2 + dd^2/(u^2 + dd^2);
closed = Simplify[(1/Pi) Integrate[(ll - u) f[u], {u, 0, ll},
    Assumptions -> dd > 0 && ll > 0], dd > 0 && ll > 0];
Print["L/mu0 = ", closed];
value = N[mu0 closed /. {dd -> d, ll -> l}, 15];
Print["d = ", d, ", l = ", l, ", c = ", c, ": L = ", value,
  " (mu0 d l = ", N[mu0 d l, 15], ")"];
