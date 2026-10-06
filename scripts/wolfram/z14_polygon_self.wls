#!/usr/bin/env wolframscript
(* Test Z14 (PHYSICS.md §2.8): reference for the self-inductance of a polygonal coil, in
   units of mu0/4pi, with its current uniformly on the wire's surface (radius b): the
   coil's own flux, i.e. the line integral of its centre line's vector potential per unit
   strength, around the closed curves on the wire's surface, averaged over the angle psi
   around the wire (as for the ring, z2_inductances.wls (3)). For the rectangle
   [0, A] x [0, B] (current counter-clockwise) the curve at psi is the rectangle
   [d, A - d] x [d, B - d] at the height h, with d = b cos psi, h = b sin psi.
   Independent of the engine's quadrature: each side's potential is the closed form
   Integral ds/|x - P - s u| = Log[(|Q - x| + (Q - x).u)/(|P - x| + (P - x).u)] (taken
   here in the numerically safe form), the line integrals by NIntegrate at 30 digits,
   and the psi average by the trapezoidal rule (exponentially convergent for this
   periodic analytic integrand), M = 64 and 128 points compared.
   Usage: wolframscript -file z14_polygon_self.wls [smoke] *)

smoke = MemberQ[$ScriptCommandLine, "smoke"];
prec = If[smoke, MachinePrecision, 30];
{aa, bb, b} = {4, 3, 1/10};
verts = {{0, 0, 0}, {aa, 0, 0}, {aa, bb, 0}, {0, bb, 0}};
segA[p_, q_, x_] := Module[{u = (q - p)/Sqrt[(q - p).(q - p)], ra = p - x, rb = q - x, la, lb, pa, pb},
  la = Sqrt[ra.ra]; lb = Sqrt[rb.rb]; pa = ra.u; pb = rb.u;
  u If[pa + pb >= 0, Log[(lb + pb)/(la + pa)], Log[(la - pa)/(lb - pb)]]];
potential[x_] := Sum[segA[verts[[k]], verts[[Mod[k, 4] + 1]], x], {k, 4}];
flux[psi_] := Module[{d = b Cos[psi], h = b Sin[psi], q},
  q = {{d, d, h}, {aa - d, d, h}, {aa - d, bb - d, h}, {d, bb - d, h}};
  Sum[With[{p0 = q[[k]], p1 = q[[Mod[k, 4] + 1]]},
      NIntegrate[potential[p0 + s (p1 - p0)].(p1 - p0), {s, 0, 1/2, 1},
        WorkingPrecision -> prec, PrecisionGoal -> If[smoke, 8, 20], MaxRecursion -> 30]],
    {k, 4}]];
avg[m_] := Sum[flux[2 Pi j/m], {j, 0, m - 1}]/m;
If[smoke,
  Print["smoke: flux(0) = ", flux[0], ", avg[8] = ", avg[8]],
  Module[{r64 = avg[64], r128},
    If[! NumericQ[r64], Print["not numeric"]; Abort[]];
    r128 = avg[128];
    Print["rectangle 4 x 3, b = 1/10: M = 64: ", ToString[N[r64, 20], CForm],
      "; M = 128: ", ToString[N[r128, 20], CForm]]]];
