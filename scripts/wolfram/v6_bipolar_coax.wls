#!/usr/bin/env wolframscript
(* Test V6 (PHYSICS.md 2.11): the space-charge-limited bipolar coaxial diode. Electrons
   (|q/m| = 1) leave the inner cylinder r = a at phi = 0, ions (q/m = mu) the outer one
   r = b at phi = V, both space-charge limited (zero field at their emitter). Steady
   radial flow, k = 1 (Poisson: laplacian phi = -4 pi rho), currents per unit length I_e,
   I_i (magnitudes):
       (r phi')' = 2 I_e / Sqrt[2 phi] - 2 I_i / Sqrt[2 mu (V - phi)],
       phi(a) = 0, phi'(a) = 0, phi(b) = V, phi'(b) = 0.
   mu enters only through I_i / Sqrt[mu]: I_e does not depend on the ions' mass.
   Near each emitter the flow starts as the planar one: phi = C_a (r - a)^(4/3) with
   C_a = (9 I_e / (2 Sqrt[2] a))^(2/3), and V - phi = C_b (b - r)^(4/3) with
   C_b = (9 I_i / (2 Sqrt[2 mu] b))^(2/3). Integrated from both ends (started at a
   distance h with the leading term) and matched at the middle (phi and phi'), solved
   for I_e, I_i by FindRoot. Checked: with I_i = 0 the inner solution alone gives the
   unipolar Langmuir-Blodgett current (V4's beta^2), and h halved changes nothing in the
   printed digits. Output: I_e, I_i, and the enhancement I_e / I_LB, for a = 1, b = 2,
   V = 1, mu = 1. Usage: wolframscript -file v6_bipolar_coax.wls [smoke] *)

smoke = MemberQ[$ScriptCommandLine, "smoke"];
{a, b, V, mu} = {1, 2, 1, 1};
prec = If[smoke, MachinePrecision, 25];
goal = If[smoke, 8, 18];

inner[ie_?NumericQ, h_, rm_] := Module[{ca, s, sol},
  ca = (9 ie/(2 Sqrt[2] a))^(2/3);
  sol = NDSolve[{
     D[r p'[r], r] == 2 ie/Sqrt[2 p[r]],
     p[a + h] == ca h^(4/3), p'[a + h] == (4/3) ca h^(1/3)},
    p, {r, a + h, rm}, WorkingPrecision -> prec, PrecisionGoal -> goal,
    AccuracyGoal -> goal, MaxSteps -> 10^6];
  {p[rm], p'[rm]} /. First[sol]];

(* The full inner equation (with the ion term), for the bipolar case. *)
innerB[ie_?NumericQ, ii_?NumericQ, h_, rm_] := Module[{ca, sol},
  ca = (9 ie/(2 Sqrt[2] a))^(2/3);
  sol = NDSolve[{
     D[r p'[r], r] == 2 ie/Sqrt[2 p[r]] - 2 ii/Sqrt[2 mu (V - p[r])],
     p[a + h] == ca h^(4/3), p'[a + h] == (4/3) ca h^(1/3)},
    p, {r, a + h, rm}, WorkingPrecision -> prec, PrecisionGoal -> goal,
    AccuracyGoal -> goal, MaxSteps -> 10^6];
  {p[rm], p'[rm]} /. First[sol]];

outerB[ie_?NumericQ, ii_?NumericQ, h_, rm_] := Module[{cb, sol},
  cb = (9 ii/(2 Sqrt[2 mu] b))^(2/3);
  sol = NDSolve[{
     D[r p'[r], r] == 2 ie/Sqrt[2 p[r]] - 2 ii/Sqrt[2 mu (V - p[r])],
     p[b - h] == V - cb h^(4/3), p'[b - h] == (4/3) cb h^(1/3)},
    p, {r, rm, b - h}, WorkingPrecision -> prec, PrecisionGoal -> goal,
    AccuracyGoal -> goal, MaxSteps -> 10^6];
  {p[rm], p'[rm]} /. First[sol]];

h0 = If[smoke, 10^-4, 10^-8];

(* Unipolar check: I_LB = (2/9) Sqrt[2] V^(3/2) / (b beta^2), beta^2 = 0.27926716071059. *)
iLB = 2/9 Sqrt[2] V^(3/2)/(b 27926716071059/100000000000000);
(* Wrappers that evaluate only for numbers: Part on an unevaluated call would take its
   first argument (a first version of this script solved ie == V that way). *)
g0[ie_?NumericQ] := inner[ie, h0, b][[1]];
uni = FindRoot[g0[ie] == V, {ie, N[iLB, prec]},
   WorkingPrecision -> prec, AccuracyGoal -> goal - 2];
Print["unipolar I_e: ", N[ie /. uni, 15], " (Langmuir-Blodgett ", N[iLB, 15], ")"];

f1[ie_?NumericQ, ii_?NumericQ, h_] :=
  innerB[ie, ii, h, (a + b)/2][[1]] - outerB[ie, ii, h, (a + b)/2][[1]];
f2[ie_?NumericQ, ii_?NumericQ, h_] :=
  innerB[ie, ii, h, (a + b)/2][[2]] - outerB[ie, ii, h, (a + b)/2][[2]];
solve[h_] := FindRoot[{f1[ie, ii, h] == 0, f2[ie, ii, h] == 0},
   {{ie, 17/10 iLB}, {ii, iLB}}, WorkingPrecision -> prec,
   AccuracyGoal -> goal - 2];
s1 = solve[h0];
s2 = solve[h0/2];
{ie1, ii1} = {ie, ii} /. s1;
If[! NumericQ[ie1], Print["not numeric"]; Abort[]];
Print["bipolar I_e: ", ToString[N[ie1, 15], CForm], "  I_i: ", ToString[N[ii1, 15], CForm]];
Print["h halved:    ", N[{ie, ii} /. s2, 15]];
Print["enhancement I_e / I_LB: ", N[ie1/iLB, 12]];
