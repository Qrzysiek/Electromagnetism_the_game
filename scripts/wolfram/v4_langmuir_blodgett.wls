#!/usr/bin/env wolframscript
(* Test V4 (PHYSICS.md 2.11; docs/TUBES.md): the space-charge-limited current of a
   coaxial diode, cathode radius a inside anode radius b, per unit length (Langmuir &
   Blodgett 1923): I/L = (8 pi eps0/9) Sqrt[2 e/m] V^(3/2) / (b beta^2), i.e. in the
   game's units (k = 1, eps0 = 1/4 pi) (2/9) Sqrt[2 |q/m|] V^(3/2) / (b beta^2), with
   beta(gamma), gamma = Log[r/a], the solution of
       3 beta beta'' + beta'^2 + 4 beta beta' + beta^2 - 1 = 0,
   beta ~ gamma near the cathode (beta(0) = 0, beta'(0) = 1). The ODE is singular at
   gamma = 0: started at gamma0 from the series beta = g - 2/5 g^2 + 11/120 g^3
   - 47/3300 g^4 (the known series, checked here against the ODE's residual).
   Output: beta^2 at b/a = 2 and 3, to 12 digits. Usage: wolframscript -file ... [smoke] *)

smoke = MemberQ[$ScriptCommandLine, "smoke"];
prec = If[smoke, MachinePrecision, 30];
ser[g_] := g - 2/5 g^2 + 11/120 g^3 - 47/3300 g^4;
g0 = If[smoke, 1/100, 1/10000];
res = 3 b[g] b''[g] + b'[g]^2 + 4 b[g] b'[g] + b[g]^2 - 1;
Print["series residual at g0: ", N[(res /. b -> ser) /. g -> g0, 6]];
sol = NDSolve[{res == 0, b[g0] == ser[g0], b'[g0] == (D[ser[x], x] /. x -> g0)}, b,
  {g, g0, Log[4]}, WorkingPrecision -> prec, PrecisionGoal -> If[smoke, 8, 20],
  AccuracyGoal -> If[smoke, 8, 20], MaxSteps -> 10^6];
beta2[r_] := (b[Log[r]] /. First[sol])^2;
Do[
  If[! NumericQ[beta2[r]], Print["not numeric"]; Abort[]];
  Print["b/a = ", r, ": beta^2 = ", ToString[N[beta2[r], 14], CForm]],
  {r, {2, 3}}];
