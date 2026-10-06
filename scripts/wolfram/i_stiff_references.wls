#!/usr/bin/env wolframscript
(* Tests I2-I4 (PHYSICS.md §5.4, the RADAU5 port): reference solutions of three classical
   stiff problems, computed with NDSolve at 30-digit working precision (independent of
   RADAU5), printed to 20 digits; then at 34 digits. Van der Pol and HIRES agree between
   the two to 6e-13 and 1e-13 (the tests use the 34-digit values); for ROBER NDSolve at 34
   digits failed its own error test at t = 7e8 (and extrapolated), so ROBER's reference is
   the 30-digit run, which agrees with the Test Set's published values to 14 digits.
   I4: HIRES (the Test Set for IVP Solvers): eight reactions, at t = 321.8122.
   I3: Robertson's chemical kinetics (ROBER, the Test Set):
       y1' = -0.04 y1 + 1e4 y2 y3, y2' = 0.04 y1 - 1e4 y2 y3 - 3e7 y2^2,
       y3' = 3e7 y2^2, y(0) = (1, 0, 0), at t = 1e11.
   I2: Van der Pol, eps = 1e-6 (Hairer & Wanner's dr1_radau5 driver):
       y1' = y2, y2' = ((1 - y1^2) y2 - y1)/eps, y(0) = (2, -0.66), at t = 2.
   (A first version at 40 digits and goal 24 had not finished Van der Pol after 25
   minutes.) *)

run[prec_, goal_] := Module[{hires, rob, vdp},
  hires = NDSolve[{
      y1'[t] == -171/100 y1[t] + 43/100 y2[t] + 832/100 y3[t] + 7/10000,
      y2'[t] == 171/100 y1[t] - 875/100 y2[t],
      y3'[t] == -1003/100 y3[t] + 43/100 y4[t] + 35/1000 y5[t],
      y4'[t] == 832/100 y2[t] + 171/100 y3[t] - 112/100 y4[t],
      y5'[t] == -1745/1000 y5[t] + 43/100 y6[t] + 43/100 y7[t],
      y6'[t] == -280 y6[t] y8[t] + 69/100 y4[t] + 171/100 y5[t] - 43/100 y6[t] +
        69/100 y7[t],
      y7'[t] == 280 y6[t] y8[t] - 181/100 y7[t],
      y8'[t] == -280 y6[t] y8[t] + 181/100 y7[t],
      y1[0] == 1, y2[0] == 0, y3[0] == 0, y4[0] == 0, y5[0] == 0, y6[0] == 0,
      y7[0] == 0, y8[0] == 57/10000},
     {y1, y2, y3, y4, y5, y6, y7, y8}, {t, 0, 3218122/10000},
     WorkingPrecision -> prec, PrecisionGoal -> goal, AccuracyGoal -> goal + 4,
     MaxSteps -> Infinity];
  Print["I4 HIRES t=321.8122 (", prec, "): ",
    NumberForm[N[{y1[t], y2[t], y3[t], y4[t], y5[t], y6[t], y7[t], y8[t]} /.
        First[hires] /. t -> 3218122/10000, 20], 20]];
  rob = NDSolve[{y1'[t] == -4/100 y1[t] + 10^4 y2[t] y3[t],
      y2'[t] == 4/100 y1[t] - 10^4 y2[t] y3[t] - 3 10^7 y2[t]^2,
      y3'[t] == 3 10^7 y2[t]^2, y1[0] == 1, y2[0] == 0, y3[0] == 0},
     {y1, y2, y3}, {t, 0, 10^11}, WorkingPrecision -> prec, PrecisionGoal -> goal,
     AccuracyGoal -> goal + 14, MaxSteps -> Infinity];
  Print["I3 ROBER t=1e11 (", prec, "): ",
    NumberForm[N[{y1[10^11], y2[10^11], y3[10^11]} /. First[rob], 20], 20]];
  vdp = NDSolve[{y1'[t] == y2[t], y2'[t] == ((1 - y1[t]^2) y2[t] - y1[t]) 10^6,
      y1[0] == 2, y2[0] == -66/100}, {y1, y2}, {t, 0, 2},
     WorkingPrecision -> prec, PrecisionGoal -> goal, AccuracyGoal -> goal + 2,
     MaxSteps -> Infinity];
  Print["I2 VdP t=2 (", prec, "): ",
    NumberForm[N[{y1[2], y2[2]} /. First[vdp], 20], 20]];
];

run[30, 16];
run[34, 19];
