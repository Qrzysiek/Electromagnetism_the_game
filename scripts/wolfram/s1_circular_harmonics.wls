#!/usr/bin/env wolframscript
(* Test S1 (PHYSICS.md §3.4): a charge on a circle (radius rho, angular frequency w0,
   beta = rho w0/c) for exactly N turns from t = 0, seen in its orbital plane along
   n = x. Jackson 14.65: d2I/dw dOmega = (q^2/4 pi^2 c) |A(w)|^2,
     A(w) = Integral_0^T f(t) e^(i w (t - n.r/c)) dt,  f = n x ((n - beta) x betadot)/kappa^2,
   (in the plane f is along y). With g(t) = f_y(t) e^(-i w (rho/c) cos(w0 t)), periodic,
   g = Sum_m g_m e^(i m w0 t), and T = N 2 pi/w0:
     A(w) = Sum_m g_m(w) (e^(i w T) - 1)/(i (w + m w0)).
   g_m by the discrete Fourier transform of one period (exponentially accurate: g is
   analytic and periodic), the band energy by quadrature in w. This includes the finite
   flight's leakage exactly; the test previously compared with the infinite-time power
   T (q^2 m^2 w0^2 beta^2/2 pi c) J'_m(m beta)^2, which it approaches as 1/T. *)

prec = MachinePrecision; (* the DFT is exponentially accurate; 1e-10 is plenty *)
q = 1; c = 10; rho = 1; w0 = 5;
beta = rho w0/c;
nPts = 128;   (* samples per period for the DFT *)
mMax = 40;    (* harmonics kept: g_m falls like beta^|m| *)

(* position, velocity, acceleration on the circle; n = x *)
fy[t_] := Module[{v, a, n = {1, 0, 0}, b, bd, kap, vec},
  v = rho w0 {-Sin[w0 t], Cos[w0 t], 0};
  a = -rho w0^2 {Cos[w0 t], Sin[w0 t], 0};
  b = v/c; bd = a/c; kap = 1 - n.b;
  vec = Cross[n, Cross[n - b, bd]]/kap^2;
  vec[[2]]];

ts = Table[N[2 Pi k/(w0 nPts)], {k, 0, nPts - 1}];
fyTab = N[fy /@ ts];
cosTab = N[Cos[w0 ts]];
expTab = Table[Exp[-I m w0 ts], {m, -mMax, mMax}];

(* g_m(w) for m = -mMax..mMax *)
gm[w_?NumericQ] := Module[{g},
  g = fyTab Exp[-I w (rho/c) cosTab];
  expTab.g/nPts];

amp[w_?NumericQ, T_] := Module[{g = gm[w]},
  Total[Table[g[[m + mMax + 1]] (Exp[I w T] - 1)/(I (w + m w0)), {m, -mMax, mMax}]]];

band[lo_, hi_, nTurns_] := Module[{T = nTurns 2 Pi/w0},
  (* split at the line centre (the harmonic inside the band) and around it, where the
     sinc^2 of width 2 pi/T sits *)
  With[{mid = (lo + hi)/2, wd = 2 Pi/T},
    q^2/(4 Pi^2 c) NIntegrate[Abs[amp[w, T]]^2,
      {w, lo, mid - 20 wd, mid - wd, mid, mid + wd, mid + 20 wd, hi},
      PrecisionGoal -> 10, AccuracyGoal -> 20, MaxRecursion -> 40,
      Method -> {"GlobalAdaptive", "MaxErrorIncreases" -> 100000}]]];

jackson[m_, nTurns_] := Module[{T = nTurns 2 Pi/w0, jp},
  jp = (BesselJ[m - 1, m beta] - BesselJ[m + 1, m beta])/2;
  N[T q^2 m^2 w0^2 beta^2/(2 Pi c) jp^2]];

Print["// (turns, harmonic m, exact band energy (m -/+ 1/2) w0, infinite-time Jackson)"];
Do[
  Do[
    e = band[(m - 1/2) w0, (m + 1/2) w0, nT];
    Print["(", nT, ", ", m, ", ", ToString[N[e, 17], CForm], ", ",
      ToString[N[jackson[m, nT], 17], CForm], "),"],
    {m, 1, 4}],
  {nT, {40, 160}}];
