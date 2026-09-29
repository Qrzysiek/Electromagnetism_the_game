(* Test R6 (crates/physics/tests/radiation_reaction.rs), Jackson Pr. 16.10-16.11: a charge
   crossing a gap of uniform field with radiation damping. The gap is two charged grids,
   smeared over a width w:
     E_x(x) = E0 g(x),  g(x) = (tanh(x/w) - tanh((x - d)/w))/2,
   crossed from x = -l to x = d + l (l = 20 w, where g < 5e-18), launched at speed v0.

   In one dimension the Lorentz-Dirac equation, written with the rapidity y (p = m c sinh y)
   and the proper time s, is exactly the Abraham-Lorentz equation (Pr. 16.8):
     m c (y' - tau y'') = f(s),  f = q E_x(x(s)),  x' = c sinh y,  t' = cosh y,
     tau = 2 q^2/(3 m c^3),
   whose physical (non-runaway) solution is the integro-differential form of Pr. 16.10(a),
     m c y'(s) = Int_0^Inf e^-u f(s + tau u) du,
   expanded in Pr. 16.10(b) as m c y' = f + tau f' + tau^2 f'' + ... (f' = df/ds). The first
   two terms are the Landau-Lifshitz equation: in one dimension, with E along the motion,
   LL is exactly m c y' = f + tau df/ds (dp/dt = qE + tau q gamma dE/dt).

   Computed here, each as the arrival time at x = d + l (t = 0 at x = -l) and the momentum
   there:
   - free: without radiation reaction;
   - LL: the Landau-Lifshitz flight (forward in s);
   - three terms: the series to tau^2 (a diagnostic: how much of LL's error the next term
     accounts for);
   - exact: integrated backward in s from x = d + l, where the solution is force-free
     (y' = f/(m c) there, 5e-18 of the gap's), which is stable (the runaway mode e^{s/tau}
     decays backward), shooting on the exit rapidity so that the rapidity at x = -l is
     y0 = artanh(v0/c); implicit Runge-Kutta (the backward problem is stiff), checked at
     another order.
   Cases: (A) c = 20, v0 = 1/4 -> 1/2 over d = 10 with w = 1/100: nonrelativistic
   (beta <= 0.025) with sharp edges (w/d = 1e-3), for Jackson's first-order formulas of
   Pr. 16.11(b); (B) c = 5, v0 = 1 -> 2.09 over d = 1 (beta 0.2 -> 0.42) with w = 1/5,
   1/10, 1/20: steep edges, where LL's error is measurable.
   Wolfram Engine 14.2, WorkingPrecision 32 (goals 20 digits); about 20 minutes;
   run: wolframscript -file r6_gap_damping.wls *)

q = 1; m = 1; wp = 32;
opts = Sequence[WorkingPrecision -> wp, PrecisionGoal -> 20, AccuracyGoal -> 20,
   MaxSteps -> 10^8];
tauOf[c_] := 2 q^2/(3 m c^3);
fOf[{c_, e0_, d_, w_, v0_}] := Function[xx, q e0 (Tanh[xx/w] - Tanh[(xx - d)/w])/2];
fpOf[{c_, e0_, d_, w_, v0_}] := Function[xx, q e0 (Sech[xx/w]^2 - Sech[(xx - d)/w]^2)/(2 w)];
fppOf[{c_, e0_, d_, w_, v0_}] := Function[xx,
   q e0 (-Sech[xx/w]^2 Tanh[xx/w] + Sech[(xx - d)/w]^2 Tanh[(xx - d)/w])/w^2];
check[v_, what_] := If[! VectorQ[v, NumericQ], Print["not numeric: ", what]; Abort[], v];

(* order 0: free; 1: Landau-Lifshitz; 2: three terms. {t at x = d + l, p there, y there}. *)
forward[cs : {c_, e0_, d_, w_, v0_}, order_] := Module[{sol, x, t, y, se, tau = tauOf[c],
    f = fOf[cs], fp = fpOf[cs], fpp = fppOf[cs], rhs, l = 20 w},
  rhs = Which[order == 0, f[x[s]]/(m c),
    order == 1, (f[x[s]] + tau fp[x[s]] c Sinh[y[s]])/(m c),
    True, (f[x[s]] + tau fp[x[s]] c Sinh[y[s]] + tau^2 fpp[x[s]] c^2 Sinh[y[s]]^2)/
      (m c - tau^2 fp[x[s]] c Cosh[y[s]])];
  sol = First@NDSolve[{x'[s] == c Sinh[y[s]], t'[s] == Cosh[y[s]], y'[s] == rhs,
      x[0] == -l, t[0] == 0, y[0] == ArcTanh[v0/c],
      WhenEvent[x[s] == d + l, "StopIntegration"]}, {x, t, y}, {s, 0, 10^4},
     Evaluate[opts], Method -> "Extrapolation"];
  se = (x /. sol)["Domain"][[1, 2]];
  check[{(t /. sol)[se], m c Sinh[(y /. sol)[se]], (y /. sol)[se]}, {cs, order}]];

(* The exact solution backward from x = d + l with exit rapidity yEnd: {y at x = -l, t_arr}. *)
back[cs : {c_, e0_, d_, w_, v0_}, yEnd_?NumericQ, meth_] := Module[{sol, x, t, y, z, ss,
    tau = tauOf[c], f = fOf[cs], l = 20 w},
  sol = First@NDSolve[{x'[s] == c Sinh[y[s]], t'[s] == Cosh[y[s]], y'[s] == z[s],
      tau z'[s] == z[s] - f[x[s]]/(m c),
      x[0] == d + l, t[0] == 0, y[0] == yEnd, z[0] == f[d + l]/(m c),
      WhenEvent[x[s] == -l, "StopIntegration"]}, {x, t, y, z}, {s, -10^4, 0},
     Evaluate[opts], Method -> meth];
  ss = (x /. sol)["Domain"][[1, 1]];
  check[{(y /. sol)[ss], -(t /. sol)[ss]}, {cs, yEnd}]];
yStart[cs_, yy_?NumericQ, meth_] := First[back[cs, yy, meth]];
exact[cs : {c_, e0_, d_, w_, v0_}, yGuess_, meth_] := Module[{yy, r, b},
  r = FindRoot[yStart[cs, yy, meth] == ArcTanh[v0/c], {yy, yGuess, yGuess + 10^-8},
    WorkingPrecision -> wp, AccuracyGoal -> wp - 8, PrecisionGoal -> wp - 8];
  b = back[cs, yy /. r, meth];
  {b[[2]], m c Sinh[yy /. r], yy /. r}];

cases = {{20, 3/320, 10, 1/100, 1/4}, {5, 2, 1, 1/5, 1}, {5, 2, 1, 1/10, 1}, {5, 2, 1, 1/20, 1}};
methods = {{"ImplicitRungeKutta", "DifferenceOrder" -> 10},
   {"ImplicitRungeKutta", "DifferenceOrder" -> 8}};
f17[v_] := ToString[N[v, 17], CForm];
rows = {};
sci[v_] := ToString[N[v, 3], CForm];
Do[Module[{fr, ll, t3, ex, ex2, tm},
    {tm, {fr, ll, t3}} = AbsoluteTiming[{forward[cs, 0], forward[cs, 1], forward[cs, 2]}];
    {tm, ex} = AbsoluteTiming[exact[cs, ll[[3]], methods[[1]]]];
    (* The check: the implicit method at another order (the backward integration is stiff:
       explicit methods stall, and "StiffnessSwitching" missed the goals by 1e-7). *)
    ex2 = exact[cs, ll[[3]], methods[[2]]];
    Print["(* ", cs, ": exact in ", tm, " s; second method differs by ",
      sci[(ex2[[1]] - ex[[1]])/(ll[[1]] - fr[[1]])], " (time), ",
      sci[(ex2[[2]] - ex[[2]])/(ll[[2]] - fr[[2]])], " (momentum) of the effect *)"];
    Print["(*   LL effect: dt = ", sci[ll[[1]] - fr[[1]]], ", dp = ", sci[ll[[2]] - fr[[2]]],
      "; exact - LL: ", sci[ex[[1]] - ll[[1]]], ", ", sci[ex[[2]] - ll[[2]]],
      "; three terms - LL: ", sci[t3[[1]] - ll[[1]]], ", ", sci[t3[[2]] - ll[[2]]], " *)"];
    AppendTo[rows, {cs, fr, ll, ex}]],
  {cs, cases}];

(* Jackson's first-order formulas (Pr. 16.11b) in case A, with T, v0, v1 from the free
   flight (extrapolating the uniform motions outside the gap to x = 0 and x = d). *)
Module[{cs = cases[[1]], fr, ll, ex, c, e0, d, w, v0, l, vel, v1, v1ll, v1ex, tt, ttll, ttex,
    tau, a},
  {fr, ll, ex} = rows[[1, 2 ;; 4]];
  {c, e0, d, w, v0} = cs; l = 20 w; tau = tauOf[c]; a = q e0/m;
  vel[p_] := p/Sqrt[m^2 + p^2/c^2];
  {v1, v1ll, v1ex} = vel /@ {fr[[2]], ll[[2]], ex[[2]]};
  {tt, ttll, ttex} = {fr[[1]] - l/v0 - l/v1, ll[[1]] - l/v0 - l/v1ll, ex[[1]] - l/v0 - l/v1ex};
  Print["(* A: T' - T: LL ", f17[ttll - tt], ", exact ", f17[ttex - tt], ", Jackson ",
    f17[-tau (1 - v0/v1)], " (exact/Jackson - 1 = ", sci[(ttex - tt)/(-tau (1 - v0/v1)) - 1], ") *)"];
  Print["(* A: v1' - v1: LL ", f17[v1ll - v1], ", exact ", f17[v1ex - v1], ", Jackson ",
    f17[-(a^2 tau/v1) tt], " (exact/Jackson - 1 = ", sci[(v1ex - v1)/(-(a^2 tau/v1) tt) - 1], ") *)"]];

Print["// (c, E0, d, w, v0, [t, p] free, [t, p] Landau-Lifshitz, [t, p] exact)"];
Do[Module[{cs, fr, ll, ex},
    {cs, fr, ll, ex} = r;
    Print["    (", Sequence @@ Riffle[f17 /@ N[cs, 20], ", "], ", [", f17[fr[[1]]], ", ",
      f17[fr[[2]]], "], [", f17[ll[[1]]], ", ", f17[ll[[2]]], "], [", f17[ex[[1]]], ", ",
      f17[ex[[2]]], "]),"]],
  {r, rows}];
