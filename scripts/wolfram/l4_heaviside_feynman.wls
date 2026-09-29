(* Test L4 (crates/physics/tests/lienard.rs), Jackson Pr. 6.2: the fields of a point charge
   in the forms of Heaviside and Feynman, against the Lienard-Wiechert fields. Units of the
   game: k = 1/(4 pi eps0) = 1, mu0/(4 pi) = 1/c^2. With R = x - r(t'), R^ = R/|R|, kappa =
   1 - v.R^/c, everything evaluated at the retarded time t' (c (t - t') = |R|), and d/dt the
   derivative in the observation time at fixed x, d/dt = (1/kappa) d/dt':
     (b) from Jefimenko's equations:
         E = q { R^/(kappa R^2) + d/dt[R^/(kappa R)]/c - d/dt[v/(kappa R)]/c^2 },
         B = (q/c^2) { v x R^/(kappa R^2) + d/dt[v x R^/(kappa R)]/c };
     (c) Feynman's electric field and Heaviside's magnetic field:
         E = q { R^/R^2 + (R/c) d/dt[R^/R^2] + d^2/dt^2[R^]/c^2 },
         B = (q/c^2) { v x R^/(kappa^2 R^2) + d/dt[v x R^/kappa]/(c R) };
   and the Lienard-Wiechert form,
         E = q { (R^ - beta)(1 - beta^2)/(kappa^3 R^2) + R^ x ((R^ - beta) x a/c)/(c kappa^3 R) },
         B = R^ x E/c.
   The world line: r(t') = (cos(6t'/5), sin(6t'/5), (3/10) sin(21t'/10)), c = 2 (beta up to
   0.68: relativistic, accelerated in three dimensions, with a varying jerk). The time
   derivatives are symbolic; the retarded time is found to 50 digits. Printed: the three
   forms' largest relative differences (they must agree to the working precision), then the
   table for the test (Feynman's E, Heaviside's B).
   Wolfram Engine 14.2; run: wolframscript -file l4_heaviside_feynman.wls *)

q = 1; c = 2;
rr[tp_] := {Cos[6 tp/5], Sin[6 tp/5], 3/10 Sin[21 tp/10]};
check[v_, what_] := If[! VectorQ[Flatten[v], NumericQ], Print["not numeric: ", what]; Abort[], v];

forms[x_, t_] := Module[{tr, rv, rn, n, v, a, kappa, dt, eF, eJ, eLW, bH, bJ, bLW, trv},
  rv = x - rr[tr]; rn = Sqrt[rv.rv]; n = rv/rn;
  v = D[rr[tr], tr]; a = D[v, tr];
  kappa = 1 - n.v/c;
  dt[f_] := D[f, tr]/kappa;
  eJ = q (n/(kappa rn^2) + dt[n/(kappa rn)]/c - dt[v/(kappa rn)]/c^2);
  bJ = (q/c^2) (Cross[v, n]/(kappa rn^2) + dt[Cross[v, n]/(kappa rn)]/c);
  eF = q (n/rn^2 + (rn/c) dt[n/rn^2] + dt[dt[n]]/c^2);
  bH = (q/c^2) (Cross[v, n]/(kappa^2 rn^2) + dt[Cross[v, n]/kappa]/(c rn));
  eLW = q ((n - v/c) (1 - v.v/c^2)/(kappa^3 rn^2) + Cross[n, Cross[n - v/c, a/c]]/(c kappa^3 rn));
  bLW = Cross[n, eLW]/c;
  trv = tr /. FindRoot[c (t - tr) == rn, {tr, t - Sqrt[x.x]/c},
     WorkingPrecision -> 60, AccuracyGoal -> 50, PrecisionGoal -> 50];
  check[N[{eF, eJ, eLW, bH, bJ, bLW} /. tr -> trv, 40], {x, t}]];

points = {{{3, 1/2, 1/5}, 11/10}, {{-6, 4, -2}, 11/10}, {{30, -25, 10}, 11/10},
   {{1/2, -3/2, 1}, 7/3}, {{0, 0, 50}, 7/3}, {{-2, -1, 0}, 13/4}};
f17[v_] := ToString[N[v, 17], CForm];
worst = 0;
rows = Table[Module[{eF, eJ, eLW, bH, bJ, bLW, de, db},
    {eF, eJ, eLW, bH, bJ, bLW} = forms[p[[1]], p[[2]]];
    de = Max[Norm[eJ - eF], Norm[eLW - eF]]/Norm[eF];
    db = Max[Norm[bJ - bH], Norm[bLW - bH]]/Norm[bH];
    worst = Max[worst, de, db];
    {p, eF, bH}], {p, points}];
Print["(* largest relative difference between the three forms: ", f17[worst], " *)"];
Print["// (x, t, E (Feynman), B (Heaviside))"];
Do[Print["    ([", Sequence @@ Riffle[f17 /@ N[r[[1, 1]], 20], ", "], "], ", f17[N[r[[1, 2]], 20]],
    ", [", Sequence @@ Riffle[f17 /@ r[[2]], ", "], "], [", Sequence @@ Riffle[f17 /@ r[[3]], ", "], "]),"],
  {r, rows}];
