#!/usr/bin/env wolframscript
(* Test B7 (PHYSICS.md §3.3, Jackson §12.6): three slow charges (q = 0.3, m = 1, in the
   plane) launched at (0, 0.5), (0, -0.5), (0, 0.1) with momenta (1, 0.2), (1, -0.1),
   (1, -0.3), flown for t = 10 at c = 50, 100, 200. The game's exact retarded interaction
   agrees with the Darwin interaction to O(1/c^2); the rest is O(1/c^3) (the mutual
   radiation fields). Here the Darwin dynamics (Landau & Lifshitz §65, Jackson §12.6), in
   the game's units (Coulomb constant 1, force q(E + v x B), B = v x E/c^2 for a slow source):
     E_j(x) = q_j [ n (1 + v_j^2/(2c^2) - 3 (n.v_j)^2/(2c^2)) / r^2
                    - (a_j + n (n.a_j)) / (2 c^2 r) ],
     B_j(x) = (q_j/c^2) v_j x n / r^2,       n = (x - x_j)/r,
     d(gamma_i m v_i)/dt = q_i Sum_{j != i} (E_j + v_i x B_j),
   with the sources' accelerations a_j inside the O(1/c^2) term taken as their Coulomb
   accelerations (the difference is O(1/c^4)). In the plane v x (v x n) stays in the
   plane. The initial velocities are the mechanical v = p/(gamma m) of the game's launch.
   Usage: wolframscript -file b7_darwin.wls [tEnd] [c...]  (defaults 10, 50 100 200;
   smoke test: 0.1 50). Output per c: the three end positions. *)

args = Rest[$ScriptCommandLine];
tEnd = If[Length[args] >= 1, ToExpression[args[[1]]], 10];
cs = If[Length[args] >= 2, ToExpression /@ Drop[args, 1], {50, 100, 200}];
q = 0.3; m = 1.0;
start = {{0, 0.5}, {0, -0.5}, {0, 0.1}};
mom = {{1., 0.2}, {1., -0.1}, {1., -0.3}};
cross2[a_, b_] := a[[1]] b[[2]] - a[[2]] b[[1]];          (* z of a x b *)
vxb[v_, bz_] := {v[[2]] bz, -v[[1]] bz};                     (* v x (bz z) *)

rhs[c_][s_?(VectorQ[#, NumericQ] &)] := Module[{x, p, g, v, a0, acc, f},
  x = Partition[s[[1 ;; 6]], 2]; p = Partition[s[[7 ;; 12]], 2];
  g = Sqrt[1 + (#.#)/(m^2 c^2)] & /@ p; v = p/(g m);
  (* Coulomb accelerations of the sources (for the O(1/c^2) term) *)
  a0 = Table[Sum[If[j == i, {0, 0}, q q (x[[i]] - x[[j]])/Norm[x[[i]] - x[[j]]]^3], {j, 3}]/m, {i, 3}];
  f = Table[Sum[If[j == i, {0, 0}, Module[{d = x[[i]] - x[[j]], r, n, e, bz},
       r = Sqrt[d.d]; n = d/r;
       e = q (n (1 + v[[j]].v[[j]]/(2 c^2) - 3 (n.v[[j]])^2/(2 c^2))/r^2
              - (a0[[j]] + n (n.a0[[j]]))/(2 c^2 r));
       bz = q cross2[v[[j]], n]/(c^2 r^2);
       q (e + vxb[v[[i]], bz])]], {j, 3}], {i, 3}];
  Join[Flatten[v], Flatten[f]]];

run[c_] := Module[{v0, sol, out},
  sol = NDSolve[{s'[t] == rhs[c][s[t]], s[0] == Join[Flatten[start], Flatten[mom]]}, s,
    {t, 0, tEnd}, PrecisionGoal -> 13, AccuracyGoal -> 14, MaxSteps -> 10^6,
    Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 8}];
  out = (s[tEnd] /. First[sol])[[1 ;; 6]];
  If[! VectorQ[out, NumericQ], Print["non-numeric result for c = ", c]; Abort[]];
  out];

Print["// t = ", tEnd, "; c, then x1, y1, x2, y2, x3, y3 (Darwin)"];
Do[
  Print["(", c, ".0, [", StringRiffle[ToString[N[#, 17], CForm] & /@ run[c], ", "], "]),"],
  {c, cs}];
