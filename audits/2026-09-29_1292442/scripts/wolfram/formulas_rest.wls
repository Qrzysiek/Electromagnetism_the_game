(* Superseded by formulas.wls. The Landau–Lifshitz line in this file used one dB
   for both unit systems and printed 0.024; that number is a script error.
   With B_G = c B and dB_G = c dB the residual is 0 (formulas.wls). *)

show[name_, diff_] := Module[{d},
  d = Quiet[Check[N[diff, 20], $Failed]];
  If[!NumericQ[d], Print[name, ": NOT NUMERIC"]; Return[]];
  Print[name, ": ", ToString[Abs[d], CForm]];
];

(* Periapsis: energy balance at u = A + B0, arbitrary precision. *)
peri = Module[{m = 1, c = 20, v = 4, b = 3, kappa = 2, g, W, L, G2, B0, A, u, resid},
  g = 1/Sqrt[1 - v^2/c^2];
  W = g m c^2;
  L = g m v b;
  G2 = 1 - kappa^2/(L^2 c^2);
  B0 = -kappa W/(L^2 c^2 G2);
  A = Sqrt[(W^2 - m^2 c^4)/(L^2 c^2 G2) + B0^2];
  u = A + B0;
  resid = (W - kappa u)^2 - m^2 c^4 - L^2 c^2 u^2;
  N[Abs[resid]/W^2, 20]];
show["periapsis energy residual", peri];

(* Landau–Lifshitz: Gaussian expression with B_G = c B, against the code. *)
ll = Module[{
    e = {1, 2, -1}, b = {1/2, -1/3, 3}, v = {1/5, -1/7, 1/11},
    dE = {1, -2, 1/2}, dB = {0, 1, -1},
    g = 2, q = 3, m = 5, c = 4, bG, gaussian, si},
  bG = c b;
  gaussian = (2 q^3/(3 m c^3)) g (dE + Cross[v, dB]/c) +
    (2 q^4/(3 m^2 c^4)) (Cross[e, bG] + Cross[bG, Cross[bG, v]]/c + e (v . e)/c) -
    (2 q^4 g^2/(3 m^2 c^5)) v ((e + Cross[v, bG]/c) . (e + Cross[v, bG]/c) - (v . e)^2/c^2);
  si = (2 q^3/(3 m c^3)) g (dE + Cross[v, dB]) +
    (2 q^4/(3 m^2 c^4)) (c Cross[e, b] + c Cross[b, Cross[b, v]] + e (v . e)/c) -
    (2 q^4 g^2/(3 m^2 c^5)) v ((e + Cross[v, b]) . (e + Cross[v, b]) - (v . e)^2/c^2);
  Norm[gaussian - si]];
show["LL gaussian to these units", ll];

mag = Module[{q = 3, m = 5, c = 7, v = 2, B = 4, g, a2, pLL, pLienard},
  g = 1/Sqrt[1 - v^2/c^2];
  pLL = (2 q^4/(3 m^2 c^3)) B^2 v^2 + (2 q^4/(3 m^2 c^5)) g^2 v^4 B^2;
  a2 = (q v B/(g m))^2;
  pLienard = (2/3) q^2 g^6 (a2 - v^2 a2/c^2)/c^3;
  Abs[N[pLL - pLienard, 20]]/N[pLienard, 20]];
show["LL magnetic power vs Lienard", mag];

(* Symbolic: 1 + gamma^2 beta^2 = gamma^2 *)
show["gamma identity", FullSimplify[1 + (1/(1 - beta^2)) beta^2 - 1/(1 - beta^2), 0 < beta < 1]];

poynt = FullSimplify[(e^2 + c^2 b^2)^2 - 4 c^2 e^2 b^2, e >= 0 && b >= 0 && c > 0];
show["poynting factor at numbers", (poynt /. {e -> 3, b -> 1, c -> 2}) - (3^2 - 4)^2];

wave = Module[{c, w, kx, ky, ex, ey, phase, ee, bb},
  phase = w (t - (kx x + ky y)/c);
  ee = {ex, ey, 0} Cos[phase];
  bb = Cross[{kx, ky, 0}, ee]/c;
  {D[ee[[1]], x] + D[ee[[2]], y],
   {D[ee[[3]], y] - D[ee[[2]], z], D[ee[[1]], z] - D[ee[[3]], x], D[ee[[2]], x] - D[ee[[1]], y]} + D[bb, t],
   {D[bb[[3]], y] - D[bb[[2]], z], D[bb[[1]], z] - D[bb[[3]], x], D[bb[[2]], x] - D[bb[[1]], y]} - D[ee, t]/c^2}
];
waveN = Flatten[wave /. {ex -> -ky, ey -> kx, kx -> 3/5, ky -> 4/5, c -> 4, w -> 5, x -> 1, y -> 2, t -> 3, z -> 0}];
show["plane wave maxwell", Norm[N[waveN, 20]]];

Print["rest finished"];
