(* Unit-consistent Landau–Lifshitz substitution, and the plane wave without localized symbols. *)

show[name_, diff_] := Module[{d},
  d = Quiet[Check[N[diff, 20], $Failed]];
  If[!NumericQ[d], Print[name, ": NOT NUMERIC ", Short[diff, 2]]; Return[]];
  Print[name, ": ", ToString[Abs[d], CForm]];
];

ll = Module[{
    e = {1, 2, -1}, b = {1/2, -1/3, 3}, v = {1/5, -1/7, 1/11},
    dE = {1, -2, 1/2}, dB = {0, 1, -1},
    g = 2, q = 3, m = 5, c = 4, bG, dBG, gaussian, si},
  bG = c*b;
  dBG = c*dB;
  gaussian = (2 q^3/(3 m c^3))*g*(dE + Cross[v, dBG]/c) +
    (2 q^4/(3 m^2 c^4))*(Cross[e, bG] + Cross[bG, Cross[bG, v]]/c + e*(v . e)/c) -
    (2 q^4*g^2/(3 m^2 c^5))*v*((e + Cross[v, bG]/c) . (e + Cross[v, bG]/c) - (v . e)^2/c^2);
  si = (2 q^3/(3 m c^3))*g*(dE + Cross[v, dB]) +
    (2 q^4/(3 m^2 c^4))*(c*Cross[e, b] + c*Cross[b, Cross[b, v]] + e*(v . e)/c) -
    (2 q^4*g^2/(3 m^2 c^5))*v*((e + Cross[v, b]) . (e + Cross[v, b]) - (v . e)^2/c^2);
  Norm[gaussian - si]];
show["LL with B_G = c B", ll];

phase = w*(t - (kx*x + ky*y)/c);
ee = {ex, ey, 0}*Cos[phase];
bb = Cross[{kx, ky, 0}, ee]/c;
maxwell = {
  D[ee[[1]], x] + D[ee[[2]], y] + D[ee[[3]], z],
  {D[ee[[3]], y] - D[ee[[2]], z], D[ee[[1]], z] - D[ee[[3]], x], D[ee[[2]], x] - D[ee[[1]], y]} + D[bb, t],
  {D[bb[[3]], y] - D[bb[[2]], z], D[bb[[1]], z] - D[bb[[3]], x], D[bb[[2]], x] - D[bb[[1]], y]} - D[ee, t]/c^2
};
vals = {ex -> -ky, ey -> kx, kx -> 3/5, ky -> 4/5, c -> 4, w -> 5, x -> 1, y -> 2, t -> 3, z -> 0};
show["plane wave maxwell", Norm[Flatten[maxwell /. vals]]];

(* omega = 0 is still a vacuum solution: constant E and B = k×E/c *)
e0 = {0, 2, 0};
b0 = Cross[{1, 0, 0}, e0]/4;
show["static wave |B|", Norm[b0] - 2/4];

Print["rest2 finished"];
