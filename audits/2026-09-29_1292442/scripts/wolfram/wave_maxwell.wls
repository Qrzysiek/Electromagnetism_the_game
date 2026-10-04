(* Vacuum Maxwell identities for the implemented plane wave, numbers first. *)
kx = 3/5; ky = 4/5; c = 4; w = 5; ex = -ky; ey = kx;
phase = w*(t - (kx*x + ky*y)/c);
ee = {ex, ey, 0}*Cos[phase];
bb = Cross[{kx, ky, 0}, ee]/c;
maxwell = {
  D[ee[[1]], x] + D[ee[[2]], y] + D[ee[[3]], z],
  {D[ee[[3]], y] - D[ee[[2]], z], D[ee[[1]], z] - D[ee[[3]], x], D[ee[[2]], x] - D[ee[[1]], y]} + D[bb, t],
  {D[bb[[3]], y] - D[bb[[2]], z], D[bb[[1]], z] - D[bb[[3]], x], D[bb[[2]], x] - D[bb[[1]], y]} - D[ee, t]/c^2
};
Print["symbolic norm: ", Norm[Flatten[Simplify[maxwell]]]];
Print["at a point: ", Norm[Flatten[N[maxwell /. {x -> 1, y -> 2, t -> 3, z -> 0}, 20]]]];
(* omega = 0, same B/E, constant fields: curls and time derivatives vanish *)
e0 = {0, 2, 0}; b0 = Cross[{1, 0, 0}, e0]/c;
Print["static |B| - |E|/c: ", Norm[b0] - Norm[e0]/c];
