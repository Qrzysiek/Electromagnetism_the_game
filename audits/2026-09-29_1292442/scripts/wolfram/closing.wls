(* Independent closed forms for the checks the Python script evaluates numerically.
   Rectangular-strip and square-duct Laplace modes, and the circular-pipe root.
   Print is Null in wolframscript; every result goes through Print.
   One kernel. A non-numeric line is reported and the script continues. *)

show[name_, x_] := Module[{y},
  y = Quiet[Check[N[x, 17], $Failed]];
  If[NumericQ[y],
    Print[name, ": ", ToString[y, CForm]],
    Print[name, ": NOT NUMERIC"]
  ]
];

w = 2;
k = Pi/w;
phi = Sin[k x] Exp[-k z];
lap = Simplify[D[phi, {x, 2}] + D[phi, {z, 2}]];
Print["strip laplace: ", lap];
show["strip laplace at a point", lap /. {x -> 3/10, z -> 2/5}];

k3 = Pi Sqrt[2]/w;
phi3 = Sin[Pi x/w] Sin[Pi y/w] Exp[-k3 z];
lap3 = Simplify[D[phi3, {x, 2}] + D[phi3, {y, 2}] + D[phi3, {z, 2}]];
Print["square duct laplace: ", lap3];
show["game k over square-duct k, minus 1", (Pi/w)/k3 - 1];

j01 = BesselJZero[0, 1];
show["J0 at first zero", BesselJ[0, j01]];
show["pi over 2 j01, minus 1", Pi/(2 j01) - 1];

(* At rest the radiation field is transverse and the angular integral is Larmor.
   S = (c^2/4 pi) E x B with B = n x E / c and E perpendicular to n, so
   S.n R^2 = (c/4 pi) |E|^2 R^2 = q^2 a^2 Sin[theta]^2 / (4 pi c^3).
   Integral of Sin[theta]^3 over the sphere is 8 pi/3, and 8 pi/3 over 4 pi is 2/3. *)
ang = Integrate[Sin[th]^3, {th, 0, Pi}, {ph, 0, 2 Pi}];
show["angular integral minus 8 pi/3", ang - 8 Pi/3];
show["larmor factor minus 2/3", ang/(4 Pi) - 2/3];
