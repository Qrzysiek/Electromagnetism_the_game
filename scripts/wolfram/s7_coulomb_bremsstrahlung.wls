(* Test S7 (crates/physics/tests/spectra.rs), Jackson Pr. 15.10: the radiation of a charge
   deflected by a fixed repulsive Coulomb charge, on the exact hyperbolic orbit

     x = a (eps + cosh xi),  y = b sinh xi,  w0 t = xi + eps sinh xi,
     a = q Q/(m v^2),  eps = sqrt(1 + b^2/a^2),  w0 = v/a,  r = a (1 + eps cosh xi),

   (units k = 1; the force on the charge is q Q r/r^3, away from the fixed charge at the
   origin). Reference: Jackson's 14.65 in full (not the dipole approximation),
     d2I/dw dOmega = q^2/(4 pi^2 c) |Int n x ((n - beta) x betadot)/kappa^2 e^{i w (t - n.r/c)} dt|^2,
   for in-plane directions n = (cos th, sin th), on the Newtonian orbit (the test's flight
   is relativistic at v/c = 1e-3: the orbits differ at O(beta^2) = 1e-6). In the plane the
   vector is along e = (-sin th, cos th), with the component
     -(e.beta)(n.betadot) - (e.betadot) kappa.
   On the real xi axis the integrand oscillates as e^{i nu eps sinh xi} (nu = w/w0); it is
   analytic in the strip 0 <= Im xi <= pi/2 (the orbit's poles, 1 + eps cosh xi = 0, lie at
   |Im xi| = pi - arccos(1/eps) > pi/2) and decays at its ends, so the path is moved to
   Im xi = pi/2, where the phase becomes the factor e^{-nu eps cosh X}: the tails are cut
   where it is below e^{-60}. Checked against the real-axis integral below.
   Also: the dipole approximation (no retardation, kappa = 1) integrated over all angles,
   dI/dw = 2 q^2/(3 pi c^3) |Int a e^{i w t} dt|^2, against Jackson's closed form of Pr. 15.10(a),
     dI/dw = (8/(3 pi)) (q a w)^2/c^3 e^{-pi nu} {[K'_{i nu}(nu eps)]^2 + (eps^2 - 1)/eps^2 [K_{i nu}(nu eps)]^2},
   (the book writes y = -b sinh xi: the mirror image of this orbit, the same energy).
   Wolfram Engine 14.2, machine precision; run: wolframscript -file s7_coulomb_bremsstrahlung.wls *)

a = 1.; vinf = 1.; b = 2.; eps = Sqrt[1 + b^2/a^2]; w0 = vinf/a; c = 10000.; q = 1.;
r[x_] := {a (eps + Cosh[x]), b Sinh[x]};
tt[x_] := (x + eps Sinh[x])/w0;
dtdx[x_] := (1 + eps Cosh[x])/w0;
vel[x_] := {a Sinh[x], b Cosh[x]}/dtdx[x];
acc[x_] := a vinf^2 r[x]/(a (1 + eps Cosh[x]))^3;
check[v_, what_] := If[! NumericQ[v], Print["not numeric: ", what]; Abort[], v];

amp[w_, th_, x_] := Module[{n = {Cos[th], Sin[th]}, e = {-Sin[th], Cos[th]}, be, bd, k},
  be = vel[x]/c; bd = acc[x]/c; k = 1 - n.be;
  (-(e.be) (n.bd) - (e.bd) k)/k^2 Exp[I w (tt[x] - n.r[x]/c)] dtdx[x]];
xmax[w_] := ArcCosh[60/((w/w0) eps)];
full[w_, th_] := Module[{A = NIntegrate[amp[w, th, X + I Pi/2], {X, -xmax[w], 0, xmax[w]},
    MaxRecursion -> 20, PrecisionGoal -> 13, AccuracyGoal -> 40]},
  check[q^2/(4 Pi^2 c) Abs[A]^2, {w, th}]];

(* The contour against the real axis at w = w0 (there the tails beyond |xi| = 9, where the
   acceleration is 1e-8 and oscillates ever faster, are cut; the real-axis integral is itself
   accurate to about 1e-6). *)
real[w_, th_] := Module[{A = NIntegrate[amp[w, th, x], {x, -9, -6, -4, -2, 0, 2, 4, 6, 9},
    MaxRecursion -> 50, PrecisionGoal -> 10]},
  check[q^2/(4 Pi^2 c) Abs[A]^2, {"real", w, th}]];
Print["(* contour vs real axis: w, th, contour, real axis, relative difference *)"];
Do[With[{f = full[w, th], g = Quiet[real[w, th]]},
    Print[w, " ", th/Degree, " ", ScientificForm[f, 13], " ", ScientificForm[g, 13], " ",
      ScientificForm[g/f - 1, 3]]], {w, {1.}}, {th, {70 Degree, 150 Degree}}];

(* The dipole approximation over all angles, against Jackson's closed form. *)
dip[w_] := Module[{A = Table[NIntegrate[acc[X + I Pi/2][[i]] Exp[I w tt[X + I Pi/2]] dtdx[X + I Pi/2],
      {X, -xmax[w], 0, xmax[w]}, MaxRecursion -> 20, PrecisionGoal -> 13], {i, 2}]},
  check[2 q^2/(3 Pi c^3) (Abs[A[[1]]]^2 + Abs[A[[2]]]^2), {"dipole", w}]];
jackson[w_] := Module[{nu = w/w0, k, kp},
  k = BesselK[I nu, nu eps];
  kp = D[BesselK[I nu, z], z] /. z -> nu eps;
  check[Re[(8 q^2 a^2 w^2)/(3 Pi c^3) Exp[-Pi nu] (kp^2 + (eps^2 - 1)/eps^2 k^2)], {"jackson", w}]];
Print["(* dipole approximation dI/dw over all angles: w, numerical, Jackson's closed form, ratio *)"];
Do[With[{d = dip[w], j = jackson[w]},
    Print[w, " ", ScientificForm[d, 13], " ", ScientificForm[j, 13], " ", ScientificForm[j/d - 1, 3]]],
  {w, {0.25, 1., 3.}}];

(* The table for the test: (w, theta in degrees, d2I/dw dOmega), and the dipole
   approximation's value at the same point (its difference is O(beta) = 1e-3). *)
dipin[w_, th_] := Module[{e = {-Sin[th], Cos[th]}, A},
  A = NIntegrate[(e.acc[X + I Pi/2]) Exp[I w tt[X + I Pi/2]] dtdx[X + I Pi/2],
    {X, -xmax[w], 0, xmax[w]}, MaxRecursion -> 20, PrecisionGoal -> 13];
  check[q^2/(4 Pi^2 c^3) Abs[A]^2, {"dipin", w, th}]];
Print["(* (w, theta [deg], d2I/dw dOmega exact) for S7; then the dipole value's ratio *)"];
Do[With[{f = full[w, th Degree], d = dipin[w, th Degree]},
    Print["    (", w, ", ", th, ".0, ", ToString[N[f, 17], CForm], "), // dipole/exact - 1 = ",
      ScientificForm[d/f - 1, 3]]],
  {w, {0.25, 1., 3.}}, {th, {0, 70, 150, 260}}];
