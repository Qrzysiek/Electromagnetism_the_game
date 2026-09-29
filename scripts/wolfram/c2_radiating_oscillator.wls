#!/usr/bin/env wolframscript
(* Test C2 (PHYSICS.md §2.1, §3.1; Jackson Pr. 16.1): a charge q = -1 (m = 1) inside a
   uniform cloud (Q = 1, radius 4, centred at the origin), c = 2, with Landau-Lifshitz
   radiation reaction, launched at (1, 0) with kinetic energy m v0^2/2, v0 = 0.6 omega0,
   along +y; flown to t = 5/Gamma, Gamma = omega0^2 tau.
   An independent integration of the full relativistic equation (Landau & Lifshitz,
   Classical Theory of Fields §76, in the game's units: Coulomb constant 1, force
   q(E + v x B)); with B = 0 and E = Q x/R^3 inside the cloud (DE/Dt = Q v/R^3):
     dp/dt = q E + (2 q^3/(3 m c^3)) gamma Q v/R^3
                 + (2 q^4/(3 m^2 c^5)) E (v.E)
                 - (2 q^4/(3 m^2 c^5)) gamma^2 v (E^2 - (E.v)^2/c^2),
     dx/dt = v = p/(gamma m).
   Output: the end state (x, y, px, py) at t = 5/Gamma, to compare with the game's. *)

prec = 30;
q = -1; m = 1; bigQ = 1; bigR = 4; c = 2;
omega0 = Sqrt[Abs[q] bigQ/(m bigR^3)];
tau = 2 q^2/(3 m c^3);
gam0 = omega0^2 tau;
tEnd = 5/gam0;
v0 = 6/10 omega0;
t0 = m v0^2/2;
p0 = Sqrt[(t0 + m c^2)^2 - m^2 c^4]/c;
ef[x_, y_] := bigQ {x, y}/bigR^3;
eqs = Module[{p = {px[t], py[t]}, xv = {x[t], y[t]}, g, v, e},
  g = Sqrt[1 + p.p/(m^2 c^2)]; v = p/(g m); e = ef @@ xv;
  Join[
    Thread[D[xv, t] == v],
    Thread[D[p, t] == q e + 2 q^3/(3 m c^3) g bigQ v/bigR^3 + 2 q^4/(3 m^2 c^5) e (v.e)
      - 2 q^4/(3 m^2 c^5) g^2 v (e.e - (e.v)^2/c^2)]]];
sol = NDSolve[Join[eqs, {x[0] == 1, y[0] == 0, px[0] == 0, py[0] == p0}],
  {x, y, px, py}, {t, 0, tEnd}, WorkingPrecision -> prec, PrecisionGoal -> 20,
  AccuracyGoal -> 22, MaxSteps -> 10^7, Method -> "ExplicitRungeKutta"];
end = {x[tEnd], y[tEnd], px[tEnd], py[tEnd]} /. First[sol];
Print["// t_end, then x, y, px, py at t_end"];
Print[ToString[N[tEnd, 17], CForm], ", ", StringRiffle[ToString[N[#, 17], CForm] & /@ end, ", "]];
