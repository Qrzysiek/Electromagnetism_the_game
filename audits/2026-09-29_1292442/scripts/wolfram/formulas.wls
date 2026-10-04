(* Independent checks of the engine formulas.
   Wolfram Engine is the reference. The script does not call the game.
   Regenerate: wolframscript -file audit/scripts/wolfram/formulas.wls
*)

Off[General::munfl];
$Assumptions = True;

show[name_, diff_] := Module[{d},
  d = Quiet[Check[N[diff, 17], $Failed]];
  If[!NumericQ[d],
    Print[name, ": NOT NUMERIC"];
    Return[]];
  d = Abs[d];
  Print[name, ": ", ToString[N[d, 17], CForm]];
];

(* ---------- units: CODATA 2018 k versus 1e-7 c^2 ---------- *)
cSI = 299792458;
kCode = 8.9875517923*^9;
mu0 = 1.25663706212*^-6; (* CODATA 2018 *)
kFromMu = mu0*cSI^2/(4*Pi);
kFromC = cSI^2*10^-7;
show["codata k vs code, relative", Abs[kFromMu - kCode]/kCode];
show["k vs 1e-7 c^2, relative", Abs[kFromC - kCode]/kCode];

(* ---------- kinetic energy identity ---------- *)
ke = FullSimplify[(Sqrt[1 + p2/(m^2 c^2)] - 1)*m*c^2 - p2/(m*(Sqrt[1 + p2/(m^2 c^2)] + 1)),
  m > 0 && c > 0 && p2 >= 0];
show["kinetic identity", Abs[ke]];

(* ---------- uniform-E world line (T6) ---------- *)
eps = Sqrt[m^2 c^4 + pp^2 c^2];
xPar = Integrate[q E0 t c^2/Sqrt[eps^2 + (q E0 c t)^2], {t, 0, T},
  Assumptions -> {q > 0, E0 > 0, c > 0, m > 0, pp >= 0, T > 0}];
xParClaim = (Sqrt[eps^2 + (q E0 c T)^2] - eps)/(q E0);
show["uniform E parallel", Abs[FullSimplify[xPar - xParClaim, q > 0 && E0 > 0 && c > 0 && T > 0]]];
yPer = Integrate[pp c^2/Sqrt[eps^2 + (q E0 c t)^2], {t, 0, T},
  Assumptions -> {q > 0, E0 > 0, c > 0, m > 0, pp > 0, T > 0}];
yClaim = (pp c/(q E0)) ArcSinh[q E0 c T/eps];
show["uniform E perpendicular", Abs[FullSimplify[yPer - yClaim, q > 0 && E0 > 0 && c > 0 && pp > 0 && T > 0]]];

(* ---------- relativistic Coulomb orbit ---------- *)
(* (W - kappa u)^2 - m^2 c^4 = L^2 c^2 (u'^2 + u^2); differentiate in phi. *)
orbit = FullSimplify[
  D[(W - kappa u[phi])^2 - m^2 c^4 - L^2 c^2 (u'[phi]^2 + u[phi]^2), phi]/(2 u'[phi])
   - (-kappa (W - kappa u[phi]) - L^2 c^2 (u''[phi] + u[phi])),
  L != 0 && c != 0];
show["orbit differentiation", Abs[orbit] /. {u[_] -> 1, u'[_] -> 1, u''[_] -> 1, W -> 1, kappa -> 1, L -> 1, c -> 1, m -> 1}];
(* u'' + Gamma^2 u + kappa W / (L^2 c^2) = 0 is the same statement as the differentiated energy. *)
lhs = uu + (1 - kappa^2/(L^2 c^2)) u;
rhs = -kappa W/(L^2 c^2);
brought = FullSimplify[(L^2 c^2) (lhs - rhs) + (-kappa (W - kappa u) - L^2 c^2 (uu + u)), L != 0 && c != 0];
show["orbit reduced form", Abs[brought]];
numOrbit = Quiet@Module[{mN = 1, cN = 20, v = 4, b = 3, kappaN = 2, g, WN, LN, G2, B0, A, uP, resid},
  g = 1/Sqrt[1 - v^2/cN^2];
  WN = g mN cN^2;
  LN = g mN v b;
  G2 = 1 - kappaN^2/(LN^2 cN^2);
  B0 = -kappaN WN/(LN^2 cN^2 G2);
  A = Sqrt[(WN^2 - mN^2 cN^4)/(LN^2 cN^2 G2) + B0^2];
  uP = A + B0;
  resid = (WN - kappaN uP)^2 - mN^2 cN^4 - LN^2 cN^2 uP^2;
  Abs[N[resid, 40]]/N[WN^2, 40]];
show["periapsis energy residual", numOrbit];

(* Rutherford limit of theta = Pi - 2/Gamma ArcCos[-B0/A], repulsive. *)
rutherford = Module[{mN = 1, v = 3/10, b = 2, kappaN = 1/2, cN, g, WN, LN, G2, B0, A, theta, claim},
  cN = 10^6;
  g = 1/Sqrt[1 - v^2/cN^2];
  WN = g mN cN^2;
  LN = g mN v b;
  G2 = 1 - kappaN^2/(LN^2 cN^2);
  B0 = -kappaN WN/(LN^2 cN^2 G2);
  A = Sqrt[(WN^2 - mN^2 cN^4)/(LN^2 cN^2 G2) + B0^2];
  theta = Abs[Pi - (2/Sqrt[G2]) ArcCos[-B0/A]];
  claim = 2 ArcTan[kappaN/(mN v^2 b)];
  Abs[theta - claim]/claim];
show["rutherford limit, relative", rutherford];

(* ---------- coil algebra ---------- *)
alpha2 = (a - rho)^2 + z^2;
beta2 = (a + rho)^2 + z^2;
mpar = 4 a rho/beta2;
gbr = beta2*((1 - mpar/2) EE - (1 - mpar) KK) - ((a^2 + rho^2 + z^2) EE - alpha2 KK);
show["coil B_rho bracket", Abs[FullSimplify[gbr, a > 0 && rho > 0 && Element[z, Reals]]]];
hbr = (4/Sqrt[mpar]) Sqrt[a/rho] - 2 Sqrt[beta2]/rho;
show["coil A_phi prefactor", Abs[FullSimplify[hbr, a > rho > 0 && z == 0]]];
axis = Module[{aa = a^2 + z^2, KK = Pi/2, EE = Pi/2},
  2 kappa/(aa Sqrt[aa]) ((a^2 - z^2) EE + aa KK) - 2 Pi kappa a^2/(a^2 + z^2)^(3/2)];
show["coil on axis", Abs[FullSimplify[axis, a > 0 && Element[z, Reals]]]];

(* AGM against EllipticK / EllipticE *)
agmKE[m_, m1_] := Module[{a = 1., b, c, sum, p2 = 0.5, k, n},
  b = Sqrt[m1]; c = Sqrt[m]; sum = 0.5 c c;
  Do[
    If[Abs[c] <= $MachineEpsilon a, Break[]];
    n = {0.5 (a + b), Sqrt[a b]};
    c = 0.5 (a - b); a = n[[1]]; b = n[[2]];
    p2 *= 2; sum += p2 c c;
  , {40}];
  k = Pi/(2 a); {k, k (1 - sum)}];
agmErr = Max@Table[Module[{m = SetPrecision[x, 32], m1 = SetPrecision[1 - x, 32], got, ref},
    got = agmKE[N[m, 20], N[m1, 20]];
    ref = {N[EllipticK[m], 20], N[EllipticE[m], 20]};
    Max[Abs[got - ref]/Abs[ref]]
  ], {x, {10^-8, 10^-3, 1/3, 0.9, 1 - 10^-8}}];
show["AGM vs EllipticK", agmErr];

gDirect = Series[(1 - m/2) EllipticE[m] - (1 - m) EllipticK[m], {m, 0, 4}];
show["g series m^2 coeff vs 3 Pi/32", Abs[SeriesCoefficient[gDirect, 2] - 3 Pi/32]];
hDirect = Series[(1 - m/2) EllipticK[m] - EllipticE[m], {m, 0, 4}];
show["h series m^2 coeff vs Pi/32", Abs[SeriesCoefficient[hDirect, 2] - Pi/32]];

(* ---------- straight segment: two logs, and curl A = B ---------- *)
len[v_] := Sqrt[v . v];
segRatio = Module[{obs = {x, y, z}, ra, rb, u = {1, 0, 0}, la, lb, pa, pb},
  ra = -obs; rb = {1, 0, 0} - obs;
  la = len[ra]; lb = len[rb]; pa = ra . u; pb = rb . u;
  FullSimplify[(lb + pb)/(la + pa) - (la - pa)/(lb - pb), y^2 + z^2 > 0]];
show["segment log forms", Abs[segRatio /. {x -> 1, y -> 1, z -> 1}]];
segCurl = Module[{obs = {x, y, z}, ra, rb, u = {1, 0, 0}, la, lb, pa, pb, ratio, AA, BB, cA, pt},
  ra = -obs; rb = {1, 0, 0} - obs;
  la = len[ra]; lb = len[rb]; pa = ra . u; pb = rb . u;
  ratio = (lb + pb)/(la + pa);
  AA = u Log[ratio];
  BB = Cross[ra, rb] (la + lb)/(la lb (la lb + ra . rb));
  cA = {D[AA[[3]], y] - D[AA[[2]], z], D[AA[[1]], z] - D[AA[[3]], x], D[AA[[2]], x] - D[AA[[1]], y]};
  pt = {x -> 2/5, y -> 3/10, z -> -1/7};
  Norm[N[(cA - BB) /. pt, 20]]];
show["segment curl A - B", segCurl];

(* ---------- Kelvin image on a grounded sphere ---------- *)
kelvin = Module[{a = 2, p = 5, q = 7, qi, pi, th, x, phi},
  qi = -q a/p; pi = a^2/p;
  th = 11/10;
  x = a {Cos[th], Sin[th], 0};
  phi = q/Norm[x - {p, 0, 0}] + qi/Norm[x - {pi, 0, 0}];
  Abs[FullSimplify[phi, a > 0]]];
show["kelvin sphere potential", kelvin];
(* image energy 1/2 q phi gives the force *)
imageForce = Module[{q = 1, a = 1, d, qi, U, fEnergy, fDirect},
  d = (p^2 - a^2)/p;
  qi = -q a/p;
  U = (1/2) q (qi/d);
  fEnergy = -D[U, p];
  fDirect = q qi/d^2;
  FullSimplify[fEnergy - fDirect, p > a > 0]];
show["image energy factor 1/2", Abs[N[imageForce /. p -> 3, 20]]];

(* ---------- Liénard a = 0 against the Heaviside field ---------- *)
heaviside = TimeConstrained[
  Module[{v, c = 5, q = 3, X, Y, R2, disc, tau, Rret, n, beta, kappa, eLW, Rnow, eH, diff},
    v = 2;
    X = 4; Y = 3;
    R2 = X^2 + Y^2;
    disc = X^2 v^2 + (c^2 - v^2) R2;
    tau = (X v - Sqrt[disc])/(v^2 - c^2); (* past root: numerator of the other sign is checked below *)
    (* pick the negative root explicitly *)
    tau = Module[{r1, r2},
      r1 = (X v + Sqrt[disc])/(v^2 - c^2);
      r2 = (X v - Sqrt[disc])/(v^2 - c^2);
      If[r1 < r2, r1, r2]];
    Rret = {X - v tau, Y, 0};
    n = Rret/Sqrt[Rret . Rret];
    beta = {v/c, 0, 0};
    kappa = 1 - n . beta;
    eLW = q (n - beta) (1 - beta . beta)/(kappa^3 (Rret . Rret));
    Rnow = {X, Y, 0};
    eH = q (1 - (v/c)^2) Rnow/(R2^(3/2) (1 - (v/c)^2 + (v X/(c Sqrt[R2]))^2)^(3/2));
    Norm[N[eLW - eH, 20]]
  ], 30, -1];
show["heaviside vs lienard", If[heaviside < 0, 1, heaviside]];

(* ---------- Landau–Lifshitz, Gaussian formula rewritten for F = q(E + v×B) ---------- *)
ll = Module[{
    e = {1, 2, -1}, b = {1/2, -1/3, 3}, v = {1/5, -1/7, 1/11},
    dE = {1, -2, 1/2}, dB = {0, 1, -1},
    g = 2, q = 3, m = 5, c = 4, bG, dBG, gaussian, si},
  (* Landau's Gaussian force, with B_G = c B and the same rule for the derivative. *)
  bG = c*b; dBG = c*dB;
  gaussian = (2 q^3/(3 m c^3))*g*(dE + Cross[v, dBG]/c) +
    (2 q^4/(3 m^2 c^4))*(Cross[e, bG] + Cross[bG, Cross[bG, v]]/c + e*(v . e)/c) -
    (2 q^4*g^2/(3 m^2 c^5))*v*((e + Cross[v, bG]/c) . (e + Cross[v, bG]/c) - (v . e)^2/c^2);
  si = (2 q^3/(3 m c^3))*g*(dE + Cross[v, dB]) +
    (2 q^4/(3 m^2 c^4))*(c*Cross[e, b] + c*Cross[b, Cross[b, v]] + e*(v . e)/c) -
    (2 q^4*g^2/(3 m^2 c^5))*v*((e + Cross[v, b]) . (e + Cross[v, b]) - (v . e)^2/c^2);
  Norm[gaussian - si]];
show["LL gaussian to these units", ll];

(* Uniform B, v perpendicular: the two velocity terms together equal the Lienard power,
   because 1 + gamma^2 beta^2 = gamma^2. *)
magPower = Module[{q = 3, m = 5, c = 7, v = 2, B = 4, g, a2, pLL, pLienard},
  g = 1/Sqrt[1 - v^2/c^2];
  pLL = (2 q^4/(3 m^2 c^3))*B^2*v^2 + (2 q^4/(3 m^2 c^5))*g^2*v^4*B^2;
  a2 = (q*v*B/(g*m))^2;
  pLienard = (2/3)*q^2*g^6*(a2 - v^2*a2/c^2)/c^3;
  Abs[pLL - pLienard]/pLienard];
show["LL magnetic power vs Lienard", magPower];

(* ---------- Poynting bound ---------- *)
(* |E × B| <= |E| |B|, and 2 |E| c |B| <= |E|^2 + c^2 |B|^2.
   The symbols stay global: Module would rename them and the substitution would miss. *)
poynt = FullSimplify[(e^2 + c^2*b^2)^2 - 4*c^2*e^2*b^2, e >= 0 && b >= 0 && c > 0];
show["poynting (u c)^2 - |S|^2 factor", (poynt /. {e -> 3, b -> 1, c -> 2}) - (3^2 - 2^2*1^2)^2];
(* equality case |E| = c |B| makes the factor 0 when they are perpendicular, already in the AM-GM *)
show["poynting equality when |E|=c|B|", Abs[FullSimplify[( (c b)^2 + c^2 b^2)^2 - 4 c^2 (c b)^2 b^2, b > 0 && c > 0]]];

(* ---------- plane wave, symbolic Maxwell ---------- *)
wave = Module[{kx = 3/5, ky = 4/5, c = 4, w = 5, ex, ey, phase, e, b},
  ex = -ky; ey = kx;
  phase = w (t - (kx x + ky y)/c);
  e = {ex, ey, 0} Cos[phase];
  b = Cross[{kx, ky, 0}, e]/c;
  Simplify[{
    D[e[[1]], x] + D[e[[2]], y] + D[e[[3]], z],
    {D[e[[3]], y] - D[e[[2]], z], D[e[[1]], z] - D[e[[3]], x], D[e[[2]], x] - D[e[[1]], y]} + D[b, t],
    {D[b[[3]], y] - D[b[[2]], z], D[b[[1]], z] - D[b[[3]], x], D[b[[2]], x] - D[b[[1]], y]} - D[e, t]/c^2
  }]
];
show["plane wave maxwell", Norm[Flatten[wave]]];

(* omega = 0 keeps B = k×E/c, and that pair is still a static solution *)
wave0 = Module[{c = 4, e = {0, 2, 0}, b},
  b = Cross[{1, 0, 0}, e]/c;
  (* curl of a constant is 0, time derivative is 0 *)
  Norm[b] - 2/c];
show["static wave |B|-|E|/c", Abs[wave0]];

(* ---------- dipole fields from the Lorenz-gauge potentials ---------- *)
dipole = Module[{p0, w = 17/10, c = 5/2, x0, t0, h = 10^-6, moment, phi, aa, ePot, bPot, eFor, bFor, tr, r, n, p, pd, pdd},
  p0 = {3/10, -1/5, 1/10};
  moment[tt_] := p0 Cos[w tt];
  phi[xx_, tt_] := Module[{rr = Sqrt[xx . xx], trr, nn},
    trr = tt - rr/c; nn = xx/rr;
    nn . moment[trr]/rr^2 + nn . moment'[trr]/(c rr)];
  aa[xx_, tt_] := Module[{rr = Sqrt[xx . xx], trr},
    trr = tt - rr/c; moment'[trr]/(c^2 rr)];
  x0 = {6/5, -2/5, 7/10}; t0 = 2/5;
  ePot = -Table[(phi[x0 + h UnitVector[3, i], t0] - phi[x0 - h UnitVector[3, i], t0])/(2 h), {i, 3}] -
    Table[(aa[x0, t0 + h][[i]] - aa[x0, t0 - h][[i]])/(2 h), {i, 3}];
  bPot = Table[
    Module[{dp = h UnitVector[3, i], dm = -h UnitVector[3, i]},
      0
    ], {i, 3}];
  bPot = {
    (aa[x0 + {0, h, 0}, t0][[3]] - aa[x0 - {0, h, 0}, t0][[3]])/(2 h) -
     (aa[x0 + {0, 0, h}, t0][[2]] - aa[x0 - {0, 0, h}, t0][[2]])/(2 h),
    (aa[x0 + {0, 0, h}, t0][[1]] - aa[x0 - {0, 0, h}, t0][[1]])/(2 h) -
     (aa[x0 + {h, 0, 0}, t0][[3]] - aa[x0 - {h, 0, 0}, t0][[3]])/(2 h),
    (aa[x0 + {h, 0, 0}, t0][[2]] - aa[x0 - {h, 0, 0}, t0][[2]])/(2 h) -
     (aa[x0 + {0, h, 0}, t0][[1]] - aa[x0 - {0, h, 0}, t0][[1]])/(2 h)
  };
  r = Sqrt[x0 . x0]; tr = t0 - r/c; n = x0/r;
  p = moment[tr]; pd = moment'[tr]; pdd = moment''[tr];
  eFor = (3 n (n . p) - p)/r^3 + (3 n (n . pd) - pd)/(c r^2) + Cross[n, Cross[n, pdd]]/(c^2 r);
  bFor = (Cross[pd, n]/r^2 + Cross[pdd, n]/(c r))/c^2;
  Norm[N[ePot - eFor, 12]] + Norm[N[bPot - bFor, 12]]
];
show["dipole from potentials", dipole];

(* static piece: E = -grad (n·p / r^2) *)
staticDip = Module[{p = {px, py, pz}, r, n, phiS, eS},
  r = Sqrt[{x, y, z} . {x, y, z}];
  n = {x, y, z}/r;
  phiS = n . p / r^2;
  eS = (3 n (n . p) - p)/r^3;
  FullSimplify[-Grad[phiS, {x, y, z}] - eS, x^2 + y^2 + z^2 > 0]
];
show["static dipole potential", Norm[staticDip]];

(* ---------- FieldMotion: Lambda^3 + omega^2 Lambda = 0, and the solution ---------- *)
fm = Module[{q, m, c, ex, ey, bz, k, L, w2, s, u0, u, du, lhs, poly},
  k = q/m;
  L = k {{0, ex/c, ey/c}, {ex/c, 0, bz}, {ey/c, -bz, 0}};
  w2 = k^2 (bz^2 - (ex^2 + ey^2)/c^2);
  poly = FullSimplify[L . L . L + w2 L];
  (* solution at a numeric point, derivative against Lambda.U *)
  {Norm[poly /. {q -> 1, m -> 2, c -> 3, ex -> 4, ey -> -1, bz -> 5}],
   Module[{qq = 1, mm = 2, cc = 3, Ex = 4, Ey = -1, Bz = 5, kk, LL, ww, ss = 1/4, U0, S1, C1, C0, UU, dU, rhs},
    kk = qq/mm;
    LL = kk {{0, Ex/cc, Ey/cc}, {Ex/cc, 0, Bz}, {Ey/cc, -Bz, 0}};
    ww = kk^2 (Bz^2 - (Ex^2 + Ey^2)/cc^2);
    U0 = {2, 1/2, -1/3};
    S1 = Sin[Sqrt[ww] ss]/Sqrt[ww];
    C1 = (1 - Cos[Sqrt[ww] ss])/ww;
    C0 = Cos[Sqrt[ww] ss];
    UU = U0 + S1 (LL . U0) + C1 (LL . LL . U0);
    dU = C0 (LL . U0) + S1 (LL . LL . U0);
    rhs = LL . UU;
    Norm[N[dU - rhs, 20]]
   ]}
];
show["field motion Lambda^3", fm[[1]]];
show["field motion dU/ds", fm[[2]]];

(* series of S1 at z = omega^2 s^2 = 1/2, the shader's switch, against Sin *)
fmSeries = Module[{z = 1/2, s = 1, w2, s1, exact},
  w2 = z/s^2;
  s1 = s (1 - z/6 (1 - z/20 (1 - z/42 (1 - z/72))));
  exact = Sin[Sqrt[w2] s]/Sqrt[w2];
  Abs[s1 - exact]/Abs[exact]];
show["field-motion series at shader switch", fmSeries];

(* ---------- Filon moment ---------- *)
filon = Module[{d = 7/5, exact, rec, e},
  exact = Table[Integrate[u^k Exp[I d u], {u, 0, 1}], {k, 0, 4}];
  e = Exp[I d];
  rec = { (e - 1)/(I d) };
  Do[AppendTo[rec, (e - k rec[[k]])/(I d)], {k, 1, 4}];
  Norm[N[exact - rec, 20]]
];
show["filon recurrence", filon];
filonSeries = Module[{d = 1/5, exact, acc, term},
  exact = Table[N[Integrate[u^k Exp[I d u], {u, 0, 1}], 24], {k, 0, 4}];
  acc = ConstantArray[0.`24, 5];
  term = 1.`24;
  Do[
    acc = acc + Table[term/(j + k + 1), {k, 0, 4}];
    term = term * (I d)/(j + 1);
  , {j, 0, 15}];
  Norm[exact - acc]
];
show["filon series", filonSeries];

(* ---------- abrupt-stop boundary term is minus the end value of F ---------- *)
(* integrand of 14.65 is d/dt of n×(n×beta)/kappa, up to the charge factor. Check d/dt identity for collinear motion. *)
stopId = Module[{c = 3, n = {1, 0, 0}, vx, ax, beta, k, f, df, rad},
  beta = {vx, 0, 0}/c;
  k = 1 - n . beta;
  f = Cross[n, Cross[n, beta]]/k;
  df = D[f, vx] ax; (* d/dt through vx only, n fixed: the far-zone direction *)
  rad = Cross[n, Cross[n - beta, {ax, 0, 0}/c]]/k^2;
  (* Jackson: d/dt [n×(n×beta)/kappa] = n×((n-beta)×beta_dot) / kappa^2 *)
  FullSimplify[df - rad, c > vx > 0]
];
show["sudden-stop integrand", Norm[stopId]];

(* ---------- quintic smoothstep used to blend the two continuations ---------- *)
ss[u_] := u^3 (10 - 15 u + 6 u^2);
show["smoothstep ends", Abs[ss[0]] + Abs[ss[1] - 1] + Abs[ss'[0]] + Abs[ss'[1]] + Abs[ss''[0]] + Abs[ss''[1]]];

(* ---------- shader Log[2] constant ---------- *)
show["shader ln2 constant", Abs[0.6931472 - Log[2]]];

(* ---------- Larmor of an oscillating dipole ---------- *)
larmor = Module[{p0, w, c, pdd, inst, avg},
  p0 = {px, py, pz}; w = omega; c = cc;
  pdd = -w^2 p0 Cos[w t];
  inst = (2/3) (pdd . pdd)/c^3;
  avg = Integrate[inst, {t, 0, 2 Pi/w}]/(2 Pi/w);
  FullSimplify[avg - (p0 . p0) w^4/(3 c^3), w > 0 && cc > 0]
];
show["dipole larmor average", Abs[larmor /. {px -> 1, py -> 0, pz -> 0, omega -> 1, cc -> 1}]];

(* Gradients of Bz used by the moment force, against a numerical derivative of Bz. *)
len[v_] := Sqrt[v . v];
segB[a_, b_, x_, kappa_] := Module[{ra = a - x, rb = b - x, la, lb},
  la = len[ra]; lb = len[rb];
  Cross[ra, rb]*kappa*(la + lb)/(la*lb*(la*lb + ra . rb))];
segG[a_, b_, x_, kappa_] := Module[{u = a - x, v = b - x, lu, lv, f},
  lu = len[u]; lv = len[v];
  f = kappa*(lu + lv)/(lu*lv*(lu*lv + u . v));
  -Cross[a - b, {0, 0, 1}]*f];
verts = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
Bsq[x_] := Sum[segB[verts[[i]], verts[[Mod[i, 4] + 1]], x, 1], {i, 4}];
Gsq[x_] := Sum[segG[verts[[i]], verts[[Mod[i, 4] + 1]], x, 1], {i, 4}];
ptG = {3/10, 1/5, 0}; hG = 10^-6;
numG = {
  (Bsq[ptG + {hG, 0, 0}][[3]] - Bsq[ptG - {hG, 0, 0}][[3]])/(2 hG),
  (Bsq[ptG + {0, hG, 0}][[3]] - Bsq[ptG - {0, hG, 0}][[3]])/(2 hG),
  (Bsq[ptG + {0, 0, hG}][[3]] - Bsq[ptG - {0, 0, hG}][[3]])/(2 hG)};
show["polygon grad Bz", Norm[Gsq[ptG] - numG]];
loopG = Module[{a = 2, rho = 1, kappa = 1, h = 10^-6, alpha2, beta2, beta, m, m1, g, claim, numd, bz},
  alpha2 = (a - rho)^2; beta2 = (a + rho)^2; beta = Sqrt[beta2];
  m = 4 a rho/beta2; m1 = alpha2/beta2;
  g = (1 - m/2) EllipticE[m] - m1 EllipticK[m];
  claim = 2 kappa beta g/(alpha2 rho);
  bz[rr_] := Module[{al = (a - rr)^2, be = (a + rr)^2, mm, KK, EE},
    mm = 4 a rr/be; KK = EllipticK[mm]; EE = EllipticE[mm];
    2 kappa/(al Sqrt[be])*((a^2 - rr^2) EE + al KK)];
  numd = (bz[rho + h] - bz[rho - h])/(2 h);
  Abs[claim - numd]];
show["loop grad Bz", loopG];

(* Dipole ∇B_z from the written formula, against Grad of B_z. *)
mvec = {1/2, -1/3, 4/5};
bzSym = 3*(mvec . {x, y, z})*z/(x^2 + y^2 + z^2)^(5/2) - mvec[[3]]/(x^2 + y^2 + z^2)^(3/2);
gSym = Grad[bzSym, {x, y, z}];
rd = {2/5, -1/7, 3/11};
gNumD = gSym /. {x -> rd[[1]], y -> rd[[2]], z -> rd[[3]]};
mrd = mvec . rd; r2d = rd . rd;
gFormD = 3*(mvec*rd[[3]] + {0, 0, mrd})/r2d^(5/2) - 15*mrd*rd[[3]]*rd/r2d^(7/2) + 3*mvec[[3]]*rd/r2d^(5/2);
show["dipole grad Bz", Norm[gFormD - gNumD]];

(* Elastic impulse: dT/dp = v, so the energy balance has slope -closing at J = 0.
   Equal masses, target at rest, head-on: the positive root is the incoming momentum.
   At large c that root approaches 2 μ (v1 - v2). *)
elasSlope = Module[{m1 = 3, m2 = 5, c = 4, p1 = 2, p2 = -1, T, v, f},
  T[m_, p_] := Sqrt[m^2*c^4 + p^2*c^2] - m*c^2;
  v[m_, p_] := p*c^2/Sqrt[m^2*c^4 + p^2*c^2];
  f[j_] := T[m1, p1 - j] + T[m2, p2 + j] - T[m1, p1] - T[m2, p2];
  Abs[f[0]] + Abs[f'[0] + (v[m1, p1] - v[m2, p2])]];
show["elastic slope at J=0", elasSlope];
elasExchange = Module[{m = 5, c = 3, p = 7, T, f},
  T[pp_] := Sqrt[m^2*c^4 + pp^2*c^2] - m*c^2;
  f[j_] := T[p - j] + T[j] - T[p];
  Abs[f[p]]];
show["elastic equal-mass exchange", elasExchange];
elasNewton = Module[{m1 = 3, m2 = 5, c = 10^6, p1 = 2, p2 = 0, T, v, f, newt, root},
  T[m_, p_] := Sqrt[m^2*c^4 + p^2*c^2] - m*c^2;
  v[m_, p_] := p*c^2/Sqrt[m^2*c^4 + p^2*c^2];
  f[j_] := T[m1, p1 - j] + T[m2, p2 + j] - T[m1, p1] - T[m2, p2];
  newt = 2*m1*m2/(m1 + m2)*(v[m1, p1] - v[m2, p2]);
  root = j /. FindRoot[f[j] == 0, {j, newt}];
  Abs[root - newt]/newt];
show["elastic root vs Newtonian, c=1e6", elasNewton];

Print["formula checks finished"];
