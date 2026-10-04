(* Polygon and loop gradients of Bz, against a numerical derivative of the same B. *)
len[v_] := Sqrt[v . v];
segB[a_, b_, x_, kappa_] := Module[{ra = a - x, rb = b - x, la, lb},
  la = len[ra]; lb = len[rb];
  Cross[ra, rb]*kappa*(la + lb)/(la*lb*(la*lb + ra . rb))];
segG[a_, b_, x_, kappa_] := Module[{u = a - x, v = b - x, lu, lv, f},
  lu = len[u]; lv = len[v];
  f = kappa*(lu + lv)/(lu*lv*(lu*lv + u . v));
  -Cross[a - b, {0, 0, 1}]*f];

verts = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
Btot[x_] := Sum[segB[verts[[i]], verts[[Mod[i, 4] + 1]], x, 1], {i, 4}];
Gtot[x_] := Sum[segG[verts[[i]], verts[[Mod[i, 4] + 1]], x, 1], {i, 4}];
pt = {3/10, 1/5, 0};
h = 10^-6;
num = {(Btot[pt + {h, 0, 0}][[3]] - Btot[pt - {h, 0, 0}][[3]])/(2 h),
       (Btot[pt + {0, h, 0}][[3]] - Btot[pt - {0, h, 0}][[3]])/(2 h),
       (Btot[pt + {0, 0, h}][[3]] - Btot[pt - {0, 0, h}][[3]])/(2 h)};
Print["polygon grad vs dBz: ", Norm[N[Gtot[pt] - num, 12]]];

(* Circular loop at z = 0: dBz/drho from the elliptic formula, against the g(m) expression. *)
loop = Module[{a = 2, rho = 1, kappa = 1, z = 0, alpha2, beta2, beta, m, m1, K, E, bz, g, claim, numd},
  alpha2 = (a - rho)^2 + z^2;
  beta2 = (a + rho)^2 + z^2;
  beta = Sqrt[beta2];
  m = 4 a rho/beta2;
  m1 = alpha2/beta2;
  K = EllipticK[m]; E = EllipticE[m];
  g = (1 - m/2) E - (1 - m) K;
  claim = 2 kappa beta g/(alpha2 rho);
  bz[rr_] := Module[{al, be, mm, KK, EE},
    al = (a - rr)^2; be = (a + rr)^2; mm = 4 a rr/be;
    KK = EllipticK[mm]; EE = EllipticE[mm];
    2 kappa/(al Sqrt[be])*((a^2 - rr^2) EE + al KK)];
  numd = (bz[rho + h] - bz[rho - h])/(2 h);
  Abs[N[claim - numd, 12]]];
Print["loop dBz/drho: ", loop];
