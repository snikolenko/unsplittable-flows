"""Independent verification of the new zero-sum (B_zm) material."""
from fractions import Fraction as F
from itertools import combinations
import sympy as sp
from sympy.solvers.simplex import lpmin

print("=" * 74)
print("[A] Theorem 'four-item zero-sum wall' (q=3, k=4)")
f = [F(13, 64), F(781, 3136), F(481, 1024), F(3919, 50176)]
d = [F(15, 16), F(975, 1024), F(61, 64), F(1)]
pi = [0, 1, 2, 3]
sig = [2, 3, 0, 1]
k, q = 4, 3
x = [1 - fi for fi in f]
print(f"    sum f = {sum(f)}  (must be k-q = {k-q}): {sum(f) == k-q}")
assert sum(f) == k - q
print(f"    sum x = {sum(x)}  (must be q = {q}): {sum(x) == q}")
print(f"    max d = {max(d)}: {max(d) == 1}")
dfm = sum(d[i]*f[i] for i in range(k))
print(f"    d.f = {dfm}   paper says 3905/4096: {dfm == F(3905,4096)}")
assert dfm == F(3905, 4096)

PREF = [pi[:j] for j in range(1, k+1)] + [sig[:j] for j in range(1, k+1)]
print(f"    ordered prefix rows: {len(PREF)} (paper says eight)")


def A_O(O):
    return sum(d[i]*x[i] for i in range(k)) - sum(d[i] for i in O)


def solve(O, moment):
    """exact LP value; moment=True adds sum_{i in S} d_i c_i = A_O"""
    S = [i for i in range(k) if i not in O]
    cv = sp.symbols(f'c0:{k}')
    t = sp.Symbol('t')
    cons = [cv[i] >= 0 for i in S] + [sp.Eq(sum(cv[i] for i in S), 1)]
    if moment:
        cons.append(sp.Eq(sum(sp.Rational(d[i])*cv[i] for i in S),
                          sp.Rational(A_O(O))))
    for P in PREF:
        e = 0
        for i in P:
            ind = 0 if i in O else 1
            e = e + sp.Rational(d[i])*(ind - (cv[i] if i in S else 0)
                                       - sp.Rational(f[i]))
        cons.append(e <= t)
    try:
        val, _ = lpmin(t, cons)
    except Exception:
        return None
    v = sp.nsimplify(val)
    return F(int(sp.fraction(v)[0]), int(sp.fraction(v)[1]))


print("\n    per-pair analysis (moment-admissible iff A_O in [min_S d, max_S d]):")
adm, vals = [], {}
for O in combinations(range(k), 2):
    S = [i for i in range(k) if i not in O]
    a = A_O(O)
    lo, hi = min(d[i] for i in S), max(d[i] for i in S)
    ok = lo <= a <= hi
    v = solve(O, True) if ok else None
    if ok:
        adm.append(O)
        vals[O] = v
    print(f"      O={O}  A_O={str(a):>18}  in [{lo},{hi}]: {str(ok):<5} "
          f"B_zm(O) = {v}")
sel = [tuple(i for i in range(k) if i not in O) for O in adm]
print(f"    admissible OMITTED pairs : {adm}")
print(f"    equivalently SELECTED S  : {sorted(sel)}   "
      f"(the paper names the selected pairs: (0,3),(1,3),(2,3))")
assert set(sel) == {(0, 3), (1, 3), (2, 3)}
radii = sorted(set(vals.values()))
print(f"    the three radii: {[str(r) for r in radii]}")
print(f"    paper lists 2085/4096, 229515/458752, 98393/196608")
assert set(vals.values()) == {F(2085, 4096), F(229515, 458752), F(98393, 196608)}
bzm = min(vals.values())
argmin = [O for O in vals if vals[O] == bzm]
selmin = [tuple(i for i in range(k) if i not in O) for O in argmin]
print(f"    B_zm = min = {bzm} = {float(bzm):.10f}")
print(f"    attained at omitted {argmin}, i.e. selected {selmin}  "
      f"(paper says the selected pair (1,3))")
assert bzm == F(229515, 458752) and selmin == [(1, 3)]
print(f"    = 1/2 + {bzm - F(1,2)}   paper says 139/458752: "
      f"{bzm - F(1,2) == F(139,458752)}")
assert bzm - F(1, 2) == F(139, 458752)
print(f"    so B_zm > 1/2: the strengthened statement is FALSE here")

print("\n    the stated optimal state at O={1,3}:")
u = [F(0), F(191, 196), F(0), F(5, 196)]
print(f"      u = {[str(t) for t in u]};  sum u = {sum(u)}  (must be k-q = 1): "
      f"{sum(u) == 1}")
r = [d[i]*(u[i] - f[i]) for i in range(k)]
print(f"      r = {[str(t) for t in r]}")
paper_r = [F(-195, 1024), F(316875, 458752), F(-29341, 65536), F(-377, 7168)]
print(f"      matches the paper: {r == paper_r}")
assert r == paper_r
print(f"      sum r = {sum(r)}  (the zero-sum identity): {sum(r) == 0}")
assert sum(r) == 0
mx = max(sum(r[i] for i in P) for P in PREF)
print(f"      max prefix value = {mx} = B_zm: {mx == bzm}")
assert mx == bzm

print("\n    the rank-only optimum on the same instance:")
B = min(v for v in (solve(O, False) for O in combinations(range(k), 2))
        if v is not None)
print(f"      B = {B} = {float(B):.12f}   paper says 1606605/67108864 = "
      f"{float(F(1606605,67108864)):.12f}")
print(f"      match: {B == F(1606605,67108864)};  B < 1/2: {B < F(1,2)}")
assert B == F(1606605, 67108864) < F(1, 2)
print("      => the rank-only HALF statement is NOT refuted by this instance")

print("=" * 74)
print("[B] the 4-2sqrt(3) envelope of the three-support cell")
a = sp.Symbol('a', positive=True)
g = a*(2-a)/(1+a)
gp = sp.simplify(sp.diff(g, a))
print(f"    g(a) = a(2-a)/(1+a),  g'(a) = {sp.factor(sp.numer(sp.together(gp)))}"
      f" / (1+a)^2")
num = sp.numer(sp.together(gp))
print(f"    paper says the numerator is 2-2a-a^2: "
      f"{sp.simplify(num - (2-2*a-a**2)) == 0}")
assert sp.simplify(num - (2 - 2*a - a**2)) == 0
roots = sp.solve(sp.Eq(2 - 2*a - a**2, 0), a)
print(f"    stationary points: {roots}")
astar = sp.sqrt(3) - 1
print(f"    in (0,1): a* = sqrt(3)-1 = {sp.N(astar,10)}: "
      f"{0 < sp.N(astar) < 1}")
val = sp.simplify(g.subs(a, astar))
print(f"    g(a*) = {sp.simplify(sp.radsimp(val))} = {sp.N(val,15)}")
print(f"    equals 4-2sqrt(3) = {sp.N(4-2*sp.sqrt(3),15)}: "
      f"{sp.simplify(val - (4-2*sp.sqrt(3))) == 0}")
assert sp.simplify(val - (4 - 2*sp.sqrt(3))) == 0
print(f"    4-2sqrt(3) > 1/2: {4-2*float(sp.sqrt(3)) > 0.5}  "
      f"(so the cell supremum exceeds one half)")

print("\nP2_V12_ZEROSUM_OK")
