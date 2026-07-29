"""Independent verification of the abstract-side claims:
 * complement-closed identity and the Hadamard bound (brute force, small n);
 * two-scale (chain-firewall) robustness;
 * the stable-set rank-separator LP values (odd holes, antiholes, K_n, Chvatal);
 * the k=7 triple system as an exact failure of the abstract lower box;
 * the lower-box deficit of the public triangle (primal LP and dual separator).
"""
from fractions import Fraction as F
from itertools import product, combinations
import sympy as sp

# ---------------------------------------------------------------------------
print("=" * 74)
print("[A] complement-closed one-sided optimum = disc/2, and Hadamard")


def sylvester(m):
    H = [[1]]
    for _ in range(m):
        H = [r + r for r in H] + [r + [-v for v in r] for r in H]
    return H


for m in (1, 2, 3):
    n = 2**m
    H = sylvester(m)
    sets = []
    for r in range(n):
        sets.append(frozenset(i for i in range(n) if H[r][i] == 1))
        sets.append(frozenset(i for i in range(n) if H[r][i] == -1))
    sets = [s for s in set(sets) if s]
    # brute force: min over chi with T>=0 of max_S chi(S) ; overload = chi(S)/2
    best = None
    for bits in product([-1, 1], repeat=n):
        if sum(bits) < 0:
            continue
        v = max(sum(bits[i] for i in S) for S in sets)
        best = v if best is None else min(best, v)
    # classical discrepancy of the same system
    disc = None
    for bits in product([-1, 1], repeat=n):
        v = max(abs(sum(bits[i] for i in S)) for S in sets)
        disc = v if disc is None else min(disc, v)
    print(f"  n={n:2d}: one-sided opt overload = {F(best,2)}"
          f"   disc/2 = {F(disc,2)}   equal={F(best,2)==F(disc,2)}"
          f"   sqrt(n)/4 = {sp.N(sp.sqrt(n)/4,6)}  bound holds={F(best,2) >= sp.Rational(1,4)*sp.sqrt(n)}")

# two-scale robustness (chain firewall) at n=8
n, m = 8, 3
H = sylvester(m)
sets = [frozenset(i for i in range(n) if H[r][i] == sgn)
        for r in range(n) for sgn in (1, -1)]
sets = [s for s in set(sets) if s]
for pattern in ([1]*4 + [0]*4, [1, 0]*4):
    d = [F(1) if p else 1 - F(1, n) for p in pattern]
    best = None
    for bits in product([0, 1], repeat=n):
        if sum(d[i]*(bits[i] - F(1, 2)) for i in range(n)) < 0:
            continue
        v = max(sum(d[i]*(bits[i] - F(1, 2)) for i in S) for S in sets)
        best = v if best is None else min(best, v)
    lb = sp.sqrt(n)/4 - F(1, 2)
    print(f"  two-scale demands {sorted(set(d))}: opt = {best} = {float(best):.6f}"
          f"   >= sqrt(n)/4-1/2 = {float(lb):.6f} : {float(best) >= float(lb)}")
    assert len(set(d)) == 2 and not any(
        a != b and (max(a, b) / min(a, b)).denominator == 1
        for a in set(d) for b in set(d))

# ---------------------------------------------------------------------------
print("=" * 74)
print("[B] stable-set rank-separator LP:  max t s.t. t<=q_i+q_j, 0<=q<=1, sum q=|V|-alpha")


def rank_lp(edges, nv, alpha):
    q = sp.symbols(f'q0:{nv}', nonnegative=True)
    tt = sp.Symbol('t')
    cons = [tt <= q[i] + q[j] for (i, j) in edges]
    cons += [q[i] <= 1 for i in range(nv)]
    cons += [q[i] >= 0 for i in range(nv)]      # sympy's lp does NOT read assumptions
    cons += [sp.Eq(sum(q), nv - alpha)]
    from sympy.solvers.simplex import lpmax
    val, _ = lpmax(tt, cons + [tt >= 0])
    return sp.nsimplify(val)


def max_stable(nv, edges):
    E = set(map(frozenset, edges))
    best = 0
    for r in range(nv, 0, -1):
        for S in combinations(range(nv), r):
            if all(frozenset(p) not in E for p in combinations(S, 2)):
                return r
    return best


print("  odd holes C_{2m+1}:  claim 1 + 1/(2m+1)")
for mm in (1, 2, 3, 4):
    nv = 2*mm + 1
    edges = [(i, (i+1) % nv) for i in range(nv)]
    a = max_stable(nv, edges)
    v = rank_lp(edges, nv, a)
    claim = 1 + sp.Rational(1, nv)
    print(f"    C_{nv}: alpha={a} (claim {mm})  LP={v}  claim={claim}  match={v==claim}")
    assert a == mm and v == claim

print("  odd antiholes:  claim 2 - 4/n")
for nv in (5, 7, 9, 11):
    cyc = {frozenset((i, (i+1) % nv)) for i in range(nv)}
    edges = [(i, j) for i, j in combinations(range(nv), 2)
             if frozenset((i, j)) not in cyc]
    a = max_stable(nv, edges)
    v = rank_lp(edges, nv, a)
    claim = 2 - sp.Rational(4, nv)
    print(f"    ~C_{nv}: alpha={a} (claim 2)  LP={v}  claim={claim}  match={v==claim}")
    assert a == 2 and v == claim

print("  complete graphs:  claim 2 - 2/n")
for nv in (3, 4, 5, 8):
    edges = list(combinations(range(nv), 2))
    v = rank_lp(edges, nv, 1)
    claim = 2 - sp.Rational(2, nv)
    print(f"    K_{nv}: LP={v}  claim={claim}  match={v==claim}")
    assert v == claim

# Chvatal graph: 12 vertices, 4-regular, alpha=4, claim 4/3
CHV = [(0,1),(0,4),(0,6),(0,9),(1,2),(1,5),(1,7),(2,3),(2,6),(2,8),
       (3,4),(3,7),(3,9),(4,5),(4,8),(5,10),(5,11),(6,10),(6,11),
       (7,8),(7,11),(8,10),(9,10),(9,11)]
a = max_stable(12, CHV)
v = rank_lp(CHV, 12, a)
print(f"    Chvatal: alpha={a} (claim 4)  LP={v}  claim=4/3  match={v==sp.Rational(4,3)}")

# ---------------------------------------------------------------------------
print("=" * 74)
print("[C] abstract lower box fails for the k=7 triple system")
k = 7
rows = [frozenset(S) for S in combinations(range(k), 3)] + [frozenset(range(k))]
good = [z for z in product([0, 1], repeat=k)
        if all(sum(z[i] - F(1, 2) for i in S) <= 1 for S in rows)]
w = [sum(z) for z in good]
print(f"  upper-good binary vectors : {len(good)}   max Hamming weight = {max(w)}")
assert max(w) == 2
print(f"  all-element row: fractional value {F(k,2)}, "
      f"max over conv(U) = {max(w)}, deficit = {F(k,2)-max(w)} > D=1 : "
      f"{F(k,2)-max(w) > 1}")
assert F(k, 2) - max(w) == F(3, 2)
print("  => exact failure of the abstract (LB); the system is non-laminar "
      "(e.g. {0,1,2} and {2,3,4} cross), so Theorem 'splice' forbids a "
      "coverage-faithful realization.")
A, B = frozenset({0, 1, 2}), frozenset({2, 3, 4})
print(f"  crossing witness: {set(A)}, {set(B)}  A&B={set(A&B)} "
      f"A\\B={set(A-B)} B\\A={set(B-A)}")

# ---------------------------------------------------------------------------
print("=" * 74)
print("[D] lower-box deficit of the public triangle (flow instance, exact)")


def simple_paths(arcs, s, tgt):
    adj = {}
    for (u, v) in arcs:
        adj.setdefault(u, []).append(v)
    out, stack = [], [(s, [s])]
    while stack:
        node, path = stack.pop()
        if node == tgt:
            out.append(tuple(zip(path, path[1:])))
            continue
        for ww in adj.get(node, []):
            if ww not in path:
                stack.append((ww, path + [ww]))
    return sorted(out, key=lambda p: (len(p), p))


arcs = [("s","t1"),("s","t2"),("s","u"),("u","t3"),("u","v"),
        ("v","t1"),("v","w"),("w","t2"),("w","t3")]
for epsv in (F(1, 100), F(1, 1000), F(1, 10**6), F(0)):
    d = [F(1), F(2, 3), F(1)]
    f = [F(1, 3), F(1, 2) - epsv, F(1, 3)]
    paths = {i: simple_paths(arcs, "s", tt) for i, tt in enumerate(["t1","t2","t3"])}
    assert [len(paths[i]) for i in range(3)] == [2, 2, 2]
    q = {a: F(0) for a in arcs}
    for i in range(3):
        early, late = paths[i][0], paths[i][1]
        for a in early: q[a] += d[i]*(1-f[i])
        for a in late:  q[a] += d[i]*f[i]
    D = F(1)
    U = []
    for z in product([0, 1], repeat=3):
        y = {a: F(0) for a in arcs}
        for i in range(3):
            for a in paths[i][z[i]]: y[a] += d[i]
        if all(y[a] <= q[a] + D for a in arcs):
            U.append((z, y))
    # primal: min over conv(U) of max_a (q_a - m_a)
    from sympy.solvers.simplex import lpmin, lpmax
    lam = sp.symbols(f'l0:{len(U)}')
    dl = sp.Symbol('delta')
    cons = [sp.Eq(sum(lam), 1)] + [l >= 0 for l in lam]
    for a in arcs:
        cons.append(sum(lam[j]*sp.Rational(U[j][1][a]) for j in range(len(U)))
                    >= sp.Rational(q[a]) - dl)
    val, sol = lpmin(dl, cons)
    val = sp.nsimplify(val)
    # dual separator: max over w >= 0, ||w||_1 = 1 of  w.q - max_{y in U} w.y
    wv = sp.symbols(f'w0:{len(arcs)}')
    th = sp.Symbol('theta')
    cons2 = [sp.Eq(sum(wv), 1)] + [w >= 0 for w in wv]
    for (z, y) in U:
        cons2.append(sum(wv[kk]*sp.Rational(y[a]) for kk, a in enumerate(arcs)) <= th)
    obj = sum(wv[kk]*sp.Rational(q[a]) for kk, a in enumerate(arcs)) - th
    val2, sol2 = lpmax(obj, cons2 + [th >= 0])
    val2 = sp.nsimplify(val2)
    print(f"  eps={epsv}: |U(q)|={len(U)}  primal deficit = {val}"
          f"   dual separator = {val2}   agree={val == val2}"
          f"   (LB holds: {val <= 1})")
    print(f"    optimal separator w: "
          f"{ {arcs[kk]: sol2[wv[kk]] for kk in range(len(arcs)) if sol2[wv[kk]] != 0} }")
    assert val == val2
    if epsv > 0:
        want = sp.Rational(1, 9) - 2*sp.Rational(epsv)/3
        assert val == want, (epsv, val, want)
        print(f"    == 1/9 - 2eps/3 = {want} : True")
    else:
        # hull discontinuity: at eps = 0 the retained-routing set jumps from 4
        # to 7 and the deficit collapses; sup over eps>0 is 1/9, not attained
        assert len(U) == 7 and val == 0
        print("    eps=0: |U| jumps 4 -> 7 and the deficit collapses to 0")
        print("    => sup_{eps>0} delta = 1/9, NOT attained (hull discontinuity)")

print("\nV3_ABSTRACT_OK")
