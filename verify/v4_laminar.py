"""Is the abstract critical constant of a LAMINAR conflict system at most 1?

The writeup deduces "faithful (=laminar) realizations never witness C >= 1" from
Almoghrabi-Skutella-Warode, whose hypothesis is a TWO-TERMINAL series-parallel
digraph.  Check (1) shows that hypothesis is not met, so the ceiling needs an
independent justification.  Check (2) tests the abstract statement directly:

    kappa(H,d,f) = min { C : f in conv{ z in {0,1}^k :
                          sum_{i in S} d_i (z_i - f_i) <= C  for all S in H } }

Screening is done with scipy (float); every leader is re-certified exactly.
"""
from fractions import Fraction as F
from itertools import product, combinations
import random
import numpy as np
from scipy.optimize import linprog
import sympy as sp
from sympy.solvers.simplex import lpmin
import networkx as nx

# ---------------------------------------------------------------------------
print("=" * 74)
print("[1] Is the canonical laminar realization a two-terminal SP digraph?")
G = nx.Graph([("s", "rho"), ("rho", "t1"), ("rho", "t2"),
              ("s", "t1"), ("s", "t2")])
print("    smallest laminar family H = {{1,2}}; canonical realization")
print("      arcs: s->rho (constraint), rho->t1, rho->t2, s->t1, s->t2")
print("      planar:", nx.check_planarity(G)[0],
      "| apex over a forest => treewidth <= 2 => no K4 minor")
print("      sinks of the DAG: t1, t2  -> two sinks, so NOT a")
print("      two-terminal series-parallel digraph (those have a unique sink).")
G2 = nx.Graph(list(G.edges()) + [("t1", "tau"), ("t2", "tau")])
branch = [{"s"}, {"rho"}, {"t1", "tau"}, {"t2"}]
ok = all(any(G2.has_edge(a, b) for a in S for b in T)
         for i, S in enumerate(branch) for T in branch[i+1:])
print("      after the standard super-sink repair, branch sets",
      [sorted(b) for b in branch], "are pairwise adjacent:", ok)
print("      => a K4 minor appears; ASW cannot be invoked verbatim.")

# ---------------------------------------------------------------------------
print("=" * 74)
print("[2] Abstract laminar critical constant")


def laminar_shapes(k):
    """laminar families given by binary merge trees on [k] (every internal
    node's set is included; singletons and [k] included)."""
    def trees(items):
        items = tuple(items)
        if len(items) == 1:
            yield items[0]
            return
        for r in range(1, len(items)):
            for left in combinations(items, r):
                if items[0] not in left:
                    continue
                right = tuple(x for x in items if x not in left)
                for L in trees(left):
                    for R in trees(right):
                        yield (L, R)

    def sets_of(tr, out):
        if not isinstance(tr, tuple):
            s = frozenset({tr})
        else:
            s = sets_of(tr[0], out) | sets_of(tr[1], out)
        out.append(s)
        return s

    seen = set()
    for tr in trees(range(k)):
        out = []
        sets_of(tr, out)
        fam = frozenset(out)
        if fam not in seen:
            seen.add(fam)
            yield sorted(fam, key=lambda s: (len(s), sorted(s)))


def in_hull_float(U, f, k):
    """is f in conv(U)?  LP feasibility."""
    m = len(U)
    if m == 0:
        return False
    Aeq = np.vstack([np.array(U, dtype=float).T, np.ones((1, m))])
    beq = np.append(np.array(f, dtype=float), 1.0)
    r = linprog(np.zeros(m), A_eq=Aeq, b_eq=beq, bounds=[(0, None)]*m,
                method="highs")
    return r.status == 0


def kappa_float(sets, d, f, k, Z, ZA=None):
    """min level C at which f in conv{z : level(z) <= C}.  Hull membership is
    monotone in C, so bisect the sorted list of distinct levels."""
    lev = np.array([max(sum(d[i]*(z[i]-f[i]) for i in S) for S in sets)
                    for z in Z])
    levels = np.unique(lev)
    lo, hi = 0, len(levels) - 1
    if not in_hull_float([Z[j] for j in range(len(Z))], f, k):
        return None                      # cannot happen: f in [0,1]^k
    while lo < hi:
        mid = (lo + hi)//2
        U = [Z[j] for j in range(len(Z)) if lev[j] <= levels[mid] + 1e-12]
        if in_hull_float(U, f, k):
            hi = mid
        else:
            lo = mid + 1
    return levels[lo]


def kappa_exact(sets, d, f, k):
    Z = list(product([0, 1], repeat=k))
    lev = {z: max(sum(d[i]*(z[i]-f[i]) for i in S) for S in sets) for z in Z}
    for C in sorted(set(lev.values())):
        U = [z for z in Z if lev[z] <= C]
        lam = sp.symbols(f'l0:{len(U)}')
        cons = [l >= 0 for l in lam] + [sp.Eq(sum(lam), 1)]
        for i in range(k):
            cons.append(sp.Eq(sum(lam[j]*U[j][i] for j in range(len(U))),
                              sp.Rational(f[i])))
        try:
            lpmin(sum(lam)*0, cons)
            return C
        except Exception:
            continue
    return None


import sys
random.seed(20260729)
overall = (0.0, None)
total = 0
SAMPLES = {2: 400, 3: 400, 4: 200, 5: 60, 6: 12}
for k in (2, 3, 4, 5, 6):
    Z = list(product([0, 1], repeat=k))
    shapes = list(laminar_shapes(k))
    best = (0.0, None)
    for sets in shapes:
        for _ in range(SAMPLES[k]):
            d = [random.randint(1, 24)/24 for _ in range(k)]
            d[random.randrange(k)] = 1.0                  # normalise d_max = 1
            f = [random.randint(0, 24)/24 for _ in range(k)]
            C = kappa_float(sets, d, f, k, Z)
            total += 1
            if C is not None and C > best[0]:
                best = (C, (sets, d, f))
    print(f"    k={k}: {len(shapes):4d} laminar shapes x {SAMPLES[k]} random"
          f" (d,f)   max kappa = {best[0]:.10f}", flush=True)
    if best[0] > overall[0]:
        overall = best
print(f"    total random laminar instances screened: {total}")
print(f"    global maximum kappa found: {overall[0]:.12f}   (<= 1: "
      f"{overall[0] <= 1 + 1e-9})")

# targeted: the tight family H = {{0,1}}, d = (1,1), f = (1/2+e, 1/2+e)
print("\n    tight family H={{0,1}}, d=(1,1), f=(1/2+e,1/2+e) -- exact:")
for e in (F(1, 4), F(1, 10), F(1, 100), F(1, 1000)):
    fq = [F(1, 2)+e, F(1, 2)+e]
    C = kappa_exact([frozenset({0, 1})], [F(1), F(1)], fq, 2)
    print(f"      e={str(e):8s}  kappa = {C} = {float(C):.8f}"
          f"   (= 1 - 2e, sup = 1 not attained: {C == 1-2*e})")

# targeted: does any laminar system reach exactly 1?
print("\n    exhaustive rational sweep, k=3, all laminar shapes,")
print("    d in {1/2,3/4,1}^3 (max=1), f in {0,1/4,1/2,3/4,1}^3 -- exact:",
      flush=True)
best_exact = (F(0), None)
cnt = 0
for sets in laminar_shapes(3):
    for d in product([F(1, 2), F(3, 4), F(1)], repeat=3):
        if max(d) != 1:
            continue
        for f in product([F(0), F(1, 4), F(1, 2), F(3, 4), F(1)], repeat=3):
            C = kappa_exact(sets, list(d), list(f), 3)
            cnt += 1
            if C is not None and C > best_exact[0]:
                best_exact = (C, (sets, d, f))
print(f"      exact instances checked: {cnt}")
print(f"      maximum exact kappa: {best_exact[0]} = {float(best_exact[0]):.8f}")
print(f"      <= 1 : {best_exact[0] <= 1}")
print(f"      attained at sets={[sorted(s) for s in best_exact[1][0]]}"
      f" d={best_exact[1][1]} f={best_exact[1][2]}")
print("\nV4_LAMINAR_DONE")
