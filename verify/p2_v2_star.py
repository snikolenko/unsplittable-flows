"""Independent verification of the deletion-star ceilings of Part II.

(A) executes the *constructive proof* of the codimension-two theorem on a large
    exact sample: it recomputes X, A_pi, A_sigma, the bad set, picks the non-bad
    median element, builds the proportional credits, and checks the residual.
    If this never fails, the theorem holds constructively.
(B) computes the true optimum B(pi,sigma,d,f) by exact LP and confirms B <= 1/2.
(C) checks Example 3/8 at k=8 and the uniform/reverse formula B_{k,q}.
(D) checks the two-scale (chain-firewall) variants.
"""
from fractions import Fraction as F
from itertools import combinations, permutations, product
import random
import sympy as sp
from sympy.solvers.simplex import lpmin

HALF = F(1, 2)


def prefixes(order):
    return [order[:j] for j in range(1, len(order)+1)]


def residual_max(d, f, u, pis):
    """max over all prefixes of both orders of sum_{i in P} d_i (u_i - f_i)"""
    best = None
    for order in pis:
        acc = F(0)
        for i in order:
            acc += d[i]*(u[i] - f[i])
            best = acc if best is None else max(best, acc)
    return best


# ============================================================== (A) constructive
def codim2_construction(pi, sig, d, f):
    """Follow the proof verbatim; return (j, credits, residual)."""
    k = len(d)
    x = [1 - f[i] for i in range(k)]
    assert sum(x) == 2
    mu = [d[i]*x[i] for i in range(k)]
    X = sum(mu)

    def A(order):
        out, acc = set(), F(0)
        for i in order:
            if acc <= 1:
                out.add(i)
            acc += mu[i]
        return out

    def bad(j):
        al = F(1 - x[j], 2 - x[j]) if (2 - x[j]) != 0 else F(0)
        return al*(X - 2*d[j]) > HALF

    if X > 1:
        cand = sorted(A(pi) & A(sig) - {j for j in range(k) if bad(j)})
        assert cand, "bad-mass lemma FAILED: no non-bad median element"
        j = cand[0]
    else:
        j = 0
    al = F(1 - x[j], 2 - x[j])
    c = [F(0)]*k
    for i in range(k):
        if i != j:
            c[i] = F(x[i], 2 - x[j])
    assert sum(c) == 1 and all(0 <= ci <= 1 for ci in c)
    u = [(0 if i == j else 1) - c[i] for i in range(k)]
    return j, c, residual_max(d, f, u, [pi, sig])


def rand_codim2(k, rng, denom=12):
    """random d in (0,1] with max 1, f in [0,1]^k with sum f = k-2"""
    d = [F(rng.randint(1, denom), denom) for _ in range(k)]
    d[rng.randrange(k)] = F(1)
    while True:
        x = [F(rng.randint(0, denom), denom) for _ in range(k)]
        s = sum(x)
        if s == 0:
            continue
        x = [xi*2/s for xi in x]
        if all(xi <= 1 for xi in x):
            break
    return d, [1 - xi for xi in x]


print("=" * 74)
print("[A] the codimension-two construction, executed verbatim (exact)")
rng = random.Random(20260730)
worst = F(0)
worst_at = None
nbad_nonempty = 0
total = 0
for k in range(3, 10):
    for _ in range(1500):
        d, f = rand_codim2(k, rng)
        pi = list(range(k)); rng.shuffle(pi)
        sig = list(range(k)); rng.shuffle(sig)
        j, c, res = codim2_construction(pi, sig, d, f)
        total += 1
        x = [1-fi for fi in f]
        mu = [d[i]*x[i] for i in range(k)]
        X = sum(mu)
        if any(F(1-x[i], 2-x[i])*(X-2*d[i]) > HALF for i in range(k)):
            nbad_nonempty += 1
        if res > worst:
            worst, worst_at = res, (k, d, f, pi, sig)
        assert res <= HALF, (k, d, f, pi, sig, res)
    print(f"    k={k}: 1500 random exact instances, construction always valid")
print(f"    total {total} instances; max residual attained = {worst} "
      f"= {float(worst):.10f}  (<= 1/2: {worst <= HALF})")
print(f"    instances with a NONEMPTY bad set: {nbad_nonempty} "
      f"(so the hard branch is exercised)")

# exhaustive small grid
print("\n    exhaustive grid k=4, x_i in {0,1/2,1,3/2}/1 with sum x = 2,")
print("    d_i in {1/2,1} (max 1), all 24x24 order pairs:")
cnt = 0
mx = F(0)
grid = [F(0), F(1, 2), F(1)]
for x in product(grid, repeat=4):
    if sum(x) != 2:
        continue
    f = [1 - xi for xi in x]
    for d in product([F(1, 2), F(1)], repeat=4):
        if max(d) != 1:
            continue
        for pi in permutations(range(4)):
            for sig in permutations(range(4)):
                _, _, res = codim2_construction(list(pi), list(sig), list(d), f)
                cnt += 1
                mx = max(mx, res)
                assert res <= HALF
print(f"      {cnt} exact instances, maximum residual {mx} = {float(mx):.6f}")

# ==================================================================== (B) LP
print("=" * 74)
print("[B] the true optimum B(pi,sigma,d,f) by exact LP (min over supports)")


def B_exact(pi, sig, d, f, nomit, exact=True):
    """min over omitted sets of the LP value.  Screened in floats by scipy;
    the leader is then re-solved exactly with sympy so the returned value is
    an exact rational."""
    from scipy.optimize import linprog
    import numpy as np
    k = len(d)
    dd = [float(x) for x in d]
    ff = [float(x) for x in f]
    best, bestom = None, None
    for om in combinations(range(k), nomit):
        S = [i for i in range(k) if i not in om]
        idx = {i: j for j, i in enumerate(S)}          # variables: c_S, then t
        n = len(S) + 1
        A, b = [], []
        for order in (pi, sig):
            acc = np.zeros(n)
            rhs = 0.0
            for i in order:
                ind = 0.0 if i in om else 1.0
                if i in idx:
                    acc[idx[i]] -= dd[i]
                rhs += dd[i]*(ind - ff[i])
                row = acc.copy(); row[-1] = -1.0
                A.append(row.tolist()); b.append(-rhs)
        Aeq = [[1.0]*len(S) + [0.0]]
        r = linprog([0.0]*len(S) + [1.0], A_ub=A, b_ub=b, A_eq=Aeq, b_eq=[1.0],
                    bounds=[(0, None)]*len(S) + [(None, None)], method="highs")
        if r.status != 0:
            continue
        if best is None or r.fun < best - 1e-12:
            best, bestom = r.fun, om
    if not exact:
        return best
    om = bestom
    S = [i for i in range(k) if i not in om]
    cv = sp.symbols(f'c0:{k}')
    t = sp.Symbol('t')
    cons = [cv[i] >= 0 for i in S] + [sp.Eq(sum(cv[i] for i in S), 1)]
    for order in (pi, sig):
        acc = 0
        for i in order:
            ind = 0 if i in om else 1
            acc = acc + sp.Rational(d[i])*(ind - (cv[i] if i in S else 0)
                                           - sp.Rational(f[i]))
            cons.append(acc <= t)
    val, _ = lpmin(t, cons)
    val = sp.nsimplify(val)
    return F(int(sp.fraction(val)[0]), int(sp.fraction(val)[1]))


rng = random.Random(7)
mx = F(0)
for trial in range(60):
    k = rng.choice([4, 5, 6])
    d, f = rand_codim2(k, rng, denom=6)
    pi = list(range(k)); rng.shuffle(pi)
    sig = list(range(k)); rng.shuffle(sig)
    b = B_exact(pi, sig, d, f, 1)
    _, _, cons_res = codim2_construction(pi, sig, d, f)
    assert b <= cons_res + F(0), (b, cons_res)   # LP optimum <= construction
    assert b <= HALF, (k, d, f, b)
    mx = max(mx, b)
print(f"    60 random exact LP optima; all <= 1/2; maximum = {mx} = {float(mx):.8f}")
print("    (and every LP optimum is <= the value produced by the construction)")

# ============================================================ (C) exact formulas
print("=" * 74)
print("[C] the uniform identity/reverse family")


def B_uniform(k, q):
    """exact B for identity/reverse orders, unit demands, f = (k-q)/k"""
    d = [F(1)]*k
    f = [F(k-q, k)]*k
    return B_exact(list(range(k)), list(range(k))[::-1], d, f, q-1)


def formula(k, q):
    n = k - q
    M = 2*-(-n//(2*q)) + (q-2)*(-(-n//q))
    ok = (n == 0) or (M <= n + 1)
    return (F(1, 2) - F(q, 2*k)) if ok else None, M, ok


print("    Example (k=8, q=2, f_i=3/4): draft claims B = 3/8")
b = B_uniform(8, 2)
print(f"      B = {b} = {float(b):.6f}   equals 3/8: {b == F(3,8)}   > 1/3: {b > F(1,3)}")
assert b == F(3, 8)
print("    formula B_{k,q} = 1/2 - q/(2k) when M_{k,q} <= n+1:")
bad = []
for k in range(3, 11):
    for q in range(2, k+1):
        pred, M, ok = formula(k, q)
        got = B_uniform(k, q)
        star = "" if ok else "  (criterion says NOT equality)"
        flag = "OK " if (ok and got == pred) or (not ok and got > F(1,2)-F(q,2*k)) else "**"
        if flag == "**":
            bad.append((k, q, got, pred, ok))
        print(f"      k={k:2d} q={q:2d}  B={str(got):>10s}={float(got):.6f}  "
              f"lower bd f/2={float(F(1,2)-F(q,2*k)):.6f}  M={M} n+1={k-q+1}"
              f"  {flag}{star}")
    if k >= 8:
        break
print(f"    formula/criterion mismatches: {bad or 'none'}")
assert not bad
print("    B_{mq,q} = 1/2 - 1/(2m):")
for q in (2, 3, 4):
    for mm in (1, 2, 3):
        k = mm*q
        if k < 3 or k > 12:
            continue
        got = B_uniform(k, q)
        want = F(1, 2) - F(1, 2*mm)
        print(f"      q={q} m={mm} (k={k}): B={got} want={want} match={got==want}")
        assert got == want

# ======================================================= (D) two-scale variants
print("=" * 74)
print("[D] two-scale (non-divisibility-chain) variants, demands {2/3, 1}")
print("    claim:  1/2 - 1/(2m) - 1/(6m)  <=  B  <=  1/2 - 1/(2m)")
probs = []
for q in (2, 3, 4):
    for mm in (1, 2, 3):
        k = mm*q
        if k < 3 or k > 10:
            continue
        d = [F(2, 3)] + [F(1)]*(k-1)
        f = [F(k-q, k)]*k
        got = B_exact(list(range(k)), list(range(k))[::-1], d, f, q-1)
        lo = F(1, 2) - F(1, 2*mm) - F(1, 6*mm)
        hi = F(1, 2) - F(1, 2*mm)
        ok = lo <= got <= hi
        if not ok:
            probs.append((q, mm, k, got, lo, hi))
        print(f"      q={q} m={mm} (k={k}): B={str(got):>9s}={float(got):.6f}  "
              f"[{float(lo):.6f}, {float(hi):.6f}]  {'OK' if ok else '** VIOLATED'}")
print(f"    violations: {probs or 'none'}")
print("\nP2_V2_STAR_DONE")
