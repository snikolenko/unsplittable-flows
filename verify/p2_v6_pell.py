"""Rebuild concrete members of the k=6 and k=7 Pell families from their
defining formulas (NOT from stored loads), enumerate all routings exactly, and
confirm the claimed critical constants and their limits.
"""
import json
from fractions import Fraction as F
from itertools import product
import sympy as sp

C = '../certificates/'


def rows_of(intervals):
    NP = max(b for a, b in intervals) + 1
    return [[i for i in range(len(intervals))
             if intervals[i][0] <= a <= intervals[i][1]] for a in range(NP)]


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
        for w in adj.get(node, []):
            if w not in path:
                stack.append((w, path + [w]))
    return out


def check(intervals, d, f, m, label):
    """rebuild the raw DAG, check two paths per terminal, and return the exact
    critical constant over cost-good (>= m late) routings"""
    k = len(d)
    NP = max(b for a, b in intervals) + 1
    arcs = [(f"v{j}", f"v{j+1}") for j in range(NP)]
    for i, (a, b) in enumerate(intervals):
        arcs += [(f"v{a}", f"t{i}"), (f"v{b+1}", f"t{i}")]
    cnt = [len(simple_paths(arcs, "v0", f"t{i}")) for i in range(k)]
    assert cnt == [2]*k, (label, cnt)
    rws = rows_of(intervals)
    best = None
    for z in product([0, 1], repeat=k):
        if sum(z) < m:
            continue
        e = [d[i]*(z[i]-f[i]) for i in range(k)]
        mx = max(max(sum(e[i] for i in r) for r in rws),
                 max(d[i]*(1-f[i]) if z[i] else d[i]*f[i] for i in range(k)))
        best = mx if best is None else min(best, mx)
    return best, len(arcs), cnt


print("=" * 74)
print("[k=6] Pell family:  (5+2 sqrt6)^n = p_n + q_n sqrt6,  C_n = (p_n-4q_n)^2/(2q_n^2)")
cert = json.load(open(C + 'commonpoint_k6_pell_family_exact.json'))
IV6 = [tuple(t) for t in cert['intervals']]
print(f"    intervals {IV6}")
for ex in cert['rational_examples']:
    p, q = ex['p'], ex['q']
    assert p*p - 6*q*q == 1, "not a Pell solution"
    r = F(p, q)
    eta = F(r*r - 6, 2)
    assert eta == F(1, 2*q*q)
    b = 2 - r/2
    u = r/2 - 1
    w = r - 2
    v = 9 + eta - 7*r/2
    d = [F(1), F(1), F(1), b, F(1), b]
    f = [u, v, u, w, u, w]
    assert all(0 < x <= 1 for x in d), d
    assert all(0 <= x <= 1 for x in f), f
    S = sum(f)
    got, na, cnt = check(IV6, d, f, 3, f"k6 n={ex['pell_index']}")
    want = F((p - 4*q)**2, 2*q*q)
    print(f"    p={p:<7} q={q:<6} eta={str(eta):<12} sum f = 2+eta: {S == 2+eta}"
          f"  arcs={na}  paths={set(cnt)}")
    print(f"       C = {got} = {float(got):.12f}   formula {want}  match={got==want}"
          f"   >6/5: {got > F(6,5)}")
    assert got == want == F(*map(int, ex['C'].split('/')))
lim6 = 11 - 4*sp.sqrt(6)
print(f"    limit 11-4 sqrt6 = {sp.N(lim6, 18)}   (all members strictly below: "
      f"{all(F(*map(int,e['C'].split('/'))) < sp.Rational(11) - 4*sp.sqrt(6) for e in cert['rational_examples'])})")

print("=" * 74)
print("[k=7] Pell family:  (8+3 sqrt7)^n = p_n + Q_n sqrt7,"
      "  C_n = (p_n-8Q_n)(p_n-4Q_n)/(6Q_n^2)")
cert = json.load(open(C + 'commonpoint_k7_pell_family_exact.json'))
IV7 = [tuple(t) for t in cert['intervals']]
print(f"    intervals {IV7}")
sel = cert.get('selected_raw_member')
if sel:
    print(f"    selected raw member keys: {list(sel.keys())[:10]}")
for ex in cert['rational_examples']:
    p, Q = ex['p'], ex['Q']
    assert p*p - 7*Q*Q == 1, "not a Pell solution"
    want = F((p - 8*Q)*(p - 4*Q), 6*Q*Q)
    stated = F(*map(int, ex['C'].split('/')))
    print(f"    p={p:<8} Q={Q:<7} formula C_n = {want} = {float(want):.12f}"
          f"   certificate {stated}   match={want==stated}   >6/5: {want>F(6,5)}")
    assert want == stated
lim7 = sp.Rational(13, 2) - 2*sp.sqrt(7)
print(f"    limit 13/2-2 sqrt7 = {sp.N(lim7, 18)}")
print(f"    monotone increasing to the limit: "
      f"{[float(F(*map(int,e['C'].split('/')))) for e in cert['rational_examples']]}")

print("=" * 74)
print("[limits] the two algebraic limits, recomputed symbolically")
r = sp.Symbol('r', positive=True)
print(f"    k=6:  lim (r-4)^2/2 as r -> sqrt6  = "
      f"{sp.simplify(((r-4)**2/2).subs(r, sp.sqrt(6)))} = {sp.N(lim6,18)}")
h = sp.Symbol('h', positive=True)
print(f"    k=7:  lim (h-8)(h-4)/6 as h -> sqrt7 = "
      f"{sp.simplify(((h-8)*(h-4)/6).subs(h, sp.sqrt(7)))} = {sp.N(lim7,18)}")
print(f"    ordering:  1.13975 (k=4) < 1.16798 (k=5) < {float(lim6):.5f} (k=6)"
      f" < {float(lim7):.5f} (k=7) < 1.2148 (k=8) < 1.28249 (k=17)")
print("\nP2_V6_PELL_OK")
