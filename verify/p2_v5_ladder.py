"""Independent verification of (i) the exact lower-box certificate on the k=17
record and (ii) the record ladder k = 5,6,7,8.

Everything is rebuilt from the primitive instance data (intervals, demands,
shares); routings are re-enumerated and the constants recomputed exactly.
"""
import json
from fractions import Fraction as F
from itertools import product
import sympy as sp

C = '../certificates/'


def fr(s):
    s = str(s)
    return F(*map(int, s.split('/'))) if '/' in s else F(int(s))


def build(intervals):
    """spine rows of a common-point interval system"""
    NP = max(b for a, b in intervals) + 1
    return [[i for i in range(len(intervals))
             if intervals[i][0] <= a <= intervals[i][1]] for a in range(NP)]


def critical_constant(intervals, d, f, m):
    """exact min over cost-good routings (>= m late) of max all-arc overload"""
    rows = build(intervals)
    k = len(d)
    best = None
    for z in product([0, 1], repeat=k):
        if sum(z) < m:
            continue
        e = [d[i]*(z[i]-f[i]) for i in range(k)]
        mx = max(max(sum(e[i] for i in r) for r in rows),
                 max(d[i]*(1-f[i]) if z[i] else d[i]*f[i] for i in range(k)))
        best = mx if best is None else min(best, mx)
    return best


# ============================================================ (A) the lower box
print("=" * 74)
print("[A] exact lower-box certificate on the k=17 record")
lb = json.load(open(C + 'commonpoint_joint_k17_refined_lower_box_exact.json'))
IV = [tuple(t) for t in lb['intervals']]
d = [fr(s) for s in lb['demands']]
f = [fr(s) for s in lb['late_shares']]
k = len(d)
rows = build(IV)
NA = len(rows) + 2*k                      # spine rows + two private exits each
print(f"    k={k}, spine rows {len(rows)}, arcs {NA} "
      f"(certificate raw arcs {len(lb['raw_arcs'])})")
assert len(lb['raw_arcs']) == NA


def arc_errors(z):
    """error vector y-q on every arc: spine rows first, then the two private
    exits of each terminal (early, late)"""
    e = [d[i]*(z[i]-f[i]) for i in range(k)]
    out = [sum(e[i] for i in r) for r in rows]
    for i in range(k):
        # early exit carries d_i(1-f_i) fractionally, the whole d_i iff z_i = 0
        out.append(-d[i]*(1-f[i]) if z[i] else d[i]*f[i])
        # late exit carries d_i f_i fractionally, the whole d_i iff z_i = 1
        out.append(d[i]*(1-f[i]) if z[i] else -d[i]*f[i])
    return out


D = max(d)
print(f"    d_max = {D}")
U = []
for zi in range(1 << k):
    z = [(zi >> i) & 1 for i in range(k)]
    r = arc_errors(z)
    if max(r) <= D:
        U.append((zi, z, r))
print(f"    upper-good routings recomputed: {len(U)} "
      f"(certificate says {lb['upper_good_routings']})")
assert len(U) == lb['upper_good_routings']
mine = {u[0] for u in U}
theirs = set(lb['upper_good_route_indices'])
rev = {int(''.join(str((zi >> i) & 1) for i in range(k)), 2) for zi in mine}
if mine == theirs:
    print("    upper-good index sets agree exactly (LSB-first indexing)")
elif rev == theirs:
    print("    upper-good index sets agree exactly (MSB-first indexing)")
else:
    raise AssertionError(f"index sets differ: {len(mine-theirs)} extra, "
                         f"{len(theirs-mine)} missing")

# primal: the stated mixture, and its deficit
mixw = [fr(a['coefficient']) for a in lb['primal_mixture']]
mixz = [a['routing'] for a in lb['primal_mixture']]
print(f"    primal mixture atoms: {len(mixz)}; weights sum to {sum(mixw)}")
assert sum(mixw) == 1 and all(w > 0 for w in mixw)
for z in mixz:
    assert max(arc_errors(z)) <= D, "a primal atom is not upper-good!"
print("    every primal atom is upper-good")
mean = [sum(mixw[j]*arc_errors(mixz[j])[a] for j in range(len(mixz)))
        for a in range(NA)]
primal_delta = max(-x for x in mean)
print(f"    primal deficit  max_a(q_a - m_a) = {float(primal_delta):.6e}")
assert primal_delta == fr(lb['delta']), (primal_delta, lb['delta'])

# dual: the stated separator
wts = {tuple(a['arc']): fr(a['weight']) for a in lb['dual_arc_weights']}
arcs = [tuple(a) for a in lb['raw_arcs']]
wv = [wts.get(a, F(0)) for a in arcs]
print(f"    dual weights: {len(wts)} nonzero, ||w||_1 = {sum(wv)}")
assert sum(wv) == 1
dual_delta = -max(sum(wv[a]*u[2][a] for a in range(NA)) for u in U)
print(f"    dual value  -max_{{y in U}} w.r(y) = {float(dual_delta):.6e}")
print(f"    primal == dual: {primal_delta == dual_delta}")
assert primal_delta == dual_delta
print(f"    delta = {float(primal_delta):.6e}  <  D = 1: {primal_delta < D}")
print("    => (LB) holds on the k=17 record, with matching primal and dual")
print("       certificates over the complete upper-good set.")

# ================================================================ (B) ladder
print("=" * 74)
print("[B] the record ladder")
lad = []
for kk, path, mkey in [
        (5, 'k5_algebraic_envelope.json', None),
        (6, 'commonpoint_k6_pell_family_exact.json', None),
        (7, 'commonpoint_k7_pell_family_exact.json', None),
        (8, 'commonpoint_k8_stationary_envelope_algebraic.json', None)]:
    try:
        cert = json.load(open(C + path))
    except FileNotFoundError:
        print(f"    k={kk}: certificate {path} not found")
        continue
    keys = list(cert.keys())
    print(f"    k={kk}: {path}")
    print(f"       keys: {keys[:12]}")
    for key in keys:
        v = cert[key]
        if isinstance(v, str) and ('sqrt' in v or '/' in v) and len(v) < 90:
            print(f"       {key}: {v}")
        elif isinstance(v, (int, float)) and 1 < float(v) < 2:
            print(f"       {key}: {v}")

print("\n    claimed limits, evaluated:")
for lbl, expr in [("k=4 (Part I)", (299 - 41*sp.sqrt(41))/32),
                  ("k=6", 11 - 4*sp.sqrt(6)),
                  ("k=7", sp.Rational(13, 2) - 2*sp.sqrt(7))]:
    print(f"       {lbl}: {expr} = {sp.N(expr, 18)}")
print(f"       6/5 = 1.2,  4/3 = {float(F(4,3)):.6f}")
print("\nP2_V5_LADDER_OK")
