"""Independent verification of the k=17 common-point record of Part II.

Reads ONLY (d, f, intervals) from the certificate; rebuilds the raw digraph,
rediscovers all simple source-terminal paths by DFS, enumerates all 2^17
routings in exact integer arithmetic, and recomputes the critical constant,
the minimizers, the rank statement, and the convex-hull certificate.
"""
import json
from fractions import Fraction as F
from itertools import combinations

CERT = '../certificates/commonpoint_joint_k17_refined_record_exact.json'
cert = json.load(open(CERT))
S = 10**9

k = cert['k']
IV = [tuple(t) for t in cert['intervals']]
d = [F(n, S) for n in cert['demand_numerators_over_1e9']]
f = [F(n, S) for n in cert['share_numerators_over_1e9']]
m = cert['rank']
claimed = F(*map(int, cert['exact_level'].split('/')))

print("=" * 74)
print("[1] instance data as read from the certificate")
print(f"    k = {k}, claimed rank m = {m}")
print(f"    max d_i = {max(d)}   (normalised d_max = 1: {max(d)==1})")
print(f"    sum f_i = {sum(f)} = {float(sum(f)):.12f}")
print(f"    m - 1 + eta with eta = 1e-9: {m-1+F(1,S)}  "
      f"match: {sum(f) == m-1+F(1,S)}")
assert max(d) == 1 and sum(f) == m - 1 + F(1, S)
assert all(0 < x <= 1 for x in d) and all(0 <= x <= 1 for x in f)

# common point?
lo = max(a for a, b in IV)
hi = min(b for a, b in IV)
print(f"    every interval contains arcs [{lo},{hi}]  -> common point: {lo<=hi}")
assert lo <= hi
NPOS = max(b for a, b in IV) + 1
print(f"    spine arcs 0..{NPOS-1}  ({NPOS} shared rows)")

# ---------------------------------------------------------------- raw digraph
print("=" * 74)
print("[2] raw digraph rebuilt from the intervals, paths found by DFS")
spine = [(f"v{j}", f"v{j+1}") for j in range(NPOS)]
arcs = list(spine)
early_arc, late_arc = {}, {}
for i, (a, b) in enumerate(IV):
    early_arc[i] = (f"v{a}", f"t{i}")
    late_arc[i] = (f"v{b+1}", f"t{i}")
    arcs += [early_arc[i], late_arc[i]]
print(f"    arcs built: {len(arcs)}   (certificate says {cert['raw_dag']['arc_count']})")
assert len(arcs) == cert['raw_dag']['arc_count']


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
    return sorted(out, key=lambda p: (len(p), p))


paths = {i: simple_paths(arcs, "v0", f"t{i}") for i in range(k)}
counts = [len(paths[i]) for i in range(k)]
print(f"    simple path counts per terminal: {set(counts)}  (all == 2: {set(counts)=={2}})")
assert counts == [2]*k
for i in range(k):
    assert paths[i][0][-1] == early_arc[i] and paths[i][1][-1] == late_arc[i]
print("    the two DFS paths per terminal are exactly the early/late pair")

# derive the interval of each terminal from the discovered paths
SPI = {a: j for j, a in enumerate(spine)}
derived = []
for i in range(k):
    e = {SPI[a] for a in paths[i][0] if a in SPI}
    l = {SPI[a] for a in paths[i][1] if a in SPI}
    assert e < l
    dl = sorted(l - e)
    derived.append((dl[0], dl[-1]))
    assert dl == list(range(dl[0], dl[-1]+1)), "not an interval!"
print(f"    intervals re-derived from the raw paths match: {derived == IV}")
assert derived == IV

R = [sorted(i for i in range(k) if a <= IV[i][1] and IV[i][0] <= a)
     for a in range(NPOS)]
print(f"    row supports: {len(R)} rows; two chains "
      f"(nested left of {lo}, nested right of it): ", end="")
left = all(set(R[a]) <= set(R[a+1]) for a in range(lo))
right = all(set(R[a]) >= set(R[a+1]) for a in range(hi, NPOS-1))
print(f"{left and right}")
assert left and right

# ------------------------------------------------------- exhaustive enumeration
print("=" * 74)
print("[3] exhaustive enumeration of all 2^17 routings (exact integers)")
Dn = cert['demand_numerators_over_1e9']
Fn = cert['share_numerators_over_1e9']
# overload*S^2 on a spine row = sum_{i in R} Dn_i*(S*z_i - Fn_i)   (exact ints)
LATE = [Dn[i]*(S - Fn[i]) for i in range(k)]      # z_i = 1
EARLY = [-Dn[i]*Fn[i] for i in range(k)]          # z_i = 0
PRIV = [max(Dn[i]*(S - Fn[i]), Dn[i]*Fn[i]) for i in range(k)]
rowidx = [R[a] for a in range(NPOS)]

best = None
minimizers = []
ngood = 0
sub_max_rank = -1
nsub = 0
lvl = claimed.numerator * S**2 // claimed.denominator   # exact: denom | S^2
assert claimed.denominator * lvl == claimed.numerator * S**2

for zmask in range(1 << k):
    z = [(zmask >> i) & 1 for i in range(k)]
    rank = sum(z)
    val = [LATE[i] if z[i] else EARLY[i] for i in range(k)]
    mx = max(sum(val[i] for i in idx) for idx in rowidx)
    pv = max(PRIV[i] if False else (LATE[i] if z[i] else -EARLY[i])
             for i in range(k))
    mx = max(mx, pv)
    if mx < lvl:
        nsub += 1
        sub_max_rank = max(sub_max_rank, rank)
    if rank < m:                       # cost-good  <=>  at least m late paths
        continue
    ngood += 1
    if best is None or mx < best:
        best, minimizers = mx, [zmask]
    elif mx == best:
        minimizers.append(zmask)

got = F(best, S**2)
print(f"    cost-good routings (rank >= {m}) : {ngood}")
print(f"    minimum all-arc overload         : {got}")
print(f"    claimed exact level              : {claimed}")
print(f"    match                            : {got == claimed}")
print(f"    decimal                          : {float(got):.18f}")
print(f"    number of minimizers             : {len(minimizers)}  "
      f"(certificate says {len(cert['minimizers'])})")
print(f"    strict-sublevel routings         : {nsub}  "
      f"(certificate says {cert['strict_sublevel_count']})")
print(f"    max rank among strict-sublevel   : {sub_max_rank}  "
      f"(certificate says {cert['strict_sublevel_maximum_rank']})")
assert got == claimed
assert len(minimizers) == len(cert['minimizers']) == 15
assert nsub == cert['strict_sublevel_count']
assert sub_max_rank == cert['strict_sublevel_maximum_rank'] == m - 1

certmin = {int(s[::-1], 2) for s in cert['minimizers']}
print(f"    minimizer sets agree             : {set(minimizers) == certmin}")
assert set(minimizers) == certmin

print(f"    > 9/8  : {got > F(9,8)}    > 6/5 : {got > F(6,5)}   "
      f"< 4/3 : {got < F(4,3)}")
print(f"    exceeds Part I record (299-41 sqrt41)/32 = 1.13974707... : "
      f"{float(got) > 1.139747070789}")

# ------------------------------------------------------------ hull certificate
print("=" * 74)
print("[4] convex-hull certificate (barycenter = fractional shares)")
mix = cert['hull_mixture']
wts = [F(*map(int, a['weight'].split('/'))) if '/' in a['weight']
       else F(int(a['weight'])) for a in mix]
print(f"    atoms: {len(mix)}   weights sum to 1: {sum(wts) == 1}")
assert sum(wts) == 1 and all(w > 0 for w in wts)
bary = [sum(w * int(a['route'][i]) for w, a in zip(wts, mix)) for i in range(k)]
print(f"    barycenter equals f              : {bary == f}")
assert bary == f
lvls = []
for a in mix:
    z = [int(ch) for ch in a['route']]
    val = [LATE[i] if z[i] else EARLY[i] for i in range(k)]
    mx = max(max(sum(val[i] for i in idx) for idx in rowidx),
             max(LATE[i] if z[i] else -EARLY[i] for i in range(k)))
    lvls.append(F(mx, S**2))
    stated = F(*map(int, a['overload'].split('/')))
    assert lvls[-1] == stated, (lvls[-1], stated)
print(f"    every atom's overload recomputed and matches the certificate")
print(f"    max atom overload                : {max(lvls)} = {float(max(lvls)):.12f}")
print(f"    all atoms at or below the level  : {max(lvls) <= got}")
assert max(lvls) <= got
print("    => f lies in conv{routings of overload <= C*}, so by the hull")
print("       characterisation C* is exactly the all-cost critical constant.")

print("=" * 74)
print("[5] firewalls")
dist = sorted(set(d))
div = [(a, b) for a in dist for b in dist if a < b and (b/a).denominator == 1]
print(f"    distinct demands: {len(dist)}; divisibility relations: {div or 'none'}")
assert not div
print(f"    cost: c^T x = {k} - sum f = {k - sum(f)}; cost-good iff #early <= "
      f"{k-m} iff rank >= {m}")
print("\nP2_V1_K17_OK")
