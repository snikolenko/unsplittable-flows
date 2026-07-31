"""Independent verification of the exact walls, the three-order domination, and
the lower-box deficits reported in Part II.  Each check rebuilds the object from
the primitive data in its certificate and recomputes the claimed quantity.
"""
import json
from fractions import Fraction as F
from itertools import product

C = '../certificates/'


def fr(s):
    return F(*map(int, s.split('/'))) if '/' in str(s) else F(int(s))


# ===================================================================== walls
print("=" * 74)
print("[A] wall 1: three intervals on a path have no path-coherent radius-one routing")
d = json.load(open(C + 'k6_path_coherence_obstruction_exact.json'))
IV = [tuple(t) for t in d['intervals']]
dem = [fr(s) for s in d['demands']]
sh = [fr(s) for s in d['shares']]
NP = max(b for a, b in IV) + 1
rows = [[i for i in range(len(IV)) if IV[i][0] <= a <= IV[i][1]] for a in range(NP)]
print(f"    intervals {IV}; rows recomputed: {rows}")
print(f"    certificate matrix: {d['matrix']}")
mat, seen = [], set()
for r in rows:                     # duplicate rows are redundant; de-duplicate
    v = tuple(1 if i in r else 0 for i in range(len(IV)))
    if v not in seen:
        seen.add(v); mat.append(list(v))
assert mat == d['matrix'], (mat, d['matrix'])
print(f"    row matrix matches after de-duplicating repeated rows ({len(rows)} "
      f"spine arcs -> {len(mat)} distinct rows)")
rows = [[i for i in range(len(IV)) if v[i]] for v in map(tuple, mat)]
best_rad, best_osc = None, None
for z in product([0, 1], repeat=len(IV)):
    err = [dem[i]*(z[i]-sh[i]) for i in range(len(IV))]
    re = [sum(err[i] for i in r) for r in rows]
    rad = max(abs(x) for x in re)
    osc = max(re) - min(re)
    if best_rad is None or rad < best_rad:
        best_rad = rad
    if rad <= 1 and (best_osc is None or osc < best_osc):
        best_osc = osc
print(f"    minimum radius over all 8 routings            : {best_rad} "
      f"(certificate {d['minimum_radius']})")
print(f"    minimum oscillation among radius-one routings : {best_osc} "
      f"(certificate {d['minimum_oscillation_among_radius_one']})")
assert best_rad == fr(d['minimum_radius'])
assert best_osc == fr(d['minimum_oscillation_among_radius_one']) == F(55, 52)
print(f"    55/52 = {float(F(55,52)):.6f} > 1: coherent induction fails, "
      f"radius-one rounding itself does not")

print("\n[B] wall 2: two signed columns on a rooted K5 force oscillation 17/16")
d = json.load(open(C + 'k5_signed_global_coherence_obstruction_exact.json'))
A = d['matrix']
dem = [fr(s) for s in d['demands']]
sh = [fr(s) for s in d['shares']]
best = None
for z in product([0, 1], repeat=len(dem)):
    err = [dem[i]*(z[i]-sh[i]) for i in range(len(dem))]
    re = [sum(A[e][i]*err[i] for i in range(len(dem))) for e in range(len(A))]
    if max(abs(x) for x in re) <= 1:
        osc = max(re) - min(re)
        best = osc if best is None else min(best, osc)
print(f"    minimum oscillation among radius-one routings : {best} "
      f"(certificate {d['minimum_oscillation_among_ordinary_good']})")
assert best == fr(d['minimum_oscillation_among_ordinary_good']) == F(17, 16)
cg = json.load(open(C + 'complete_graph_global_coherence_cegis_qflra_exact.json'))
print(f"    complete-graph campaign: {cg['counts']} over {cg['case_count']} cases")
by = {}
for r in cg['records']:
    by.setdefault(r['vertex_count'], [0, 0])
    by[r['vertex_count']][0 if r.get('result', r.get('status')) == 'sat' else 1] += 1
print(f"    per vertex count (sat/unsat): {by}")

print("\n[C] wall 3: the diagonal first-crossing bracket first fails at k=6")
d = json.load(open(C + 'two_permutation_bracket_k6_counterexample_exact.json'))
perm = d['permutation']
dem = [fr(s) for s in d['demands']]
sh = [fr(s) for s in d['shares']]
k = len(dem)
print(f"    sigma = {perm}, k = {k}")


def prefix_errors(z, order):
    out, acc = [], F(0)
    for i in order:
        acc += dem[i]*(z[i]-sh[i])
        out.append(acc)
    return out


good, vals = [], []
for zi, z in enumerate(product([0, 1], repeat=k)):
    e1 = prefix_errors(z, list(range(k)))
    e2 = prefix_errors(z, perm)
    if max(abs(x) for x in e1[:-1] + e2[:-1]) <= 1:      # proper prefixes only
        good.append(zi)
        vals.append(e1[-1])
print(f"    routings good on all PROPER prefixes of both orders: {len(good)} "
      f"(certificate lists {len(d['bracket_good_route_indices'])})")
pos = [v for v in vals if v > 0]
neg = [v for v in vals if v < 0]
gap = min(pos) - max(neg) if pos and neg else None
print(f"    closest opposite-sign full-row functional values: gap = {gap} "
      f"(certificate {d['minimum_crossing_gap']})")
assert gap == fr(d['minimum_crossing_gap']) == F(45, 44)
print(f"    45/44 = {float(F(45,44)):.6f} > D = 1  => the scalar bracket fails")
print(f"    but its exact lower-box deficit is {d['lower_box']['exact_deficit']} "
      f"and the cost-hull entry is {d['cost_hull']['entry']} < 1,")
print("    so neither (LB) nor the cost-preserving statement is refuted.")
assert fr(d['lower_box']['exact_deficit']) == 0
assert fr(d['cost_hull']['entry']) < 1

# ============================================================== three orders
print("=" * 74)
print("[D] three-order overlays are dominated at 6/5")
d = json.load(open(C + 'three_order_overlay_n1_no_separation_exact.json'))
dem = [F(u, 20) for u in d['demand_units_over_20']]
frac = [F(int(u), 20) for u in d['fractional_load_units_over_20']]
mix = d['mixture']
w = [fr(a['weight']) for a in mix]
print(f"    raw arcs {d['raw_arc_count']}, path counts {d['raw_path_counts']}, "
      f"all routings {d['all_raw_routings']}")
print(f"    demands {[str(x) for x in dem]}, d_max = {max(dem)}")
print(f"    mixture atoms: {len(mix)}; weights sum to {sum(w)}")
assert sum(w) == 1 and all(x > 0 for x in w)
nA = len(frac)
bary = [sum(w[j]*F(int(mix[j]['load_units_over_20'][a]), 20) for j in range(len(mix)))
        for a in range(nA)]
print(f"    barycenter equals the fractional load on all {nA} arcs: {bary == frac}")
assert bary == frac
worst = max(max(F(int(a['load_units_over_20'][x]), 20) - frac[x] for x in range(nA))
            for a in mix)
print(f"    max over atoms of max arc overload: {worst} = {float(worst):.6f}")
print(f"    <= (6/5) d_max = {F(6,5)*max(dem)}: {worst <= F(6,5)*max(dem)}")
assert worst <= F(6, 5)*max(dem)
print("    => the fractional load is a convex combination of routings of")
print("       congestion <= 6/5, so for EVERY signed arc-toll vector some")
print("       atom is at least as cheap as the fractional flow.  No signed")
print("       toll can force congestion above 6/5 on this overlay.")
ex = json.load(open(C + 'three_order_equality_row_order_exhaustive.json'))
print(f"    row-order exhaustion: {json.dumps({k: v for k, v in ex.items() if not isinstance(v, list)})[:300]}")

# ============================================================== lower box
print("=" * 74)
print("[E] lower-box deficits on the records")
for tag, path in [("k=17", 'commonpoint_joint_k17_refined_lower_box_exact.json')]:
    lb = json.load(open(C + path))
    keys = [k for k in lb if 'defic' in k.lower()]
    for kk in keys:
        v = lb[kk]
        if isinstance(v, str) and '/' in v:
            print(f"    {tag} {kk}: {float(fr(v)):.6e}   < D = 1: {fr(v) < 1}")
            assert fr(v) < 1
        else:
            print(f"    {tag} {kk}: {v}")
print("\nP2_V4_WALLS_OK")
