"""Recompute |V^up| and |B| for the k=17 record directly from the instance, and
locate the zero-sum q=3 numbers quoted in the paper."""
import json
import glob
import os
from fractions import Fraction as F

C = '../certificates/'
S = 10**9
rec = json.load(open(C + 'commonpoint_joint_k17_refined_record_exact.json'))
k = rec['k']
IV = [tuple(t) for t in rec['intervals']]
Dn = rec['demand_numerators_over_1e9']
Fn = rec['share_numerators_over_1e9']
NP = max(b for a, b in IV) + 1
rows = [[i for i in range(k) if IV[i][0] <= e <= IV[i][1]] for e in range(NP)]
LATE = [Dn[i]*(S - Fn[i]) for i in range(k)]
EARLY = [-Dn[i]*Fn[i] for i in range(k)]
UNIT = S*S                       # overload 1 in the scaled integers


def rho(z):
    """max overload over all arcs, in units of 1/S^2"""
    val = [LATE[i] if z >> i & 1 else EARLY[i] for i in range(k)]
    m = max(sum(val[i] for i in r) for r in rows)
    p = max(LATE[i] if z >> i & 1 else -EARLY[i] for i in range(k))
    return max(m, p)


print("=" * 74)
print("[A] |B| and |V^up| for the k=17 record, recomputed from the instance")
B = [z for z in range(1 << k) if rho(z) <= UNIT]
print(f"    |B| = #\\{{z : rho(z) <= d_max\\}} = {len(B)}")
cert = json.load(open(C + 'commonpoint_sandwich_phase2_exact.json'))
row17 = [r for r in cert['certified_record_chain'] if r['k'] == 17][0]
print(f"    certificate 'dgg_bases'                  = {row17['dgg_bases']}")
lb = json.load(open(C + 'commonpoint_joint_k17_refined_lower_box_exact.json'))
print(f"    lower-box cert 'upper_good_routings'     = {lb['upper_good_routings']}")
assert len(B) == row17['dgg_bases'] == lb['upper_good_routings'] == 4799
print("    all three agree: B is exactly the upper-good set of Proposition 4.4")

Bs = set(B)
V = set(B)
for z in B:
    for i in range(k):
        if not (z >> i & 1):
            V.add(z | (1 << i))
print(f"    |V^up| = |B u {{z+e_i : z in B, z_i=0}}| = {len(V)}")
print(f"    paper's table says                        19566")
print(f"    match: {len(V) == 19566}")

print()
print("=" * 74)
print("[B] where do the zero-sum q=3 numbers live?")
NEEDLES = {
    '5470154188105031': 'modified d_4 numerator',
    '6017794294600000': 'modified d_4 denominator',
    '178129713169798083': 'scalar proportional wall numerator',
    '308033000000000000': 'scalar proportional wall denominator',
    '431852198179': 'B_zm numerator',
    '24113213169798083': 'wall minus 1/2 numerator',
}
found = {n: [] for n in NEEDLES}
for path in glob.glob(C + '*.json'):
    try:
        blob = open(path, encoding='utf-8', errors='ignore').read()
    except Exception:
        continue
    for n in NEEDLES:
        if n in blob:
            found[n].append(os.path.basename(path))
for n, what in NEEDLES.items():
    fs = found[n]
    print(f"    {what:<38} {'FOUND in ' + fs[0] if fs else 'NOT FOUND in any certificate'}"
          + (f" (+{len(fs)-1} more)" if len(fs) > 1 else ""))

print()
print("=" * 74)
print("[C] arithmetic of the quoted zero-sum values")
a = F(178129713169798083, 308033000000000000)
print(f"    wall = {a} = {float(a):.15f}")
print(f"    wall - 1/2 = {a - F(1,2)}  (paper says 24113213169798083/308033000000000000)")
print(f"    exact match: {a - F(1,2) == F(24113213169798083, 308033000000000000)}")
b = F(431852198179, 2000000000000)
print(f"    B_zm = {b} = {float(b):.15f} < 1/2: {b < F(1,2)}")
d4 = F(5470154188105031, 6017794294600000)
print(f"    modified d_4 = {float(d4):.15f}  in (0,1]: {0 < d4 <= 1}")
print("\nP2_V11_DONE")
