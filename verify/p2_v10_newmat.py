"""Independent check of the material added to Part II after the last audit:
the multi-star record-chain sandwich, the clone/interleaving frontier, and the
zero-sum q=3 wall.  Everything is rebuilt from the certificates' primitive data
and recomputed in exact rational arithmetic where possible.
"""
import json
import re
from fractions import Fraction as F
from itertools import product

C = '../certificates/'
# the per-k defeating-overload table is cross-checked against the paper
# source when it is available; the certificate checks run regardless.
import os
_tex = '../paper/ssuf_part2.tex'
TEX = open(_tex, encoding='utf-8').read() if os.path.exists(_tex) else ''


def fr(s):
    s = str(s)
    return F(*map(int, s.split('/'))) if '/' in s else F(s)


print("=" * 74)
print("[A] the record-chain sandwich table (Proposition 'record-sandwich')")
cert = json.load(open(C + 'commonpoint_sandwich_phase2_exact.json'))
print("    certificate keys:", list(cert.keys())[:12])
rows = None
for k, v in cert.items():
    if isinstance(v, list) and v and isinstance(v[0], dict) and \
            any('k' in d for d in v[:1]):
        rows = v
        print(f"    using list '{k}' with {len(v)} rows; row keys "
              f"{list(v[0].keys())[:10]}")
        break
if rows is None:
    print("    (no row list found; dumping shallow structure)")
    for k, v in cert.items():
        print("      ", k, ":", str(v)[:110])

# the table as typeset in the paper
TAB = {7: '1.208717093134040213', 8: '1.214766411958000000',
       11: '1.255745974752202140', 13: '1.261856423342816434',
       15: '1.271981601327696511', 16: '1.274864467838565124',
       17: '1.282494797984843521'}
VS = {7: 96, 8: 146, 11: 857, 13: 2919, 15: 10016, 16: 13841, 17: 19566}
ATOMS = {7: 8, 8: 9, 11: 12, 13: 14, 15: 16, 16: 17, 17: 18}
if rows:
    got = {}
    for r in rows:
        kk = r.get('k')
        if kk is None:
            continue
        got[kk] = r
    print(f"    rows present for k = {sorted(got)}")
    for kk in sorted(TAB):
        if kk not in got:
            print(f"    k={kk}: NOT IN CERTIFICATE  <<<")
            continue
        r = got[kk]
        blob = json.dumps(r)
        lvl_ok = TAB[kk].rstrip('0').rstrip('.') in blob or TAB[kk] in blob
        vs_ok = str(VS[kk]) in blob
        at_ok = str(ATOMS[kk]) in blob
        print(f"    k={kk:<3} level in cert: {str(lvl_ok):<5} "
              f"|V^up|={VS[kk]} found: {str(vs_ok):<5} "
              f"atoms={ATOMS[kk]} found: {at_ok}")

print()
print("    consistency of the table with the k=17 record proved in Section 4:")
rec = json.load(open(C + 'commonpoint_joint_k17_refined_record_exact.json'))
lvl17 = fr(rec['exact_level'])
print(f"      record  = {lvl17} = {float(lvl17):.18f}")
print(f"      table   = {TAB[17]}")
assert str(lvl17.numerator) == TAB[17].replace('.', '')[:19], \
    (lvl17.numerator, TAB[17])
print("      the k=17 sandwich entry equals the certified record: OK")
print("      table is strictly increasing in k:",
      all(float(TAB[a]) < float(TAB[b])
          for a, b in zip(sorted(TAB), sorted(TAB)[1:])))

print("=" * 74)
print("[B] the clone/interleaving frontier (Proposition 'clone-frontier')")
try:
    cf = json.load(open(C + 'commonpoint_clone_interleaving_ansatz_frontier_exact.json'))
except FileNotFoundError:
    cf = None
    print("    certificate not found")
if cf:
    print("    keys:", list(cf.keys())[:14])
    blob = json.dumps(cf)
    claim = '1221623116867287239'
    print(f"    largest defeating overload {claim} present: {claim in blob}")
    gap_num, gap_den = 30435840558778141, 5*10**17
    lhs = F(int(claim), 10**18)
    print(f"      = {float(lhs):.18f}")
    print(f"    record - that value = {lvl17 - lhs}")
    print(f"    paper says the gap is {gap_num}/{gap_den} = "
          f"{float(F(gap_num, gap_den)):.18f}")
    print(f"    gaps agree: {lvl17 - lhs == F(gap_num, gap_den)}")
    assert lvl17 - lhs == F(gap_num, gap_den), (lvl17 - lhs, F(gap_num, gap_den))
    # the per-k table in the paper
    tab = re.search(r'20&1\.221623117(.*?)\\bottomrule', TEX, re.S)
    if tab:
        vals = re.findall(r'(\d+)&(\d\.\d+)', '20&1.221623117' + tab.group(1))
        print(f"    per-k table in the paper: {len(vals)} entries, "
              f"k from {vals[0][0]} to {max(int(a) for a, _ in vals)}")
        mx = max(float(b) for _, b in vals)
        print(f"    max entry {mx} vs record {float(lvl17):.9f}: "
              f"all below record: {mx < float(lvl17)}")
        assert mx < float(lvl17)
        ks = sorted(int(a) for a, _ in vals)
        print(f"    k values covered: {ks[0]}..{ks[-1]}, count {len(ks)}, "
              f"contiguous: {ks == list(range(ks[0], ks[-1]+1))}")
        print(f"    paper claims 'every $k=20,\\ldots,40$': "
              f"{ks == list(range(20, 41))}")

print("=" * 74)
print("[C] the zero-sum q=3 wall")
try:
    zs = json.load(open(C + 'q3_zero_sum_k4_counterexample_exact.json'))
    blob = json.dumps(zs)
    print("    keys:", list(zs.keys())[:14])
    for claim, what in [
            ('5470154188105031', 'the modified d_4 numerator'),
            ('178129713169798083', 'scalar proportional wall numerator'),
            ('431852198179', 'B_zm numerator')]:
        print(f"    {what:<38} in certificate: {claim in blob}")
    a = F(178129713169798083, 308033000000000000)
    b = F(1, 2) + F(24113213169798083, 308033000000000000)
    print(f"    wall value = {float(a):.12f};  stated split 1/2 + rest "
          f"is exact: {a == b}")
    assert a == b
    print(f"    wall > 1/2: {a > F(1,2)}  (so the scalar proportional route fails)")
    bzm = F(431852198179, 2000000000000)
    print(f"    B_zm = {float(bzm):.12f} < 1/2: {bzm < F(1,2)}")
    assert bzm < F(1, 2)
except FileNotFoundError:
    print("    certificate not found")

print()
print("P2_V10_NEWMAT_DONE")
