"""Independent verification of the two-order upper-bound theorems of Part II:
the Knuth transfer radius B(d), the deque-order class, and the switched-interval
prefix-representation norm.
"""
from fractions import Fraction as F
from itertools import product, permutations, combinations
import random

print("=" * 74)
print("[A] the Knuth transfer radius  B(d) = min_{lambda>=0} (lambda + sum|d_i-lambda|)")


def B(d):
    """exact minimum; the objective is piecewise linear with breakpoints at
    0 and the d_i, so the minimum is attained at one of them"""
    cand = [F(0)] + list(d)
    return min(l + sum(abs(di - l) for di in d) for l in cand)


def lower_median(d):
    s = sorted(d)
    k = len(s)
    return s[(k-1)//2] if k % 2 else s[k//2 - 1]


rng = random.Random(11)
bad = []
for trial in range(4000):
    k = rng.randint(1, 9)
    d = [F(rng.randint(1, 12), 12) for _ in range(k)]
    D = max(d)
    b = B(d)
    # claim 1: minimiser is a lower median
    lm = lower_median(d)
    if lm + sum(abs(di - lm) for di in d) != b:
        bad.append(('median', d))
    # claim 2: D <= B <= ceil(k/2) D
    if not (D <= b <= -(-k//2)*D):
        bad.append(('range', d, b))
    # claim 3: B = D  iff  the k-1 smallest demands are equal
    s = sorted(d)
    eq = (k == 1) or (s[0] == s[k-2])
    if (b == D) != eq:
        bad.append(('locus', d, b, D, eq))
print(f"    4000 random exact demand vectors: violations = {bad[:3] or 'none'}")
assert not bad
print("    lower median is a minimiser                       : confirmed")
print("    D <= B(d) <= ceil(k/2) D                          : confirmed")
print("    B(d) = D  <=>  d_(1) = ... = d_(k-1)              : confirmed")

print("\n    two-scale formula  B/D = min{a + h(1-a), 1 + l(1-a)}  for d in {aD, D}:")
bad2 = []
for a_num in range(1, 12):
    a = F(a_num, 12)
    for lo in range(0, 6):
        for hi in range(0, 6):
            if lo + hi == 0 or hi == 0:
                continue
            d = [a]*lo + [F(1)]*hi
            got = B(d)
            want = min(a + hi*(1-a), 1 + lo*(1-a))
            if got != want:
                bad2.append((a, lo, hi, got, want))
print(f"    violations: {bad2[:3] or 'none'}")
assert not bad2
print("    exact 5/4 and 4/3 low-dispersion values:")
for (a, lo, hi) in [(F(1,2), 1, 1), (F(2,3), 1, 1), (F(1,2), 1, 2), (F(3,4), 1, 2)]:
    d = [a]*lo + [F(1)]*hi
    print(f"      d = {[str(x) for x in d]}: B/D = {B(d)} = {float(B(d)):.6f}")

print("=" * 74)
print("[B] the deque-order class")
print("    claim: sigma is a deque order of pi  <=>  every prefix of sigma is")
print("    an interval of pi;  there are exactly 2^(k-1) such orders.")


def deque_orders(k):
    """orders built by repeatedly appending the immediate unused predecessor
    or successor in the identity order"""
    out = set()
    for start in range(k):
        stack = [(start, start, (start,))]
        while stack:
            lo, hi, seq = stack.pop()
            if len(seq) == k:
                out.add(seq)
                continue
            if lo > 0:
                stack.append((lo-1, hi, seq + (lo-1,)))
            if hi < k-1:
                stack.append((lo, hi+1, seq + (hi+1,)))
    return out


for k in range(1, 8):
    dq = deque_orders(k)
    # independent characterisation: every prefix is an interval of identity
    def prefix_interval(sig):
        for j in range(1, len(sig)+1):
            p = sorted(sig[:j])
            if p != list(range(p[0], p[0]+len(p))):
                return False
        return True
    chk = {s for s in permutations(range(k)) if prefix_interval(s)}
    print(f"    k={k}: |deque orders| = {len(dq)}   2^(k-1) = {2**(k-1)}   "
          f"matches interval characterisation: {dq == chk}")
    assert len(dq) == 2**(k-1) and dq == chk

print("\n    every sigma-prefix is a difference of two pi-prefixes, so a")
print("    two-sided D/2 bound on pi-prefixes gives D on sigma-prefixes:")
k = 6
ok = True
for sig in sorted(deque_orders(k)):
    for j in range(1, k+1):
        p = sorted(sig[:j])
        # p = identity prefix ending at p[-1] minus identity prefix ending at p[0]-1
        ok &= (p == list(range(p[0], p[-1]+1)))
print(f"    k=6, all {len(deque_orders(k))} deque orders: {ok}")
assert ok
print("    Liu-Reis with m=2 agents gives constant 1 - 1/(2m-2) = 1/2, i.e. D/2.")

print("=" * 74)
print("[C] the switched-interval prefix-representation norm")
print("    V_tau(a) = |a_k| + sum_j |a_j - a_{j+1}|;  error bound (D/2) V_tau(a)")


def V(a):
    return abs(a[-1]) + sum(abs(a[j]-a[j+1]) for j in range(len(a)-1))


def is_single_interval(a):
    nz = [j for j, v in enumerate(a) if v != 0]
    if not nz:
        return True
    if nz != list(range(nz[0], nz[-1]+1)):
        return False
    return len({a[j] for j in nz}) == 1        # constant sign on the interval


for k in (3, 4, 5, 6):
    rows = [r for r in product([-1, 0, 1], repeat=k) if any(r)]
    v2 = [r for r in rows if V(r) <= 2]
    si = [r for r in rows if is_single_interval(r)]
    print(f"    k={k}: rows with V<=2: {len(v2)};  constant-sign single-interval "
          f"rows: {len(si)};  equal: {set(v2)==set(si)}")
    assert set(v2) == set(si)
    two = [r for r in rows if not is_single_interval(r)]
    print(f"          minimum V over NON-single-interval rows: {min(V(r) for r in two)}")
print("    => V<=2 characterises constant-sign single-interval rows exactly.")
print("    NOTE: the draft says rows with two separated intervals have norm")
print("    'at least four'.  The true minimum is THREE, attained e.g. by")
print(f"    a = (1,0,...,0,1): V = {V((1,0,0,1))} for k=4.  The characterisation")
print("    of V<=2 is unaffected, but the constant 4 should read 3.")

# verify the prefix-representation identity itself
print("\n    identity  a = sum_j c_j chi(prefix_j)  with  c_j = a_j - a_{j+1},")
print("    c_k = a_k, and sum_j |c_j| = V_tau(a):")
rng = random.Random(3)
for _ in range(500):
    k = rng.randint(1, 8)
    a = [rng.randint(-2, 2) for _ in range(k)]
    c = [a[j]-a[j+1] for j in range(k-1)] + [a[-1]]
    rec = [sum(c[j] for j in range(i, k)) for i in range(k)]
    assert rec == a
    assert sum(abs(x) for x in c) == V(a)
print("    500 random rows: identity and norm both confirmed")

print("\nP2_V3_UPPER_OK")
