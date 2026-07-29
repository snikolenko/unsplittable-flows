"""Search small DAGs for a genuine 'monotone local repair is impossible'
witness: a routing that is lower-good but not upper-good, every single-terminal
reroute restoring the upper box of which breaks the lower box.
"""
from fractions import Fraction as F
from itertools import product, combinations
import random

random.seed(11)


def simple_paths(arcs, s, tgt):
    adj = {}
    for (a, b) in arcs:
        adj.setdefault(a, []).append(b)
    out, stack = [], [(s, [s])]
    while stack:
        node, path = stack.pop()
        if node == tgt:
            out.append(tuple(zip(path, path[1:])))
            continue
        for x in adj.get(node, []):
            if x not in path:
                stack.append((x, path + [x]))
    return sorted(out, key=lambda p: (-len(p), p))


def analyse(arcs, terms, d, shares, D):
    P = {t: simple_paths(arcs, "s", t) for t in terms}
    if any(len(P[t]) != len(shares[t]) for t in terms):
        return None
    q = {a: F(0) for a in arcs}
    for t in terms:
        for p, sh in zip(P[t], shares[t]):
            for a in p:
                q[a] += d[t]*sh
    R = []
    for ch in product(*[range(len(P[t])) for t in terms]):
        y = {a: F(0) for a in arcs}
        for t, j in zip(terms, ch):
            for a in P[t][j]:
                y[a] += d[t]
        R.append((ch, y,
                  all(y[a] <= q[a] + D for a in arcs),
                  all(y[a] >= q[a] - D for a in arcs)))
    wit = []
    for (ch, y, up, lo) in R:
        if up or not lo:
            continue
        reps = []
        for ti, t in enumerate(terms):
            for j in range(len(P[t])):
                if j == ch[ti]:
                    continue
                c2 = list(ch); c2[ti] = j
                r2 = next(x for x in R if x[0] == tuple(c2))
                if r2[2]:
                    reps.append(r2[3])
        if reps and not any(reps):
            wit.append((ch, len(reps)))
    return P, q, R, wit


# ---- the recorded ten-arc instance, restated unambiguously -----------------
ARCS = [("s","u"),("s","w"),("u","v"),("u","tC"),("v","w"),
        ("v","tA"),("v","tB"),("v","tC"),("w","tA"),("w","tB")]
TERM = ["tA","tB","tC"]
d = {"tA": F(5,6), "tB": F(1), "tC": F(1,2)}
shares = {"tA": [F(1,6), F(2,3), F(1,6)],     # suvwtA, suvtA, swtA
          "tB": [F(1,2), F(1,3), F(1,6)],     # suvwtB, suvtB, swtB
          "tC": [F(35,36), F(1,36)]}          # suvtC, sutC
out = analyse(ARCS, TERM, d, shares, F(1))
P, q, R, wit = out
print("recorded ten-arc instance:")
print("  q =", {f'{a[0]}{a[1]}': str(q[a]) for a in ARCS})
print("  routings:", len(R), " upper-good:", sum(1 for r in R if r[2]),
      " lower-good:", sum(1 for r in R if r[3]))
print("  witnesses of 'every repair breaks the lower box':", wit,
      "(none: the recorded claim does not hold)")

# ---- random search over small layered DAGs --------------------------------
NODES = ["s", "u", "v", "w"]
TERMS3 = ["t1", "t2", "t3"]
CAND = [(a, b) for i, a in enumerate(NODES) for b in NODES[i+1:]]
CANDT = [(a, t) for a in NODES for t in TERMS3]

best = None
trials = 0
for trial in range(40000):
    ne = random.randint(3, len(CAND))
    core = random.sample(CAND, ne)
    exits = []
    for t in TERMS3:
        opts = [e for e in CANDT if e[1] == t]
        exits += random.sample(opts, random.randint(2, 3))
    arcs = core + exits
    if not all(any(a == "s" for a, _ in arcs) for _ in [0]):
        continue
    P = {t: simple_paths(arcs, "s", t) for t in TERMS3}
    if any(len(P[t]) < 2 or len(P[t]) > 4 for t in TERMS3):
        continue
    dd = {t: F(random.randint(4, 12), 12) for t in TERMS3}
    if max(dd.values()) != 1:
        continue
    sh = {}
    ok = True
    for t in TERMS3:
        n = len(P[t])
        cuts = sorted(random.randint(0, 12) for _ in range(n-1))
        vals = [F(b-a, 12) for a, b in zip([0]+cuts, cuts+[12])]
        if sum(vals) != 1:
            ok = False
        sh[t] = vals
    if not ok:
        continue
    trials += 1
    res = analyse(arcs, TERMS3, dd, sh, F(1))
    if res is None:
        continue
    _, qq, RR, ww = res
    if ww:
        best = (arcs, dd, sh, qq, RR, ww)
        break

print(f"\nrandom search: {trials} valid small-DAG instances tested")
if best:
    arcs, dd, sh, qq, RR, ww = best
    print("  FOUND a witness:")
    print("   arcs   =", arcs)
    print("   demands=", {k: str(v) for k, v in dd.items()})
    print("   shares =", {k: [str(x) for x in v] for k, v in sh.items()})
    print("   witness routings:", ww)
else:
    print("  no instance found in which every single-terminal upper-box repair")
    print("  from a lower-good routing breaks the lower box.")
print("\nV6_REPAIR_SEARCH_DONE")
