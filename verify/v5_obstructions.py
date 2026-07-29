"""Verify the ten-arc DAG obstruction (Proposition 'obs3' in the paper):
a lower-good routing whose every first upper-box repair breaks the lower box,
on an instance that nevertheless satisfies (LB).
"""
from fractions import Fraction as F
from itertools import product

ARCS = [("s","u"),("s","w"),("u","v"),("u","tC"),("v","w"),
        ("v","tA"),("v","tB"),("v","tC"),("w","tA"),("w","tB")]
TERM = ["tA","tB","tC"]
D = F(1)
d = {"tA": F(5,6), "tB": F(1), "tC": F(1,2)}


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
    return sorted(out, key=lambda p: (len(p), p))


P = {t: simple_paths(ARCS, "s", t)[::-1] for t in TERM}   # longest first
print("simple path counts:", {t: len(P[t]) for t in TERM},
      " total routings:", len(P["tA"])*len(P["tB"])*len(P["tC"]))
for t in TERM:
    for p in P[t]:
        print("   ", t, "->", "".join([p[0][0]] + [b for _, b in p]))

shares = {"tA": [F(1,6), F(2,3), F(1,6)],
          "tB": [F(1,2), F(1,3), F(1,6)],
          "tC": [F(35,36), F(1,36)]}
for t in TERM:
    assert len(shares[t]) == len(P[t]), (t, len(shares[t]), len(P[t]))
    assert sum(shares[t]) == 1

q = {a: F(0) for a in ARCS}
for t in TERM:
    for p, sh in zip(P[t], shares[t]):
        for a in p:
            q[a] += d[t]*sh
print("\nq =", {f"{a[0]}{a[1]}": str(q[a]) for a in ARCS})
CLAIM = {("s","u"): F(73,36), ("u","v"): F(145,72), ("v","w"): F(23,36),
         ("s","w"): F(11,36), ("u","tC"): F(1,72), ("v","tA"): F(5,9),
         ("v","tB"): F(1,3), ("v","tC"): F(35,72), ("w","tA"): F(5,18),
         ("w","tB"): F(2,3)}
mismatch = {a: (q[a], CLAIM[a]) for a in ARCS if q[a] != CLAIM[a]}
print("matches the recorded q:", not mismatch, mismatch if mismatch else "")

R = []
for choice in product(*[range(len(P[t])) for t in TERM]):
    y = {a: F(0) for a in ARCS}
    for t, j in zip(TERM, choice):
        for a in P[t][j]:
            y[a] += d[t]
    up = all(y[a] <= q[a] + D for a in ARCS)
    lo = all(y[a] >= q[a] - D for a in ARCS)
    R.append((choice, y, up, lo))
print("\ntotal routings:", len(R),
      " upper-good:", sum(1 for r in R if r[2]),
      " lower-good:", sum(1 for r in R if r[3]),
      " both:", sum(1 for r in R if r[2] and r[3]))

# the distinguished lower-good routing: A and B on their LONG paths, C short.
def path_index(t, name):
    for j, p in enumerate(P[t]):
        if "".join([p[0][0]] + [b for _, b in p]) == name:
            return j
    raise KeyError(name)


start = (path_index("tA", "suvwtA"), path_index("tB", "suvwtB"),
         path_index("tC", "sutC"))
ys = dict(next(r[1] for r in R if r[0] == start))
print("\nstart routing", start,
      " lower-good:", next(r[3] for r in R if r[0] == start),
      " upper-good:", next(r[2] for r in R if r[0] == start))
viol = [a for a in ARCS if ys[a] > q[a] + D]
print("  upper-box violations:", [f"{a[0]}{a[1]}" for a in viol],
      "  excess:", [str(ys[a]-q[a]-D) for a in viol])

print("\n  single-terminal reroutes that repair the upper box:")
any_keeps_lower = False
for ti, t in enumerate(TERM):
    for j in range(len(P[t])):
        if j == start[ti]:
            continue
        ch = list(start); ch[ti] = j
        r = next(x for x in R if x[0] == tuple(ch))
        if r[2]:
            low = r[3]
            any_keeps_lower |= low
            bad = [f"{a[0]}{a[1]}" for a in ARCS if r[1][a] < q[a] - D]
            print(f"    move {t} -> path {j}: upper-good, lower-good={low}"
                  f"  breaks {bad}")
print("  some first repair keeps the lower box:", any_keeps_lower,
      "  <-- the source note claimed this should be False")

# systematic search: is there ANY lower-good, upper-bad routing all of whose
# upper-box-restoring single-terminal reroutes break the lower box?
print("\nsystematic search over all lower-good, upper-bad routings:")
found = []
for (ch, y, up, lo) in R:
    if up or not lo:
        continue
    repairs = []
    for ti, t in enumerate(TERM):
        for j in range(len(P[t])):
            if j == ch[ti]:
                continue
            c2 = list(ch); c2[ti] = j
            r2 = next(x for x in R if x[0] == tuple(c2))
            if r2[2]:
                repairs.append((t, j, r2[3]))
    tag = ("no repair at all" if not repairs else
           ("ALL repairs break the lower box"
            if not any(l for _, _, l in repairs) else
            "some repair keeps the lower box"))
    print(f"  routing {ch}: {len(repairs)} single-terminal repairs -> {tag}")
    if repairs and not any(l for _, _, l in repairs):
        found.append(ch)
print("  routings witnessing the claimed obstruction:", found)

# (LB) itself holds: exhibit the recorded four-term decomposition
import sympy as sp
from sympy.solvers.simplex import lpmin
U = [r for r in R if r[2]]
lam = sp.symbols(f'l0:{len(U)}')
dl = sp.Symbol('delta')
cons = [sp.Eq(sum(lam), 1)] + [l >= 0 for l in lam]
for a in ARCS:
    cons.append(sum(lam[j]*sp.Rational(U[j][1][a]) for j in range(len(U)))
                >= sp.Rational(q[a]) - dl)
val, sol = lpmin(dl, cons)
print("\n  exact lower-box deficit of this instance:", sp.nsimplify(val),
      " (LB) holds:", val <= 1)
print("  support of an optimal mixture:")
for j in range(len(U)):
    if sol[lam[j]] != 0:
        print(f"     lambda={sol[lam[j]]}  routing {U[j][0]}")
print("\nV5_OBSTRUCTIONS_OK")
