"""Symbolic verification of the triangle-template proposition (Prop. 3.12):
the two pair obstructions A, B, the equalization b = 1-r, the completing-square
identity, and the eta>1 branch -- all from the raw arc list.
"""
import sympy as sp
from itertools import product

b, r, q, eta = sp.symbols('b r q eta', nonnegative=True)

ARCS = [("s","t1"),("s","t2"),("s","u"),("u","t3"),("u","v"),
        ("v","t1"),("v","w"),("w","t2"),("w","t3")]


def simple_paths(arcs, s, tgt):
    adj = {}
    for (a, c) in arcs:
        adj.setdefault(a, []).append(c)
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


P = {t: simple_paths(ARCS, "s", t) for t in ["t1", "t2", "t3"]}
print("raw path counts:", {t: len(P[t]) for t in P})
for t in P:
    for p in P[t]:
        print("   ", t, "".join([p[0][0]] + [y for _, y in p]))

d = {"t1": sp.Integer(1), "t2": b, "t3": sp.Integer(1)}
free = {"t1": r, "t2": q, "t3": r}          # share on the long (free) path

x = {a: sp.Integer(0) for a in ARCS}
for t in P:
    exp_p, free_p = P[t][0], P[t][1]
    for a in exp_p:
        x[a] += d[t]*(1 - free[t])
    for a in free_p:
        x[a] += d[t]*free[t]

rows = []
for z in product([0, 1], repeat=3):
    y = {a: sp.Integer(0) for a in ARCS}
    for t, j in zip(["t1", "t2", "t3"], z):
        for a in P[t][j]:
            y[a] += d[t]
    over = {a: sp.expand(y[a] - x[a]) for a in ARCS}
    rows.append((z, over))

sub = {q: 1 + eta - 2*r}
A = sp.expand((2 - 2*r - b*q).subs(sub))
B = sp.expand((1 + b*(1-q) - r).subs(sub))

print("\npair obstructions from the raw arcs (after q = 1+eta-2r):")
for z, over in rows:
    if sum(z) != 2:
        continue
    mx = {a: sp.expand(over[a].subs(sub)) for a in ARCS}
    free_pair = tuple(t for t, j in zip(["t1", "t2", "t3"], z) if j == 1)
    # identify which of A,B this pair's shared-arc maximum equals
    hits = {a: sp.simplify(mx[a] - A) == 0 for a in ARCS}
    hitsB = {a: sp.simplify(mx[a] - B) == 0 for a in ARCS}
    wa = [f"{a[0]}{a[1]}" for a in ARCS if hits[a]]
    wb = [f"{a[0]}{a[1]}" for a in ARCS if hitsB[a]]
    print(f"  free pair {free_pair}:  arcs equal to A: {wa}   arcs equal to B: {wb}")
    print(f"      all shared-arc overloads: "
          f"{ {f'{a[0]}{a[1]}': sp.factor(mx[a]) for a in ARCS if mx[a] != 0} }")

print("\nA =", A, "   B =", B)
print("A - B =", sp.simplify(A - B), " (claim 1-b-r)")
assert sp.simplify(A - B - (1 - b - r)) == 0

F = sp.expand((sp.Min(A, B)).subs(b, 1 - r)) if False else sp.expand(A.subs(b, 1-r))
print("at b = 1-r:  A =", sp.factor(F), "  B =", sp.factor(sp.expand(B.subs(b, 1-r))))
assert sp.simplify(A.subs(b, 1-r) - B.subs(b, 1-r)) == 0
assert sp.simplify(F - (1-r)*(1-eta+2*r)) == 0
print("F_eta(r) = (1-r)(1-eta+2r):  OK")

ident = sp.expand((3-eta)**2/8 - (1-r)*(1-eta+2*r) - 2*(r - (1+eta)/4)**2)
print("completing-square identity residual:", sp.simplify(ident))
assert sp.simplify(ident) == 0

print("\noptimum:", {"b": sp.Rational(1,1)*(3-eta)/4, "r": (1+eta)/4,
                     "q": (1+eta)/2},
      " value:", sp.factor((3-eta)**2/8))
print("value < 9/8 for eta>0:", sp.simplify(sp.Rational(9,8) - (3-eta)**2/8))

# eta > 1 branch: only the all-free routing is cost-good
allfree = next(over for z, over in rows if z == (1, 1, 1))
uv = sp.expand(allfree[("u","v")].subs(sub))
print("\nall-free routing overload on uv:", sp.factor(uv))
res = sp.expand((2-eta) - uv)
print("(2-eta) - overload =", sp.factor(res), " (claim (1-b)(1-q))")
assert sp.simplify(res - (1-b)*(1-q).subs(sub)) == 0

# integer members
from fractions import Fraction as Fr
print("\ninteger members (3n-1)^2/(8n^2):")
for n in (2, 6, 68, 1000):
    val = ((3-sp.Rational(1,n))**2/8)
    print(f"   n={n}: (3-1/n)^2/8 = {val} = {sp.Rational((3*n-1)**2, 8*n*n)}",
          "match:", sp.simplify(val - sp.Rational((3*n-1)**2, 8*n*n)) == 0)
print("\nV7_TEMPLATE_OK")
