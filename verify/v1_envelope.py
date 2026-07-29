"""Independent symbolic verification of the four-interval envelope.

Written from the prose description in writeup/ssuf_lower_bound_and_barriers.tex
only; builds the raw arc list, discovers paths by DFS, and does every algebraic
step with sympy over Q[t,eps].
"""
import sympy as sp
from itertools import product

t, e = sp.symbols('t epsilon')

# ---------------------------------------------------------------- raw digraph
SPINE = [(f"v{j}", f"v{j+1}") for j in range(5)]
EXITS = {0: (("v0", "t0"), ("v3", "t0")),
         1: (("v0", "t1"), ("v5", "t1")),
         2: (("v1", "t2"), ("v5", "t2")),
         3: (("v2", "t3"), ("v4", "t3"))}
ARCS = list(SPINE) + [a for pair in EXITS.values() for a in pair]


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


PATHS = {}
for i in range(4):
    ps = simple_paths(ARCS, "v0", f"t{i}")
    assert len(ps) == 2, (i, ps)
    PATHS[i] = tuple(sorted(ps, key=len))          # (early, late)
print("[1] blind DFS path counts:", [len(simple_paths(ARCS, 'v0', f't{i}')) for i in range(4)])

# interval = spine arcs on late path minus spine arcs on early path
SPINE_IDX = {a: j for j, a in enumerate(SPINE)}
I = []
for i in range(4):
    early, late = PATHS[i]
    Ie = {SPINE_IDX[a] for a in early if a in SPINE_IDX}
    Il = {SPINE_IDX[a] for a in late if a in SPINE_IDX}
    assert Ie <= Il
    I.append(Il - Ie)
print("[2] intervals derived from raw paths:", [sorted(s) for s in I])
assert [sorted(s) for s in I] == [[0, 1, 2], [0, 1, 2, 3, 4], [1, 2, 3, 4], [2, 3]]

SUPPORT = [sorted(i for i in range(4) if a in I[i]) for a in range(5)]
print("[3] spine row supports:", SUPPORT)

# ---------------------------------------------------------- family parameters
d = [sp.Integer(1), (1 - t)**2, 1 - t, sp.Integer(1)]
f = [2*t - t**2, 1 - 4*t + t**2, t, t + e]
C = 1 + 2*t - 8*t**2 + 6*t**3 - t**4

assert sp.expand(sum(f) - (1 + e)) == 0
print("[4] sum of shares = 1 + eps: OK")

# parameter rectangle:  t = 7/50 + s/50,  eps = tau/68,  (s,tau) in [0,1]^2
s_, tau = sp.symbols('s tau', nonnegative=True)
BOX = {t: sp.Rational(7, 50) + s_/50, e: tau/sp.Integer(68)}


def box_nonneg(expr, name):
    """Certify expr >= 0 on the rectangle via a Bernstein/Handelman certificate:
    expand in the basis s^i (1-s)^j tau^k (1-tau)^l and check all coefficients
    are >= 0.  Equivalent to: substitute s = a/(a+b), tau = c/(c+dd) and check
    the homogenised numerator has nonnegative coefficients."""
    a, b, cc, dd = sp.symbols('a b c dd', nonnegative=True)
    p = sp.expand(expr.subs(BOX))
    degs = sp.Poly(p, s_, tau).total_degree()
    ds = sp.Poly(p, s_).degree() if p.has(s_) else 0
    dt = sp.Poly(p, tau).degree() if p.has(tau) else 0
    hom = sp.expand(sp.together(p.subs({s_: a/(a+b), tau: cc/(cc+dd)}))
                    * (a+b)**ds * (cc+dd)**dt)
    hom = sp.expand(sp.simplify(hom))
    poly = sp.Poly(hom, a, b, cc, dd)
    bad = [(m, c) for m, c in zip(poly.monoms(), poly.coeffs()) if c < 0]
    ok = not bad
    print(f"    {name:<28} Bernstein-nonneg: {ok}" + ("" if ok else f"  bad={bad[:3]}"))
    return ok


def spine_over(z, arc):
    return sum(d[i]*(sp.Integer(z[i]) - f[i]) for i in range(4) if arc in I[i])


def private_over(z, i):
    # chosen exit arc gains the whole demand; it carried d_i*f_i (late) or
    # d_i*(1-f_i) (early) fractionally
    return d[i]*(1 - f[i]) if z[i] else d[i]*f[i]


PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
CLAIM = {(0, 1): (C, 0), (0, 2): (C, 1), (0, 3): (C - e, 2),
         (1, 2): (C, 4), (1, 3): (C - e, 3), (2, 3): (C + t*(1-t) - e, 3)}

print("\n[5] exact pair overload vectors (all five spine arcs):")
for p in PAIRS:
    z = [1 if i in p else 0 for i in range(4)]
    row = [sp.factor(sp.expand(spine_over(z, a))) for a in range(5)]
    print(f"    late {p}: " + " | ".join(str(r) for r in row))

print("\n[6] witness identities and rectangle dominance certificates:")
all_ok = True
for p, (val, arc) in CLAIM.items():
    z = [1 if i in p else 0 for i in range(4)]
    idt = sp.expand(spine_over(z, arc) - val)
    assert idt == 0, (p, idt)
    print(f"  late {p}: witness arc {arc}, value {sp.factor(val)}  [identity OK]")
    for a in range(5):
        if a == arc:
            continue
        gap = sp.expand(val - spine_over(z, a))
        print(f"      gap to arc {a}: {sp.factor(gap)}")
        all_ok &= box_nonneg(gap, f"pair{p} arc{a}")
    for i in range(4):
        gap = sp.expand(val - private_over(z, i))
        all_ok &= box_nonneg(gap, f"pair{p} priv{i}")
print("[6] all dominance certificates nonnegative on rectangle:", all_ok)

# ------------------------------------------- C(t) - eps > 1 on the rectangle
print("\n[7] C(t) - eps > 1 on the rectangle:")
print("    C(7/50) =", sp.nsimplify(C.subs(t, sp.Rational(7, 50))),
      "=", sp.N(C.subs(t, sp.Rational(7, 50)), 12))
print("    C(4/25) =", sp.nsimplify(C.subs(t, sp.Rational(4, 25))),
      "=", sp.N(C.subs(t, sp.Rational(4, 25)), 12))
margin = sp.expand(C - e - 1 - sp.Rational(1, 10))     # claim C-eps-1 > 1/10
print("    certificate for C - eps - 1 > 1/10:")
box_nonneg(margin, "C-eps-1-1/10")

# ------------------------------------------------------ full 16-routing check
print("\n[8] min over cost-good routings of max all-arc overload:")
for tv in (sp.Rational(7, 50), sp.Rational(3, 20), sp.Rational(4, 25),
           sp.Rational(149, 1000)):
    for ev in (sp.Rational(1, 68), sp.Rational(1, 200), sp.Rational(1, 10**6)):
        sub = {t: tv, e: ev}
        best = None
        for z in product([0, 1], repeat=4):
            if sum(z) < 2:            # cost-good <=> at least two late
                continue
            vals = [spine_over(z, a).subs(sub) for a in range(5)]
            vals += [private_over(z, i).subs(sub) for i in range(4)]
            m = max(vals)
            best = m if best is None else min(best, m)
        want = (C - e).subs(sub)
        assert best == want, (tv, ev, best, want)
    print(f"    t={tv}: min over cost-good routings == C(t)-eps for all tested eps  OK")

# ------------------------------------------------------------------- calculus
print("\n[9] optimization of C(t):")
Cp = sp.expand(sp.diff(C, t))
assert sp.expand(Cp + 2*(t-1)*(2*t**2 - 7*t + 1)) == 0
print("    C'(t) =", sp.factor(Cp))
roots = sp.solve(sp.Eq(Cp, 0), t)
print("    critical points:", roots, "=", [sp.N(r, 12) for r in roots])
tstar = (7 - sp.sqrt(41))/4
assert sp.simplify(Cp.subs(t, tstar)) == 0
assert sp.Rational(7, 50) < sp.N(tstar) < sp.Rational(4, 25)
val = sp.radsimp(sp.expand(C.subs(t, tstar)))
target = (299 - 41*sp.sqrt(41))/32
assert sp.simplify(val - target) == 0
print("    t*      =", tstar, "=", sp.N(tstar, 15))
print("    C(t*)   =", target, "=", sp.N(target, 18))
print("    9/8     =", sp.N(sp.Rational(9, 8), 18))
# second derivative test
print("    C''(t*) =", sp.N(sp.diff(C, t, 2).subs(t, tstar), 8), "(<0 => max)")
# C is strictly concave-with-single-max on the rectangle: show C' > 0 for
# t<t*, <0 for t>t* inside [7/50,4/25]
print("    C'(7/50) =", sp.N(Cp.subs(t, sp.Rational(7, 50)), 8),
      " C'(4/25) =", sp.N(Cp.subs(t, sp.Rational(4, 25)), 8))

# ----------------------------------------------------- integer subfamily t=3/20
print("\n[10] integer subfamily t=3/20:")
C320 = C.subs(t, sp.Rational(3, 20))
print("    C(3/20) =", C320, "=", sp.N(C320, 12), " (claim 182359/160000)")
assert C320 == sp.Rational(182359, 160000)
n = sp.Symbol('n', positive=True, integer=True)
thr = sp.solve(sp.Eq(sp.Rational(182359, 160000) - 1/n, sp.Rational(9, 8)), n)
print("    C_n > 9/8  <=>  n >", thr, "=", [sp.N(x, 10) for x in thr])
C68 = sp.Rational(182359, 160000) - sp.Rational(1, 68)
print("    C_68 =", C68, "= 9/8 +", sp.nsimplify(C68 - sp.Rational(9, 8)))
assert C68 == sp.Rational(3060103, 2720000)
assert C68 - sp.Rational(9, 8) == sp.Rational(103, 2720000)
C67 = sp.Rational(182359, 160000) - sp.Rational(1, 67)
print("    C_67 =", C67, "= 9/8 +", sp.nsimplify(C67 - sp.Rational(9, 8)), "(must be < 9/8)")
assert C67 < sp.Rational(9, 8)
print("\nV1_ENVELOPE_OK")
