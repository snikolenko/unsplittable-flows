"""Independent raw-arc-list verification of every concrete instance in the
writeup: the public 16/15 seven-vertex gadget, the 9/8 triangle family, the
8867/7800 four-interval point, and the integer four-interval subfamily.

Nothing is assumed about which paths exist: every simple s-t_i path is found by
DFS from the raw arc list, arc loads are rebuilt from path amounts, and all
arithmetic is exact (fractions.Fraction).
"""
from fractions import Fraction as F
from itertools import product


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


def analyse(arcs, s, terms, demands, path_amounts, cost, label,
            require_paths=None):
    """path_amounts[i] : dict  path-tuple -> amount   (must sum to demands[i])"""
    paths = {i: simple_paths(arcs, s, tt) for i, tt in enumerate(terms)}
    counts = [len(paths[i]) for i in range(len(terms))]
    if require_paths is not None:
        assert counts == require_paths, (counts, require_paths)
    # rebuild fractional loads from the supplied path decomposition
    x = {a: F(0) for a in arcs}
    for i, amounts in enumerate(path_amounts):
        assert sum(amounts.values()) == demands[i], (i, sum(amounts.values()))
        for p, amt in amounts.items():
            assert p in paths[i], (i, p, "not a discovered simple path")
            for a in p:
                x[a] += amt
    cx = sum(x[a]*cost.get(a, 0) for a in arcs)
    D = max(demands)
    best, wit = None, []
    ngood = 0
    for choice in product(*[range(c) for c in counts]):
        y = {a: F(0) for a in arcs}
        for i, j in enumerate(choice):
            for a in paths[i][j]:
                y[a] += demands[i]
        cy = sum(y[a]*cost.get(a, 0) for a in arcs)
        if cy > cx:
            continue
        ngood += 1
        over = max(y[a] - x[a] for a in arcs)
        arg = max(arcs, key=lambda a: y[a] - x[a])
        if best is None or over < best:
            best, wit = over, [(choice, arg)]
        elif over == best:
            wit.append((choice, arg))
    print(f"  {label}")
    print(f"    raw simple-path counts   : {counts}")
    print(f"    total routings           : {1}"
          f" x ".join([""]+[str(c) for c in counts])[3:], end="")
    print(f"  = {eval('*'.join(map(str, counts)))}")
    print(f"    c^T x                    : {cx}")
    print(f"    cost-good routings       : {ngood}")
    print(f"    d_max                    : {D}")
    print(f"    exact critical constant  : {best}/{D} = {F(best, D)}"
          f" = {float(F(best,D)):.12f}")
    print(f"    tight routings (choice, witness arc): {wit}")
    return F(best, D)


print("=" * 74)
print("[A] public seven-vertex instance   (claim 16/15)")
arcs = [("s","t1"),("s","t2"),("s","u"),("u","t3"),("u","v"),
        ("v","t1"),("v","w"),("w","t2"),("w","t3")]
cost = {("s","t1"): 2, ("s","t2"): 3, ("u","t3"): 2}
d = [15, 10, 15]
P = [{(("s","t1"),): F(10), (("s","u"),("u","v"),("v","t1")): F(5)},
     {(("s","t2"),): F(6),  (("s","u"),("u","v"),("v","w"),("w","t2")): F(4)},
     {(("s","u"),("u","t3")): F(10),
      (("s","u"),("u","v"),("v","w"),("w","t3")): F(5)}]
got = analyse(arcs, "s", ["t1","t2","t3"], d, P, cost, "seven-vertex gadget",
              require_paths=[2, 2, 2])
assert got == F(16, 15)

print("=" * 74)
print("[B] triangle family  C_n = (3n-1)^2/(8n^2) -> 9/8")
for n in (2, 6, 68, 97, 1000):
    D = 8*n*n
    d = [D, 2*n*(3*n-1), D]
    f1 = F(n+1, 4*n); f2 = F(n+1, 2*n)
    cost = {("s","t1"): 3*n-1, ("s","t2"): 4*n, ("u","t3"): 3*n-1}
    P = [{(("s","t1"),): d[0]*(1-f1),
          (("s","u"),("u","v"),("v","t1")): d[0]*f1},
         {(("s","t2"),): d[1]*(1-f2),
          (("s","u"),("u","v"),("v","w"),("w","t2")): d[1]*f2},
         {(("s","u"),("u","t3")): d[2]*(1-f1),
          (("s","u"),("u","v"),("v","w"),("w","t3")): d[2]*f1}]
    got = analyse(arcs, "s", ["t1","t2","t3"], d, P, cost,
                  f"triangle family n={n}", require_paths=[2, 2, 2])
    assert got == F((3*n-1)**2, 8*n*n), (n, got)
    print(f"    matches (3n-1)^2/(8n^2)  : True   (9/8 - {F(9,8)-got})")

print("=" * 74)
print("[C] four-interval instance, exact rational point   (claim 8867/7800)")
spine = [(f"v{j}", f"v{j+1}") for j in range(5)]
ex = {0: (("v0","t0"), ("v3","t0")), 1: (("v0","t1"), ("v5","t1")),
      2: (("v1","t2"), ("v5","t2")), 3: (("v2","t3"), ("v4","t3"))}
arcs4 = list(spine) + [a for pr in ex.values() for a in pr]
LATE_VTX = {0: 3, 1: 5, 2: 5, 3: 4}
EARLY_VTX = {0: 0, 1: 0, 2: 1, 3: 2}


def four_interval_paths(i):
    early = tuple(spine[:EARLY_VTX[i]]) + (ex[i][0],)
    late = tuple(spine[:LATE_VTX[i]]) + (ex[i][1],)
    return early, late


d = [7800, 5772, 6825, 7800]
u = [2028, 2664, 1001, 1040]
cu = [259, 350, 296, 259]
cost = {ex[i][0]: cu[i] for i in range(4)}
P = []
for i in range(4):
    early, late = four_interval_paths(i)
    P.append({early: F(d[i]-u[i]), late: F(u[i])})
got = analyse(arcs4, "v0", [f"t{i}" for i in range(4)], d, P, cost,
              "four-interval 8867/7800", require_paths=[2, 2, 2, 2])
assert got == F(8867, 7800), got
print(f"    9/8 + 23/1950            : {F(9,8)+F(23,1950)}  match={got==F(9,8)+F(23,1950)}")

print("=" * 74)
print("[D] four-interval integer subfamily   (claim 182359/160000 - 1/n)")
for n in (67, 68, 100, 1000, 10**6):
    d = [160000*n, 115600*n, 136000*n, 160000*n]
    u = [44400*n, 48841*n, 20400*n, 24000*n + 160000]
    cu = [289, 400, 340, 289]
    cost = {ex[i][0]: cu[i] for i in range(4)}
    P = []
    for i in range(4):
        early, late = four_interval_paths(i)
        P.append({early: F(d[i]-u[i]), late: F(u[i])})
    got = analyse(arcs4, "v0", [f"t{i}" for i in range(4)], d, P, cost,
                  f"integer subfamily n={n}", require_paths=[2, 2, 2, 2])
    want = F(182359, 160000) - F(1, n)
    assert got == want, (n, got, want)
    print(f"    == 182359/160000 - 1/n   : True")
    print(f"    > 9/8                    : {got > F(9,8)}"
          f"   (excess {got - F(9,8)})")
    W = d[0]*289
    assert all(d[i]*cu[i] == W for i in range(4))

print("=" * 74)
print("[E] firewall checks on the integer instance")
d = [160000, 115600, 136000, 160000]
dist = sorted(set(d))
print("    distinct demands         :", dist, " ratios 400:289:340:400")
div = [(a, b) for a in dist for b in dist if a < b and b % a == 0]
print("    divisibility relations   :", div, "(chain firewall respected)" if not div else "FAIL")
assert not div
# planarity + K4 minor of the underlying undirected graph
import networkx as nx
G = nx.Graph()
G.add_edges_from([(a, b) for (a, b) in arcs4])
print("    planar (undirected)      :", nx.check_planarity(G)[0])
branch = [{"v0"}, {"v1"}, {"v2","v3","t0"}, {"v4","v5","t1","t2"}]
for S in branch:
    assert nx.is_connected(G.subgraph(S)), S
adj = all(any(G.has_edge(a, b) for a in S for b in T)
          for i, S in enumerate(branch) for T in branch[i+1:])
print("    K4 branch sets connected + pairwise adjacent:", adj, "=> K4 minor, not series-parallel")
assert adj
print("    is series-parallel        : False (K4 minor)")
print("\nV2_INSTANCES_OK")
