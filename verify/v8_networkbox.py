"""Independently re-verify the seven-terminal signed-tree-network lower-box
certificate (claimed deficit 31/250), both abstractly and as a raw flow DAG.
"""
import json
from fractions import Fraction as F
from itertools import product
import sympy as sp
from sympy.solvers.simplex import lpmin, lpmax

import os
HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.path.join(HERE, "..", "certificates", "network_lower_box_exact.json")
cert = json.load(open(CERT))
parent = cert["parent"]
assign = cert["assignment"]
d = [F(s) for s in cert["demands"]]
u = [F(s) for s in cert["amounts"]]
k = len(d)
D = max(d)
f = [u[i]/d[i] for i in range(k)]
print("tree parent array :", parent)
print("endpoint pairs    :", assign)
print("demands           :", [str(x) for x in d], "  D =", D)
print("shares f          :", [str(x) for x in f])

# ---- rebuild the signed network matrix from the tree, independently --------
n = len(parent)
root = parent.index(-1)


def anc_path(v):
    out = []
    while parent[v] != -1:
        out.append(v)          # edge (parent[v], v) identified by its lower end
        v = parent[v]
    return out


def signed_column(qv, pv):
    """signed incidence of the undirected tree path from qv to pv, oriented
    towards pv: +1 on edges of the pv-branch, -1 on edges of the qv-branch."""
    A, B = anc_path(qv), anc_path(pv)
    sA, sB = set(A), set(B)
    col = {}
    for e in B:
        if e not in sA:
            col[e] = col.get(e, 0) + 1
    for e in A:
        if e not in sB:
            col[e] = col.get(e, 0) - 1
    return col


EDGES = [v for v in range(n) if parent[v] != -1]
M = [[0]*k for _ in EDGES]
for i, pair in enumerate(assign):
    qv, pv = pair["minus_endpoint"], pair["plus_endpoint"]
    col = signed_column(qv, pv)
    for e, s in col.items():
        M[EDGES.index(e)][i] = s
print("\nreconstructed signed matrix (rows = tree edges):")
for r, e in zip(M, EDGES):
    print("   edge(par=%2d -> %2d):" % (parent[e], e), r)
claim = cert["matrix"]
print("matches the stored matrix:", M == claim)

# ---- lower box, abstractly -------------------------------------------------
def err(z, row):
    return sum(row[i]*d[i]*(F(z[i]) - f[i]) for i in range(k))


U = [z for z in product([0, 1], repeat=k) if all(err(z, r) <= D for r in M)]
print("\nbinary routings:", 2**k, " upper-good:", len(U),
      " (certificate says", cert["upper_good_routings"], ")")

lam = sp.symbols(f'l0:{len(U)}')
dl = sp.Symbol('delta')
cons = [sp.Eq(sum(lam), 1)] + [l >= 0 for l in lam]
for r in M:
    cons.append(sum(lam[j]*sp.Rational(err(U[j], r)) for j in range(len(U)))
                >= -sp.Rational(D) + dl*0 - dl + dl)   # placeholder, replaced below
cons = [sp.Eq(sum(lam), 1)] + [l >= 0 for l in lam]
for r in M:
    cons.append(sum(lam[j]*sp.Rational(err(U[j], r)) for j in range(len(U)))
                >= -dl)
val, _ = lpmin(dl, cons)
val = sp.nsimplify(val)
print("primal deficit delta =", val, "=", float(val),
      "  matches 31/250:", val == sp.Rational(31, 250))

wv = sp.symbols(f'w0:{len(M)}')
th = sp.Symbol('theta')
cons2 = [sp.Eq(sum(wv), 1)] + [w >= 0 for w in wv]
for z in U:
    cons2.append(sum(wv[e]*sp.Rational(err(z, M[e])) for e in range(len(M)))
                 >= -th)
val2, sol2 = lpmax(-th + 0*sum(wv), cons2 + [th >= 0])
val2 = sp.nsimplify(-val2) if False else sp.nsimplify(val2)
print("dual separator value =", val2, "  agrees:", val2 == val)
print("dual weights         :",
      [str(sol2[wv[e]]) for e in range(len(M))])
print("certificate weights  :", cert["dual_row_weights"])
print("(LB) holds here (delta <= D):", val <= D)

# ---- realize as a raw DAG and re-check every simple path -------------------
arcs = []
for v in range(n):
    if parent[v] != -1:
        arcs.append((f"n{parent[v]}", f"n{v}"))
for i, pair in enumerate(assign):
    qv, pv = pair["minus_endpoint"], pair["plus_endpoint"]
    arcs.append((f"n{qv}", f"t{i}"))
    arcs.append((f"n{pv}", f"t{i}"))


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
    return out


counts = [len(simple_paths(arcs, f"n{root}", f"t{i}")) for i in range(k)]
print("\nraw DAG: arcs =", len(arcs), " simple-path counts =", counts,
      " (certificate says", cert["raw_path_counts"], ")")
print("\nV8_NETWORKBOX_DONE")
