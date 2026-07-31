"""Draw the ladder instances k = 5,6,7,8 as digraphs and as interval families.

For each k the script (i) loads the intervals and, where a concrete rational
member exists, the demands and shares; (ii) VERIFIES the instance from scratch
-- common-point property, exactly two simple s->t_i paths from the raw arc
list, the two chains, and the exact critical constant; (iii) emits TikZ.
"""
import json
from fractions import Fraction as F
from itertools import product

C = '../../certificates/'


def fr(s):
    s = str(s)
    return F(*map(int, s.split('/'))) if '/' in s else F(int(s))


def load(k):
    if k == 5:
        d = json.load(open(C + 'k5_core_fixed_raw_exact.json'))
        return ([tuple(t) for t in d['intervals']],
                [fr(x) for x in d['demands']], [fr(x) for x in d['shares']],
                fr(d['critical']), None)
    if k == 6:
        d = json.load(open(C + 'commonpoint_k6_pell_family_exact.json'))
        m = d['selected_raw_member']
        return ([tuple(t) for t in d['intervals']],
                [fr(x) for x in m['demands']],
                [fr(x) for x in m['late_shares']],
                fr(m['exact_lower_bound']), r'$11-4\sqrt6$')
    if k == 7:
        d = json.load(open(C + 'commonpoint_k7_pell_family_exact.json'))
        m = d['selected_raw_member']
        return ([tuple(t) for t in d['intervals']],
                [fr(x) for x in m['demands']], [fr(x) for x in m['shares']],
                fr(m['exact_lower_bound']), r'$\tfrac{13}{2}-2\sqrt7$')
    if k == 8:
        d = json.load(open(C + 'commonpoint_k8_stationary_envelope_algebraic.json'))
        return [tuple(t) for t in d['intervals']], None, None, None, r'$\approx1.2148$'
    raise KeyError(k)


def simple_paths(arcs, s, tgt):
    adj = {}
    for (u, v) in arcs:
        adj.setdefault(u, []).append(v)
    out, stack = [], [(s, [s])]
    while stack:
        n, path = stack.pop()
        if n == tgt:
            out.append(tuple(zip(path, path[1:])))
            continue
        for w in adj.get(n, []):
            if w not in path:
                stack.append((w, path + [w]))
    return out


def verify(k, IV, d, f, claimed):
    NP = max(b for a, b in IV) + 1
    lo, hi = max(a for a, b in IV), min(b for a, b in IV)
    assert lo <= hi, (k, 'not common-point')
    arcs = [(f"v{j}", f"v{j+1}") for j in range(NP)]
    for i, (a, b) in enumerate(IV):
        arcs += [(f"v{a}", f"t{i}"), (f"v{b+1}", f"t{i}")]
    cnt = [len(simple_paths(arcs, "v0", f"t{i}")) for i in range(len(IV))]
    assert cnt == [2]*len(IV), (k, cnt)
    rows = [[i for i in range(len(IV)) if IV[i][0] <= e <= IV[i][1]]
            for e in range(NP)]
    assert all(set(rows[e]) <= set(rows[e+1]) for e in range(lo))
    assert all(set(rows[e]) >= set(rows[e+1]) for e in range(hi, NP-1))
    msg = (f"  k={k}: common arcs [{lo},{hi}], {NP} spine arcs, "
           f"{len(arcs)} arcs, paths {set(cnt)}, rows {NP}")
    if d is None:
        print(msg + "  (structure only)")
        return lo, hi, NP, len(arcs), None
    S = sum(f)
    m = 0
    while m - 1 + F(1) < S or (m - 1 >= S):
        m += 1
    best = None
    for z in product([0, 1], repeat=len(IV)):
        if sum(z) < m:
            continue
        e = [d[i]*(z[i]-f[i]) for i in range(len(IV))]
        mx = max(max(sum(e[i] for i in r) for r in rows),
                 max(d[i]*(1-f[i]) if z[i] else d[i]*f[i]
                     for i in range(len(IV))))
        best = mx if best is None else min(best, mx)
    ok = (best == claimed)
    print(msg + f", sum f={float(S):.6f}, m={m}, C*={best} "
          f"({float(best):.9f})  matches certificate: {ok}")
    assert ok, (k, best, claimed)
    return lo, hi, NP, len(arcs), m


CB, CR, CG = 'cblue', 'cred', 'cgreen'
OUT = []
for k in (5, 6, 7, 8):
    IV, d, f, claimed, limit = load(k)
    lo, hi, NP, na, m = verify(k, IV, d, f, claimed)
    xs = 10.6/max(NP + 1, 12)          # keep every panel the same total width
    L = [r"\newcommand{\FigLadder%s}{%%" % 'FGHI'[k-5],
         r"\begin{tikzpicture}[x=%.3fcm,y=0.62cm,>=stealth," % xs,
         r"    sp/.style={circle,draw=black!70,fill=black!4,inner sep=0pt,"
         r"minimum size=3.6mm,font=\tiny},",
         r"    tm/.style={circle,draw=cred!80,fill=cred!8,inner sep=0pt,"
         r"minimum size=4.0mm,font=\tiny},",
         r"    spine/.style={->,line width=0.85pt,draw=black!75},",
         r"    early/.style={->,line width=0.6pt,draw=cred!85},",
         r"    late/.style={->,line width=0.6pt,draw=cblue!85}]"]
    # --- digraph
    L.append(r"  \fill[cred!12] (%.2f,-0.34) rectangle (%.2f,0.34);"
             % (lo + .06, hi + .94))
    for v in range(NP + 1):
        L.append(r"  \node[sp] (v%d) at (%d,0) {};" % (v, v))
    for e in range(NP):
        L.append(r"  \draw[spine] (v%d) -- (v%d);" % (e, e + 1))
    for v in (0, NP):
        L.append(r"  \node[font=\tiny,black!60] at (%d,-0.62) {$v_{%d}$};"
                 % (v, v))
    L.append(r"  \node[font=\tiny,cred!85] at (%.2f,0.68) {$e^\ast$};"
             % ((lo + hi + 1)/2))
    order = sorted(range(len(IV)), key=lambda i: IV[i][1] - IV[i][0])
    seen = {}
    for slot, i in enumerate(order):
        a, b = IV[i]
        sgn = 1 if slot % 2 == 0 else -1
        y = sgn * (1.55 + 1.15 * (slot // 2))
        dx = 0.62 * seen.get(a, 0)          # separate terminals sharing a left end
        seen[a] = seen.get(a, 0) + 1
        L.append(r"  \node[tm] (t%d) at (%.2f,%.2f) {$%d$};" % (i, a - dx, y, i))
        L.append(r"  \draw[early] (v%d) to[out=%d,in=%d] (t%d);"
                 % (a, (90 if sgn > 0 else -90) + (25 if dx else 0),
                    (-90 if sgn > 0 else 90), i))
        L.append(r"  \draw[late] (v%d) to[out=%d,in=%d] (t%d);"
                 % (b + 1, 103*sgn, -15*sgn, i))
    # --- interval family, drawn below the digraph
    base = -(1.55 + 1.15 * ((len(IV)-1)//2)) - 2.1
    L.append(r"  \fill[cred!10] (%.2f,%.2f) rectangle (%.2f,%.2f);"
             % (lo, base + .55, hi + 1, base - len(IV) - .35))
    L.append(r"  \draw[->,line width=0.8pt,black!70] (-0.1,%.2f) -- (%.2f,%.2f);"
             % (base + .25, NP + .4, base + .25))
    for e in (0, NP):
        L.append(r"  \node[font=\tiny,black!60] at (%d,%.2f) {%d};"
                 % (e, base + .62, e))
    for i, (a, b) in enumerate(IV):
        y = base - i - .55
        L.append(r"  \fill[cblue!20,rounded corners=1.2pt] (%.2f,%.2f) "
                 r"rectangle (%.2f,%.2f);" % (a + .07, y - .27, b + .93, y + .27))
        L.append(r"  \draw[cblue!70,rounded corners=1.2pt,line width=0.5pt] "
                 r"(%.2f,%.2f) rectangle (%.2f,%.2f);"
                 % (a + .07, y - .27, b + .93, y + .27))
        L.append(r"  \node[font=\tiny,cblue!85] at (%.2f,%.2f) {$[%d,%d]$};"
                 % ((a + b + 1)/2, y, a, b))
        L.append(r"  \node[font=\tiny,anchor=east,black!70] at (%.2f,%.2f) "
                 r"{$t_{%d}$};" % (a - .06, y, i))
        if d is not None:
            L.append(r"  \node[font=\tiny,anchor=west,black!55] at (%.2f,%.2f) "
                     r"{$d=%.3f$, $f=%.3f$};"
                     % (NP + .55, y, float(d[i]), float(f[i])))
    L.append(r"\end{tikzpicture}}")
    OUT.append(chr(10).join(L))
    print(f"     -> FigLadder{'FGHI'[k-5]}")

open('fig_ladder_instances.tex', 'w', encoding='utf-8').write(
    "% ---- ladder instances k = 5,6,7,8 ----" + chr(10)
    + (chr(10)*2).join(OUT) + chr(10))
print('wrote fig_ladder_instances.tex')
