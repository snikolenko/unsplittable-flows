"""Generate TikZ digraph figures for Part II.

  fig_cpschema.tex  -- schematic of a general common-point interval system
  fig_parti.tex     -- the Part I four-interval instance as a raw digraph,
                       with fractional arc loads and costs
  fig_k17graph.tex  -- the whole k=17 record digraph (67 arcs)

Exit arcs are nested by interval length so that they do not cross.
"""
import json

OUT = {}
BS = chr(92)

# =====================================================================
# 1. schematic: what a common-point interval system IS
# =====================================================================
sch = [r"""% ---- schematic of a common-point interval system ----
\newcommand{\FigCPSchema}{%
\begin{tikzpicture}[x=1.06cm,y=0.95cm,>=stealth,
    sp/.style={circle,draw=black!70,fill=black!4,inner sep=0pt,
               minimum size=5.2mm,font=\scriptsize},
    tm/.style={circle,draw=cred!80,fill=cred!8,inner sep=0pt,
               minimum size=5.2mm,font=\scriptsize},
    spine/.style={->,line width=1.05pt,draw=black!75},
    early/.style={->,line width=0.75pt,draw=cred!85},
    late/.style={->,line width=0.75pt,draw=cblue!85}]
"""]
NS = 8
# arc j occupies x in [j, j+1]; the four intervals below meet exactly in arc 3
sch.append(r"  \fill[cred!12] (3.06,-0.40) rectangle (3.94,0.40);")
for j in range(NS + 1):
    sch.append(r"  \node[sp] (v%d) at (%d,0) {$v_{%d}$};" % (j, j, j))
for j in range(NS):
    sch.append(r"  \draw[spine] (v%d) -- (v%d);" % (j, j + 1))
sch.append(r"  \node[font=\scriptsize,cred!85] at (3.5,0.72) {common arc $e^\ast=3$};")
demo = [(2, 3), (1, 5), (3, 7), (0, 6)]
assert set(range(2, 4)) & set(range(1, 6)) & set(range(3, 8)) & set(range(0, 7)) == {3}
order = sorted(range(len(demo)), key=lambda n: demo[n][1] - demo[n][0])
for slot, n in enumerate(order):
    a, b = demo[n]
    sgn = 1 if slot % 2 == 0 else -1
    lane = slot // 2
    y = sgn * (1.70 + 1.42 * lane)
    sch.append(r"  \node[tm] (t%d) at (%.2f,%.2f) {$t_{%d}$};" % (n, a, y, n))
    sch.append(r"  \draw[early] (v%d) -- (t%d);" % (a, n))
    sch.append(r"  \draw[late]  (v%d) to[out=%d,in=%d] (t%d);"
               % (b + 1, 100 * sgn, -18 * sgn, n))
    sch.append(r"  \node[font=\tiny,black!60,anchor=west] at (%.2f,%.2f) "
               r"{$I_{%d}=[%d,%d]$};" % (a + 0.30, y, n, a, b))
sch.append(r"""  \node[anchor=west,font=\scriptsize,align=left] at (9.1,1.7)
    {\textcolor{cred}{$\longrightarrow$}\; early exit $v_{\ell_i}t_i$,\\
     \phantom{xx}per-unit cost $1/d_i$\\[3pt]
     \textcolor{cblue}{$\longrightarrow$}\; late exit $v_{r_i+1}t_i$,\\
     \phantom{xx}cost $0$\\[3pt]
     spine arcs: cost $0$};
  \node[anchor=west,font=\scriptsize,align=left,black!65] at (9.1,-1.9)
    {Every $s\to t_i$ path is a spine\\
     prefix plus one exit, so there\\
     are exactly two of them.\\[3pt]
     Switching $t_i$ from early to\\
     late adds $d_i$ on exactly $I_i$.};
\end{tikzpicture}}""")
OUT['fig_cpschema.tex'] = "\n".join(sch)

# =====================================================================
# 2. the Part I four-interval instance, as a real digraph with loads
# =====================================================================
d = [7800, 5772, 6825, 7800]
u = [2028, 2664, 1001, 1040]
IV = [(0, 2), (0, 4), (1, 4), (2, 3)]
NP = 5
spine_load = [0] * NP
for i, (a, b) in enumerate(IV):
    for e in range(0, a):
        spine_load[e] += d[i] - u[i]
    for e in range(0, b + 1):
        spine_load[e] += u[i]
assert spine_load == [19317, 13493, 6733, 4705, 3665], spine_load

g = [r"""% ---- the Part I four-interval instance as a raw digraph ----
\newcommand{\FigPartI}{%
\begin{tikzpicture}[x=1.95cm,y=1.28cm,>=stealth,
    sp/.style={circle,draw=black!70,fill=black!4,inner sep=0pt,
               minimum size=6.4mm,font=\small},
    tm/.style={circle,draw=cred!80,fill=cred!8,inner sep=0pt,
               minimum size=6.4mm,font=\small},
    spine/.style={->,line width=1.15pt,draw=black!75},
    early/.style={->,line width=0.85pt,draw=cred!85},
    late/.style={->,line width=0.85pt,draw=cblue!85}]
"""]
g.append(r"  \fill[cred!12] (2.06,-0.26) rectangle (2.94,0.26);")
for j in range(NP + 1):
    g.append(r"  \node[sp] (v%d) at (%d,0) {$v_{%d}$};" % (j, j, j))
g.append(r"  \node[font=\scriptsize,black!60] at (-0.45,0) {$s=$};")
for j in range(NP):
    g.append(r"  \draw[spine] (v%d) -- node[above,font=\scriptsize,black!65,"
             r"inner sep=1.6pt]{$%d$} (v%d);" % (j, spine_load[j], j + 1))
g.append(r"  \node[font=\scriptsize,cred!85] at (2.5,-0.52) {$e^\ast=2$};")
# nest by interval length: t3 (len 2) innermost, t1 (len 5) outermost
LANE = {3: (1, 1.35, 0.0), 0: (-1, 1.35, 0.30), 2: (1, 2.70, 0.0), 1: (-1, 2.70, -0.55)}
for i, (a, b) in enumerate(IV):
    sgn, hgt, dx = LANE[i]
    y = sgn * hgt
    g.append(r"  \node[tm] (t%d) at (%.2f,%.2f) {$t_{%d}$};" % (i, a + dx, y, i))
    g.append(r"  \draw[early] (v%d) to[out=%d,in=%d] node[left,font=\tiny,"
             r"cred!85,inner sep=1.4pt]{$%d$} (t%d);"
             % (a, (90 if sgn > 0 else -90) + (-25 if dx < 0 else (18 if dx > 0 else 0)),
                (-90 if sgn > 0 else 90), d[i] - u[i], i))
    g.append(r"  \draw[late] (v%d) to[out=%d,in=%d] "
             r"node[%s,font=\tiny,cblue!85,pos=0.52,inner sep=1.4pt]{$%d$} (t%d);"
             % (b + 1, 105 * sgn, -16 * sgn,
                "above" if sgn > 0 else "below", u[i], i))
    g.append(r"  \node[font=\tiny,black!55,anchor=west] at (%.2f,%.2f) "
             r"{$I_{%d}=[%d,%d]$};" % (a + dx + 0.16, y + sgn * 0.32, i, a, b))
g.append(r"""  \node[anchor=west,font=\scriptsize,align=left] at (5.45,1.55)
    {$\vd=(7800,\,5772,\,6825,\,7800)$\\[3pt]
     early per-unit costs $(259,350,296,259)$,\\
     every other arc free; each full early\\
     path then costs $W=2\,020\,200$\\[3pt]
     arc labels are the fractional loads $x_a$};
  \node[anchor=west,font=\scriptsize,align=left,black!65] at (5.45,-1.75)
    {$\vc^{\mathsf T}\vx=6\,057\,492<3W$, so a routing\\
     is cost-good iff at least two terminals\\
     go late: first cost-good rank $m=2$,\\
     complement mass $q=3$.\\[3pt]
     $C^\star=8867/7800=1.1368\ldots$};
\end{tikzpicture}}""")
OUT['fig_parti.tex'] = "\n".join(g)

# =====================================================================
# 3. the k=17 record digraph
# =====================================================================
cert = json.load(open('../../certificates/'
                      'commonpoint_joint_k17_refined_record_exact.json'))
IV17 = [tuple(t) for t in cert['intervals']]
K, N17 = len(IV17), max(b for a, b in IV17) + 1
h = [r"""% ---- the k=17 record digraph ----
\newcommand{\FigKSeventeen}{%
\begin{tikzpicture}[x=0.46cm,y=0.46cm,>=stealth,
    sp/.style={circle,draw=black!65,fill=black!5,inner sep=0pt,
               minimum size=3.0mm},
    tm/.style={circle,draw=cred!75,fill=cred!10,inner sep=0.6pt,
               minimum size=3.6mm,font=\tiny},
    spine/.style={->,line width=0.85pt,draw=black!75},
    early/.style={->,line width=0.45pt,draw=cred!75},
    late/.style={->,line width=0.45pt,draw=cblue!75}]
"""]
h.append(r"  \fill[cred!14] (16.06,-0.5) rectangle (16.94,0.5);")
for j in range(N17 + 1):
    h.append(r"  \node[sp] (v%d) at (%d,0) {};" % (j, j))
for j in range(N17):
    h.append(r"  \draw[spine] (v%d) -- (v%d);" % (j, j + 1))
for j in (0, 8, 16, 24, 33):
    h.append(r"  \node[font=\tiny,black!60] at (%d,-0.85) {$v_{%d}$};" % (j, j))
h.append(r"  \node[font=\tiny,cred!85] at (16.5,1.0) {$e^\ast$};")
# nest by interval length within each side so that late arcs do not cross
lengths = sorted(range(K), key=lambda i: IV17[i][1] - IV17[i][0])
for slot, i in enumerate(lengths):
    a, b = IV17[i]
    sgn = 1 if slot % 2 == 0 else -1
    lane = slot // 2
    y = sgn * (2.3 + 1.35 * lane)
    h.append(r"  \node[tm] (t%d) at (%d,%.2f) {$%d$};" % (i, a, y, i))
    h.append(r"  \draw[early] (v%d) -- (t%d);" % (a, i))
    h.append(r"  \draw[late] (v%d) to[out=%d,in=%d] (t%d);"
             % (b + 1, 104 * sgn, -14 * sgn, i))
h.append(r"""  \node[font=\scriptsize,anchor=west] at (0.4,14.9)
    {\textcolor{cred}{red}: early exit $v_{\ell_i}t_i$, per-unit cost $1/d_i$
     \quad\textcolor{cblue}{blue}: late exit $v_{r_i+1}t_i$, free
     \quad spine arcs: free};
\end{tikzpicture}}""")
OUT['fig_k17graph.tex'] = "\n".join(h)

for name, body in OUT.items():
    open(name, 'w', encoding='utf-8').write(body + "\n")
    print('wrote', name)
