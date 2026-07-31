"""The instance behind the running Example (deletion stars), drawn twice:
as a raw digraph and as a staircase interval system with the deletion star.

k=6, unit demands, f_i = 2/3, identity and reverse orders; realized by the
staircase I_i = [i, i+5] on a spine of 11 arcs, common arc e* = 5.
"""
from fractions import Fraction as F

k = 6
IV = [(i, i + 5) for i in range(k)]
NP = max(b for a, b in IV) + 1            # 11 spine arcs
CP = max(a for a, b in IV)                # common arc 5
assert min(b for a, b in IV) == CP

# sanity: the rows really are the two chains
rows = [[i for i in range(k) if IV[i][0] <= e <= IV[i][1]] for e in range(NP)]
assert rows[:CP+1] == [list(range(e+1)) for e in range(CP+1)]
assert rows[CP:] == [list(range(e-CP, k)) for e in range(CP, NP)]
print("rows:", rows)

j = 2
S = [i for i in range(k) if i != j]
c = {0: F(0), 1: F(1, 3), 3: F(0), 4: F(1, 3), 5: F(1, 3)}
h = {i: (F(1, 3) - c[i] if i in c else -F(2, 3)) for i in range(k)}
assert sum(c.values()) == 1 and sum(h.values()) == 0
print("credits:", {i: str(v) for i, v in c.items()})
print("residual:", {i: str(v) for i, v in h.items()})

L = [r"""% ---- the running-example instance for deletion stars ----
\newcommand{\FigStarInstance}{%
\begin{tikzpicture}[x=0.95cm,y=0.80cm,>=stealth,
    sp/.style={circle,draw=black!70,fill=black!4,inner sep=0pt,
               minimum size=4.6mm,font=\tiny},
    tm/.style={circle,draw=cred!80,fill=cred!8,inner sep=0pt,
               minimum size=4.8mm,font=\scriptsize},
    om/.style={circle,draw=cgreen!85,fill=cgreen!14,line width=0.9pt,
               inner sep=0pt,minimum size=4.8mm,font=\scriptsize},
    spine/.style={->,line width=1.0pt,draw=black!75},
    early/.style={->,line width=0.7pt,draw=cred!85},
    late/.style={->,line width=0.7pt,draw=cblue!85}]
"""]
# ---- panel A: the digraph
L.append(r"  \fill[cred!12] (%.2f,-0.30) rectangle (%.2f,0.30);" % (CP + .06, CP + .94))
for v in range(NP + 1):
    L.append(r"  \node[sp] (v%d) at (%d,0) {$%d$};" % (v, v, v))
for e in range(NP):
    L.append(r"  \draw[spine] (v%d) -- (v%d);" % (e, e + 1))
L.append(r"  \node[font=\scriptsize,cred!85] at (%.2f,0.62) {$e^\ast=%d$};"
         % (CP + .5, CP))
L.append(r"  \node[font=\scriptsize,black!60] at (-0.75,0) {$s=v_0$};")
for i, (a, b) in enumerate(IV):
    sgn = 1 if i % 2 == 0 else -1
    y = sgn * (1.45 + 0.95 * (i // 2))
    st = 'om' if i == j else 'tm'
    L.append(r"  \node[%s] (t%d) at (%d,%.2f) {$%d$};" % (st, i, a, y, i))
    L.append(r"  \draw[early] (v%d) -- (t%d);" % (a, i))
    L.append(r"  \draw[late] (v%d) to[out=%d,in=%d] (t%d);"
             % (b + 1, 105 * sgn, -16 * sgn, i))
L.append(r"""  \node[font=\scriptsize,anchor=west] at (12.05,0)
    {\textcolor{cred}{$\rightarrow$} early $v_it_i$\ \
     \textcolor{cblue}{$\rightarrow$} late $v_{i+6}t_i$};""")
L.append(r"\end{tikzpicture}}")
open('fig_starinst.tex', 'w', encoding='utf-8').write("\n".join(L) + "\n")
print('wrote fig_starinst.tex')

# ---- panel B: the staircase with the star
M = [r"""% ---- the staircase interval system and the deletion star ----
\newcommand{\FigStarStaircase}{%
\begin{tikzpicture}[x=0.95cm,y=0.62cm,>=stealth]
"""]
M.append(r"  \fill[cred!10] (%.2f,0.65) rectangle (%.2f,%.2f);"
         % (CP, CP + 1, -k - 0.55))
for e in range(NP + 1):
    M.append(r"  \draw[black!18,line width=0.5pt] (%d,0.55) -- (%d,%.2f);"
             % (e, e, -k - 0.35))
M.append(r"  \draw[->,line width=0.9pt,black!70] (-0.15,0.25) -- (%.2f,0.25);"
         % (NP + 0.4))
for e in (0, CP, NP):
    M.append(r"  \node[font=\tiny,black!65] at (%d,0.62) {$%d$};" % (e, e))
for i, (a, b) in enumerate(IV):
    y = -1 - i
    col = 'cgreen' if i == j else 'cblue'
    M.append(r"  \fill[%s!22,rounded corners=1.4pt] (%.2f,%.2f) rectangle "
             r"(%.2f,%.2f);" % (col, a + .07, y - .30, b + .93, y + .30))
    M.append(r"  \draw[%s!75,rounded corners=1.4pt,line width=0.6pt] "
             r"(%.2f,%.2f) rectangle (%.2f,%.2f);"
             % (col, a + .07, y - .30, b + .93, y + .30))
    M.append(r"  \node[font=\tiny,%s!85] at (%.2f,%.2f) {$I_%d=[%d,%d]$};"
             % (col, (a + b + 1) / 2, y, i, a, b))
    lab = (r"\textbf{omitted}" if i == j
           else (r"$c_%d=%s$" % (i, ('0' if c[i] == 0 else r'\tfrac13'))))
    M.append(r"  \node[font=\tiny,anchor=west,black!70] at (%.2f,%.2f) {%s};"
             % (NP + 1.15, y, lab))
    M.append(r"  \node[font=\tiny,anchor=east,%s!85] at (%.2f,%.2f) {$t_%d$};"
             % (col, a - .05, y, i))
M.append(r"""  \draw[->,cgreen!75,line width=0.8pt] (0.1,-7.05) -- (5.4,-7.05);
  \node[font=\tiny,cgreen!70!black,anchor=west] at (0.1,-7.55)
    {rows $e\le5$: prefixes of $\pi$};
  \draw[->,orange!85!black,line width=0.8pt] (10.9,-7.05) -- (5.6,-7.05);
  \node[font=\tiny,orange!70!black,anchor=east] at (10.9,-7.55)
    {rows $e\ge5$: prefixes of $\sigma$};
  \node[font=\tiny,anchor=west,black!70] at (%.2f,-7.05)
    {$\sum_ic_i=1$};""" % (NP + 1.15))
M.append(r"\end{tikzpicture}}")
open('fig_starstair.tex', 'w', encoding='utf-8').write("\n".join(M) + "\n")
print('wrote fig_starstair.tex')
