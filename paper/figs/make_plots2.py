"""Figures for Part II."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif", "font.size": 9,
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})
BLUE, RED, GREEN, GREY, ORANGE = "#2166ac", "#b2182b", "#1b7837", "#666666", "#d95f02"

CERT = '../../certificates/commonpoint_joint_k17_refined_record_exact.json'
cert = json.load(open(CERT))
IV = cert['intervals']
DEM = [x/1e9 for x in cert['demand_numerators_over_1e9']]
SH = [x/1e9 for x in cert['share_numerators_over_1e9']]
K, NP = len(IV), max(b for a, b in IV) + 1
CP = max(a for a, b in IV)                      # common point

# ================================================= Fig: the k=17 interval system
fig, ax = plt.subplots(figsize=(6.5, 3.3))
order = sorted(range(K), key=lambda i: IV[i][0])
for row, i in enumerate(order):
    a, b = IV[i]
    y = -row
    ax.plot([a, b + 1], [y, y], lw=5.2, color=BLUE, alpha=0.28,
            solid_capstyle="butt", zorder=2)
    ax.plot([a, b + 1], [y, y], lw=0.8, color=BLUE, alpha=0.85,
            solid_capstyle="butt", zorder=3)
    ax.text(a - 0.35, y, f"$t_{{{i}}}$", fontsize=6.6, ha="right", va="center",
            color=BLUE)
    ax.text(b + 1.35, y, f"$d={DEM[i]:.3f}$, $f={SH[i]:.3f}$", fontsize=6.0,
            ha="left", va="center", color=GREY)
ax.axvspan(CP, CP + 1, color=RED, alpha=0.16, zorder=1)
ax.text(CP + 0.5, 1.1, "common arc", fontsize=7.2, ha="center", color=RED)
ax.set_xlim(-2.6, NP + 8.5)
ax.set_ylim(-K + 0.3, 1.9)
ax.set_yticks([])
ax.set_xticks([0, 8, CP, 24, NP])
ax.set_xlabel("spine arc", fontsize=8.5)
ax.spines[["top", "right", "left"]].set_visible(False)
ax.set_title(r"The $k=17$ record: 17 crossing intervals on a 33-arc spine, "
             r"all through one common arc", fontsize=8.8, pad=6)
fig.savefig("fig_k17.pdf")
plt.close(fig)
print("wrote fig_k17.pdf")

# =============================================== Fig: two chains of a common point
fig, ax = plt.subplots(figsize=(6.4, 2.5))
demo = [(0, 6), (1, 8), (2, 5), (3, 9), (4, 7)]
n = 10
for row, (a, b) in enumerate(demo):
    y = -row
    ax.plot([a, b + 1], [y, y], lw=4.6, color=BLUE, alpha=0.26,
            solid_capstyle="butt", zorder=2)
    ax.text(a - 0.3, y, f"$t_{{{row}}}$", fontsize=7.5, ha="right", va="center",
            color=BLUE)
ax.axvspan(4, 5, color=RED, alpha=0.16, zorder=1)
for e in range(n + 1):
    ax.plot([e, e], [-5.4, 0.7], lw=0.5, color=GREY, alpha=0.35, zorder=0)
ax.annotate("", xy=(4.2, 1.15), xytext=(0.2, 1.15),
            arrowprops=dict(arrowstyle="->", lw=1.0, color=GREEN))
ax.text(2.2, 1.32, r"left rows: $R_e=\{i:\,\ell_i\leq e\}$ — a chain, "
        r"the prefixes of $\pi$", fontsize=7.4, ha="center", color=GREEN)
ax.annotate("", xy=(4.8, -6.2), xytext=(10.2, -6.2),
            arrowprops=dict(arrowstyle="->", lw=1.0, color=ORANGE))
ax.text(7.5, -6.6, r"right rows: $R_e=\{i:\,r_i\geq e\}$ — a chain, "
        r"the prefixes of $\sigma$", fontsize=7.4, ha="center", color=ORANGE,
        va="top")
ax.set_xlim(-1.4, 11.2); ax.set_ylim(-7.6, 1.9)
ax.axis("off")
fig.savefig("fig_twochains.pdf")
plt.close(fig)
print("wrote fig_twochains.pdf")

# ==================================================== Fig: the sawtooth construction
fig, ax = plt.subplots(figsize=(6.3, 2.5))
k, q = 12, 3
f = (k - q)/k
n = k - q
# segments: edge ceil(n/2q), internal ceil(n/q)
import math
e_len = math.ceil(n/(2*q))
i_len = math.ceil(n/q)
segs = [e_len] + [i_len]*(q - 2) + [e_len]
extra = (n + 1) - sum(segs)
segs[0] += extra
xs, ys = [0], [0]
cur = 0.0
pos = 0
for si, L in enumerate(segs):
    T = f/2 if si in (0, len(segs)-1) else f
    step = T/L if L else 0
    for _ in range(L):
        pos += 1
        cur += step
        xs.append(pos); ys.append(cur)
    if si < len(segs) - 1:
        pos += 1
        cur -= f
        xs.append(pos); ys.append(cur)
ax.step(xs, ys, where="post", color=BLUE, lw=1.5, zorder=3)
ax.plot(xs, ys, "o", ms=2.6, color=BLUE, zorder=4)
ax.axhline(f/2, color=GREEN, ls="--", lw=0.9)
ax.axhline(-f/2, color=GREEN, ls="--", lw=0.9)
ax.text(0.15, f/2 + 0.015, r"$+f/2$", fontsize=7.6, color=GREEN)
ax.text(11.85, -f/2 + 0.02, r"$-f/2$", fontsize=7.6, color=GREEN, ha="right")
ax.axhline(0, color="k", lw=0.6)
ax.set_xlabel(r"position in the identity order", fontsize=8.5)
ax.set_ylabel(r"cumulative residual $H_j$", fontsize=8.5)
ax.set_title(rf"The sawtooth deletion star ($k={k}$, $q={q}$): every prefix error "
             rf"is at most $f/2=\frac{{1}}{{2}}-\frac{{q}}{{2k}}$", fontsize=8.6, pad=5)
ax.spines[["top", "right"]].set_visible(False)
fig.savefig("fig_sawtooth.pdf")
plt.close(fig)
print("wrote fig_sawtooth.pdf")

# ========================================================== Fig: the record ladder
fig, ax = plt.subplots(figsize=(6.2, 2.9))
ks = [4, 5, 6, 7, 8, 17]
vals = [(299 - 41*np.sqrt(41))/32, 1.1679802157706246, 11 - 4*np.sqrt(6),
        6.5 - 2*np.sqrt(7), 1.2148, 1.282494797984843521]
labels = [r"$\frac{299-41\sqrt{41}}{32}$", r"$T_*$", r"$11-4\sqrt{6}$",
          r"$\frac{13}{2}-2\sqrt{7}$", r"$1.2148$", r"$1.28249\ldots$"]
offs = [(-6, -17), (0, 9), (-14, 8), (8, -17), (4, 8), (0, 9)]
ax.plot(ks, vals, "o-", color=BLUE, lw=1.3, ms=5, zorder=4,
        markeredgecolor="white", markeredgewidth=0.8)
for x, y, l, o in zip(ks, vals, labels, offs):
    ax.annotate(l, (x, y), textcoords="offset points", xytext=o,
                fontsize=7.4, ha="center", color=BLUE, zorder=5,
                bbox=dict(fc="white", ec="none", pad=0.6, alpha=0.85))
for v, lab, col, ls in [(9/8, r"$9/8$ (Part I)", GREY, ":"),
                        (6/5, r"$6/5$", ORANGE, "--"),
                        (4/3, r"$4/3$", RED, "--"),
                        (3/2, r"$3/2$", GREEN, "-.")]:
    ax.axhline(v, color=col, ls=ls, lw=0.9, zorder=1)
    ax.text(3.45, v, lab, fontsize=7.6, color=col, va="center", ha="left",
            bbox=dict(fc="white", ec="none", pad=0.8))
ax.text(11.5, 1.508, r"unconditional ceiling on every rank-$(k{-}1)$ cell",
        fontsize=7.2, color=GREEN, ha="center")
ax.set_xlim(3.3, 18.0); ax.set_ylim(1.10, 1.56)
ax.set_xlabel(r"number of terminals $k$"); ax.set_ylabel(r"certified constant")
ax.set_xticks(ks)
ax.spines[["top", "right"]].set_visible(False)
fig.savefig("fig_ladder.pdf")
plt.close(fig)
print("wrote fig_ladder.pdf")

# ============================================ Fig: B_{k,q} sharpness at every q
fig, ax = plt.subplots(figsize=(5.8, 2.6))
for q, col in zip((2, 3, 5), (BLUE, ORANGE, GREEN)):
    ms = np.arange(1, 26)
    ax.plot(ms, 0.5 - 1/(2*ms), "o-", ms=2.8, lw=1.1, color=col,
            label=rf"$q={q}$:  $B_{{{q}m,{q}}}=\frac{{1}}{{2}}-\frac{{1}}{{2m}}$")
ax.axhline(0.5, color=RED, ls="--", lw=1.0)
ax.text(25.4, 0.5, r"$\frac{1}{2}$", color=RED, fontsize=8.5, va="center")
ax.set_xlabel(r"$m$  (so $k=qm$)"); ax.set_ylabel(r"residual radius $B_{k,q}$")
ax.set_ylim(0, 0.56); ax.set_xlim(0, 26.5)
ax.legend(frameon=False, fontsize=7.6, loc="lower right")
ax.spines[["top", "right"]].set_visible(False)
ax.set_title(r"The deletion-star half frontier is approached at every fixed "
             r"codimension $q$", fontsize=8.6, pad=5)
fig.savefig("fig_allq.pdf")
plt.close(fig)
print("wrote fig_allq.pdf")
