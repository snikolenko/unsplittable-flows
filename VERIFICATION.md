# Independent verification report

Covers both papers: Paper I (Zenodo 10.5281/zenodo.21701162, audited
2026-07-29) and Paper II (Zenodo 10.5281/zenodo.21716713, audited 2026-07-31).

Every claim was re-derived from scratch by the programs in `verify/`, written
from the prose statements rather than from the manuscripts' own auxiliary code.
All graph claims rebuild the raw arc list and discover paths by depth-first
search; all arithmetic is `fractions.Fraction` or SymPy over `Q`.

The Paper I audit is immediately below; the Paper II audit is the final
section.

| script | what it certifies |
|---|---|
| `v1_envelope.py` | four-interval envelope: blind path discovery, intervals, row supports, feasibility, the six witness identities, **exact Bernstein certificates for all 36 dominance gaps on the parameter rectangle**, full 16-routing exactness at sample points, `C'` factorization, `t*`, `C(t*)` |
| `v2_instances.py` | raw-arc verification of the public 16/15 gadget, the 9/8 triangle family for n = 2, 6, 68, 97, 1000, the 8867/7800 four-interval point, and the integer subfamily for n = 67, 68, 100, 1000, 10^6; planarity, K4 minor, non-divisibility of the demands |
| `v3_abstract.py` | complement-closed identity and the Hadamard bound by brute force (n ≤ 8), the two-scale variant, the stable-set rank LPs (odd holes, antiholes, K_n, Chvátal), the k = 7 abstract (LB) failure, and the lower-box deficit of the public triangle by matching primal and dual exact LPs |
| `v4_laminar.py` | whether the canonical laminar realization is a two-terminal series-parallel digraph (it is not); the laminar abstract ceiling over 22 240 random rational instances on all 1069 laminar shapes with k ≤ 6, plus an exhaustive exact rational sweep of 7125 instances for k = 3 (largest value 3/4) |
| `v5_obstructions.py` | the ten-arc "monotone local repair" instance |
| `v6_repair_search.py` | search for a genuine local-repair obstruction |
| `v7_template.py` | the triangle-template proposition symbolically from the raw arcs: A, B, the equalizer b = 1−r, the completing-square identity, the η > 1 branch, the integer members |
| `v8_networkbox.py` | the seven-terminal signed-tree-network lower-box certificate δ = 31/250 |

## Confirmed exactly

* **Theorem (main lower bound).** `C* = (299 − 41√41)/32 = 1.139747070789…`
  The maximizer is `t* = (7 − √41)/4`, the unique critical point of
  `C(t) = 1 + 2t − 8t² + 6t³ − t⁴` in the feasible interval `(0, 2−√3)`.
  Verified: `C ≡ 3/8 + (41/8)t (mod 2t² − 7t + 1)`, giving the closed form by
  hand.
* **Lemma (six pair obstructions).**  All six witness identities hold exactly.
  All 36 dominance gaps were computed symbolically and are strictly positive
  for `0 < t < 1`, `0 < ε < (1−t)²` — a strictly larger region than the
  rectangle in the original note.  Bernstein/Handelman certificates confirm
  nonnegativity on the rectangle without sampling.
* **Corollary (integer subfamily).** `C_n = 182359/160000 − 1/n`, exceeding 9/8
  exactly for `n ≥ 68`; `C_67 = 9/8 − 1947/10720000 < 9/8`.  Every instance has
  exactly two simple paths per terminal, 11 of 16 routings cost-good, tight
  late sets `{0,3}` and `{1,3}`.
* **Proposition (triangle template).** Exact supremum 9/8, attained in the
  limit at `b = (3−η)/4, r = (1+η)/4, σ = (1+η)/2`; integer members
  `(3n−1)²/(8n²)`.  All formulas re-derived from the raw arc list.
* **Proposition (public gadget).** `16/15`, four cost-good routings.
* **Splicing theorem, arborescence form, apex-forest, two-arm decomposition,
  network characterization, clean-pair forest.**  Proofs checked line by line;
  the k = 7 triple system is non-laminar with an explicit crossing pair.
* **Abstract discrepancy.** Complement-closed identity and Hadamard bound
  confirmed by brute force; three-permutation constant confirmed against the
  actual statement of [NN, Lemma 2] (`disc_{L+}^k` is a max over *independent*
  prefix choices in the three permutations, hence the sum of the three maxima
  — so `(k+T+2)/3` is right).
* **Stable-set values.** `C_{2ℓ+1} → 1 + 1/(2ℓ+1)`, `~C_n → 2 − 4/n`,
  `K_n → 2 − 2/n`, Chvátal `4/3` — all exact LP optima.
* **Lower-box deficits.** `1/9 − 2ε/3` on the triangle (primal = dual), and
  `31/250` on the seven-terminal signed-network instance (primal six-term
  mixture and dual weights `1/4, 1/2, 1/4` both reproduce it; raw path counts
  `(2,2,2,2,2,2,2)`).

## Corrections made to the earlier draft

1. **Ring-loading constants were wrong.** The draft said the ring-loading
   upper bound is `+1.4 d_max` and attributed the `+1.1 d_max` lower bound to
   Däubel.  Correct: Däubel's upper bound is **`+1.3 d_max`** (SIAM J. Discrete
   Math. 36(2):867–887, 2022), improving Skutella's `+19/14 d_max` (2016) and
   the classical `+1.5 d_max` of Schrijver–Seymour–Winkler; the `+1.1 d_max`
   lower bound is **Skutella's** (2016).

2. **The Liu–Reis citation is real.** The draft removed it as unlocatable.
   Siyue Liu and Victor Reis, *Weighted Chairman Assignment and Flow-Time
   Scheduling*, ITCS 2026, LIPIcs 362, Article 98, arXiv:2511.18546.  Their
   Theorem 1 gives exactly `(1 − 1/(2m−2)) max_j d_j`, and they note it
   confirms the Morell–Skutella conjecture in the common-order case.

3. **Corollary "no amplification" over-claimed.**  The draft deduced
   "faithful realizations never witness `C ≥ 1`" from Almoghrabi–Skutella–
   Warode.  That theorem's hypothesis is a *two-terminal* series-parallel
   digraph, which has a unique sink; the canonical laminar realization has one
   sink per terminal, and the standard super-sink repair creates a `K4` minor
   already for `F = {{1,2}}`.  The paper now proves what is actually available
   — laminar incidence, an out-arborescence, an apex-forest of treewidth ≤ 2,
   hence the planar bound `2 d_max` — and states the sharper `≤ 1` claim as
   a conjecture, with the equal-demand case known and 22 240 unequal-demand
   laminar instances verified (max `23/24`, supremum 1 approached by an
   explicit two-element family).

4. **The "monotone local repair is false" claim does not hold.**  The ten-arc
   DAG recorded in the project notes reproduces exactly (the recorded `q` matches
   to the last fraction once the path order is fixed as longest-first), but its
   conclusion is wrong: rerouting `t_A` or `t_B` from `suvwt` to `suvt`
   restores the upper box *and* keeps the lower box.  Only the two
   source-arc exchanges break it.  A search over 4855 further small
   three-commodity instances found no genuine witness.  The paper now states
   this accurately as a remark and flags the question as open.

5. **Sign error in the hull separation.**  The separating set is
   `conv U_C(x) + R_{≥0}^A` (upward closed, forcing the normal to be
   nonnegative), not `+ R_{≤0}^A`.

6. **The largest certified lower-box deficit is `31/250`, not `1/9`.**  The
   `31/250` signed-network certificate is realizable as an exact-two-path DAG
   by the network characterization, so it is a genuine flow deficit and
   dominates the `1/9` triangle value.

7. **Ambiguities pinned down.**  Path shares in the ten-arc instance are now
   given by naming the paths rather than by an implicit ordering; the
   four-interval feasibility region is stated as `0 < t < 2−√3`,
   `0 < ε < (1−t)²` rather than as an arbitrary rectangle.

## Paper versions

There are two versions of the manuscript.  The submission version --- the one
in `paper/` --- contains the lower bound and the template optimum, plus a
discussion section stating the conjecture and the open problems.  An extended
version, not yet included here, additionally contains the realization barrier,
the signed-tree-network characterization, the abstract discrepancy dichotomy,
and the lower-box material.  The two share a preamble, the bibliography and the
figures.

Only the submission version's content is covered by the "confirmed exactly"
list above at full strength: the barrier and lower-box material is proved but
has had less independent scrutiny, and its status is recorded in the two
sections that follow.

All vectors are set in bold throughout both versions; components carrying an
index (`x_a`, `d_i`, `f_i`, `z_i`) stay upright, and `c(P)` denotes the scalar
total cost of a path.  Symbols that previously collided were renamed: the
triangle template's middle share is now `sigma`, the pairwise-ceiling slack is
`g_i`, the splicing theorem's terminal indices are `i, j, r`, the odd-hole
index is `ell`, and the Newman-Nikolov prefix indices are `i, j, ell`.

## Not verified / stated as open

* The laminar ceiling conjecture (`≤ 1` with unequal demands).  This is the
  laminar case of Goemans' conjecture; only the equal-demand case is known.
* Whether any instance obstructs general monotone local repair.
* The conjecture that the exact supremum of forced constants is 2.

## Paper II (Zenodo 10.5281/zenodo.21716713)

The second paper was audited the same way: every claim re-derived from the
certificates by programs written from the statements rather than from the
search code.

Confirmed exactly:

* the k = 17 record 1282494797984843521/10^18, rebuilt from its intervals
  alone -- 67-arc list, two simple paths per terminal by DFS, all 131072
  routings in exact integers, 109294 cost-good, 15 minimizers, 8885
  strict-sublevel routings all of rank at most 6, and the 18-atom hull mixture
  whose barycenter is exactly f;
* the lower box on that record: 4799 upper-good routings, primal mixture and
  dual separator both giving delta = 6.3627205289e-10 < d_max;
* the codimension-two theorem, by executing its construction verbatim on
  10 500 random exact instances (3 <= k <= 9) and an exhaustive grid of
  164 160 at k = 4, with the bad set nonempty in a substantial fraction of
  cases so that the hard branch is genuinely exercised;
* B_{k,q} and its equality criterion for all 3 <= k <= 8 and 2 <= q <= k by
  exact LP, including the cases where the criterion predicts strict inequality;
* the Knuth radius facts, the deque count 2^(k-1), and the V_tau <= 2
  characterization by exhaustion;
* the three walls (55/52, 17/16, 45/44), the K4/K5/K6 coherence campaign
  (43 satisfiable, 16 unsatisfiable), and the three-order 6/5 domination
  mixture;
* the k = 6 Pell family rebuilt from its defining formulas, with constants
  recomputed by full enumeration; the k = 7 formula and both algebraic limits;
* the four-item zero-sum wall 229515/458752 = 1/2 + 139/458752 and the
  4 - 2*sqrt(3) cell envelope;
* |B| = 4799 and |V-up| = 19566 recomputed directly from the instance;
* the ladder instances k = 5,6,7,8 of Appendix C: the common-point property,
  two simple paths per terminal, the two chains, and the exact critical
  constant of each drawn member.

Corrections made during the audit:

1. **The k = 16 entry of the sandwich table was wrong** in its last three
   digits. The certificate gives 25497289356771303/20000000000000000 =
   1.27486446783856515, not ...565124.
2. **An equation label inside a display** (eq:gamma-recovery) had no equation
   counter, so the reference to it in the proof picked up an unrelated number.
3. **Two mis-attributions to Paper I.** Neither the lower box (LB) nor the hull
   characterization appears in the published Paper I; both were cut when it was
   reduced to its first three sections. Paper II now states and proves both.
4. **Undefined terminology**: support path, support routing, lower-good, cell,
   crossed two-arm broom, dummy node, and signed tree network -- the last of
   these used in a theorem *hypothesis*, so a reader could not check whether
   the theorem applied.
5. **Table 1 conflated two axes.** "Unconditional" was being read as "holds for
   arbitrary graphs". The table now separates the instances covered from what
   the statement rests on, and says explicitly that the class ceilings bound
   the constant restricted to the common-point / two-order class, not the
   universal constant C.
