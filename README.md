# unsplittable-flows

Companion repository for two papers by Sergey Nikolenko on cost-preserving
single-source unsplittable flow, both written in close collaboration with
GPT 5.6 Sol, Claude Fable 5, and Claude Opus 5.

| | paper | DOI |
|---|---|---|
| I | [*A Planar Lower Bound 1.1397… > 9/8 for Cost-Preserving Single-Source Unsplittable Flows*](paper/planar-lower-bound-ssuf.pdf) | [10.5281/zenodo.21701162](https://doi.org/10.5281/zenodo.21701162) |
| II | [*Deletion-Star Ceilings and a 1.2825… Lower Bound for Cost-Preserving Single-Source Unsplittable Flows*](paper/Lower_Bound_1_28249_SSUF_Deletion_Star_Nikolenko.pdf) | [10.5281/zenodo.21716713](https://doi.org/10.5281/zenodo.21716713) |

The rest of the repository holds the machine-checkable artifacts behind them:
exact rational certificates, and independent verifiers that rebuild every
instance from its raw arc list and re-derive every numerical claim.

## The problem

An instance of *single-source unsplittable flow* consists of a digraph, a
source `s`, terminals `t_i` with demands `d_i > 0`, and nonnegative arc costs
`c`. Given a feasible fractional flow `x`, one asks for an unsplittable flow
`y` — one path per terminal, carrying the whole demand — that is both
*congestion-good* and *cost-preserving*:

```
y_a  ≤  x_a + C·d_max   for every arc a,          cᵀy  ≤  cᵀx.
```

Dinitz, Garg and Goemans (1999) proved the congestion half alone with `C = 1`.
Goemans conjectured that cost preservation comes for free at the same constant;
that was disproved in July 2026 by a seven-vertex instance whose critical
constant is `16/15`. The question these papers attack is how large a universal
`C` must be, if one exists at all.

## Paper I — the 1.1397… lower bound

```
C  ≥  C*  =  (299 − 41·√41) / 32  =  1.139747070789…  >  9/8,
```

with an integer subfamily of exact critical constant `182359/160000 − 1/n` for
every `n ≥ 68`, and a proof that `9/8` is the exact supremum of the
three-terminal template underlying the known counterexample — so a genuinely
different conflict system was needed.

The construction is a directed spine `v₀ → v₁ → ⋯ → v₅` with two private exit
arcs per terminal. Every terminal therefore has exactly **two** simple
source–terminal paths, so the analysis cannot be invalidated by unintended path
switching — the failure mode that breaks naive gadget constructions. Switching
terminal `i` from its early to its late exit adds `d_i` on a consecutive
interval of spine arcs; the four terminals induce the *crossing* intervals
`[0,2]`, `[0,4]`, `[1,4]`, `[2,3]`, and the crossing is paid for by the
unavoidable baseline loads of the expensive routes rather than by faithful
incidence.

## Paper II — the 1.28249… record and the first ceilings

Three fronts.

**A new record.** A seventeen-terminal common-point interval instance certifies

```
C  ≥  1282494797984843521 / 10^18  =  1.282494797984843521…
```

verified by enumerating all `2^17` routings from the raw 67-arc list, with an
18-atom exact convex-hull certificate. The ladder of families now reads
`k = 4: 1.13975`, `5: 1.16798`, `6: 11−4√6`, `7: 13/2−2√7`, `8: ≈1.2148`,
`17: 1.28249`.

**Ceilings.** Writing the class as a weighted two-permutation prefix system, the
paper introduces the *deletion-star residual* `B`, shows that `B ≤ β` yields
congestion at most `1 + β`, and proves that on the codimension-two cost face
`Σf_i = k − 2` one always has `B ≤ 1/2` — for arbitrary orders, arbitrary
demands (zero demands allowed) and arbitrary shares. Hence congestion at most
`3/2` on every common-point cell of complement mass two. The constant `1/2` is
asymptotically sharp at *every* fixed codimension `q`, with exact value
`B_{mq,q} = 1/2 − 1/(2m)`, and survives a non-divisible two-scale perturbation.

**Structure.** The naive atomless continuum limit of the class has value
*zero*, so the class supremum is a genuinely microscopic quantity; a
cost-preserving transfer through Knuth's two-way rounding network has exact
radius `B(d) = min_λ (λ + Σ|d_i − λ|)`; the lower-box statement holds for all
`2^{k−1}` deque orders and for a maximal switched-interval class of signed tree
networks; and literal three-order overlays are dominated at `6/5`.

> **A scoping note that matters.** The ceilings above bound the constant
> *restricted to the common-point / two-order class*, not the universal `C`. A
> universal constant must serve every instance, and that class does not exhaust
> them. Since every lower bound in both papers comes from that same class, the
> ceilings say how far this construction programme can go — they are limits on
> the method, not on the problem. Only the planar bound of Traub–Vargas Koch–
> Zenklusen (unconditional) and the `(LB) ⇒ 2` implication (conditional) bound
> `C` itself. Table 1 of Paper II states this explicitly.

## Contents

| path | what it is |
|---|---|
| `paper/` | both papers, and the figure generators under `paper/figs/` |
| `verify/` | the independent verifiers (see below) |
| `certificates/` | exact rational certificates consumed by the verifiers and named in the papers |
| `VERIFICATION.md` | verification report: what was checked, how, and which claims in earlier drafts had to be corrected |

## Verification

Every numerical claim is certified in exact rational (`fractions.Fraction`) or
symbolic (SymPy over ℚ) arithmetic. Nothing is sampled and nothing is assumed
about the graph: each program rebuilds the instance from its **raw arc list**,
discovers all simple source–terminal paths by depth-first search rather than
receiving the intended paths as input, and enumerates all integral routings.

### Paper I

| script | certifies |
|---|---|
| `v1_envelope.py` | the two-parameter envelope: blind path discovery, the intervals and row supports, feasibility, the six witness identities, exact Bernstein certificates for all 36 dominance gaps over the whole parameter rectangle, full 16-routing exactness, the factorization of `C′`, the location of `t*`, and `C(t*) = (299−41√41)/32` |
| `v2_instances.py` | raw-arc verification of the seven-vertex `16/15` gadget, the `9/8` triangle family (`n = 2, 6, 68, 97, 1000`), the `8867/7800` four-interval point, and the integer subfamily (`n = 67, 68, 100, 1000, 10⁶`), plus planarity, the `K₄` minor, and non-divisibility of the demands |
| `v7_template.py` | the triangle-template proposition symbolically from the raw arcs |
| `v3_abstract.py` | complement-closed one-sided discrepancy, the Hadamard bound, the two-scale variant, the stable-set rank LPs, the `k = 7` abstract lower-box failure, and the lower-box deficit of the calibration family by matching primal and dual exact LPs |
| `v4_laminar.py` | whether the canonical laminar realization is a two-terminal series-parallel digraph (it is not), and the laminar abstract ceiling |
| `v5_obstructions.py`, `v6_repair_search.py` | the ten-arc local-repair instance, and a search for a genuine local-repair obstruction |
| `v8_networkbox.py` | the seven-terminal signed-tree-network lower-box certificate `δ = 31/250` |

### Paper II

| script | certifies |
|---|---|
| `p2_v1_k17.py` | the `k = 17` record rebuilt from its intervals alone: 67-arc list, DFS path discovery, all `131072` routings in exact integers, the 15 minimizers, the 8885 strict-sublevel routings and their maximum rank, and the 18-atom hull certificate |
| `p2_v5_ladder.py` | the exact lower-box certificate on the record — primal mixture and dual separator over the complete set of `4799` upper-good routings — and the ladder certificates |
| `p2_v6_pell.py` | the `k = 6` Pell family rebuilt from its defining formulas, with exact constants recomputed by full enumeration; the `k = 7` formula and both algebraic limits |
| `p2_v2_star.py` | the constructive proof of the codimension-two theorem executed verbatim on 10 500 random exact instances and an exhaustive grid of 164 160 at `k = 4`; the `B_{k,q}` formula and equality criterion by exact LP; the two-scale bracket |
| `p2_v3_upper.py` | the Knuth radius `B(d)` (median, range, `= d_max` locus, two-scale formula), the deque count and interval characterization, and the `V_τ ≤ 2` characterization |
| `p2_v4_walls.py` | the three exact walls (`55/52`, `17/16`, `45/44`), the complete-graph coherence campaign, and the three-order `6/5` domination mixture |
| `p2_v12_zerosum.py` | the four-item zero-sum wall and the `4 − 2√3` cell envelope |
| `p2_v11_vup.py` | `|ℬ| = 4799` and `|V↑| = 19566` recomputed directly from the instance |
| `p2_v10_newmat.py` | the record-chain sandwich table and the clone/interleaving frontier against their certificates |

### Running them

```bash
python3 -m pip install -r verify/requirements.txt
cd verify && ./run_all.sh
```

Each script prints an `_OK` / `_DONE` sentinel on success and raises on any
mismatch, so a non-zero exit means a claim failed to reproduce. `v1` and
`p2_v1` take a few minutes (the arithmetic is exact throughout); the search
scripts `v4_laminar`, `v6_repair_search` and `p2_v2_star` are skipped unless
you pass `--with-searches`.

Tested with Python 3.10, SymPy 1.14, NumPy 1.26, SciPy 1.13, NetworkX 3.2.

## Citing

```bibtex
@misc{NikolenkoPlanarLower26,
  author       = {Sergey Nikolenko},
  title        = {A Planar Lower Bound $1.1397\ldots > 9/8$ for
                  Cost-Preserving Single-Source Unsplittable Flows},
  year         = {2026},
  howpublished = {Zenodo preprint, \url{https://doi.org/10.5281/zenodo.21701162}},
  doi          = {10.5281/zenodo.21701162},
  url          = {https://doi.org/10.5281/zenodo.21701162}
}

@misc{NikolenkoDeletionStar26,
  author       = {Sergey Nikolenko},
  title        = {Deletion-Star Ceilings and a $1.2825\ldots$ Lower Bound for
                  Cost-Preserving Single-Source Unsplittable Flows},
  year         = {2026},
  howpublished = {Zenodo preprint, \url{https://doi.org/10.5281/zenodo.21716713}},
  doi          = {10.5281/zenodo.21716713},
  url          = {https://doi.org/10.5281/zenodo.21716713}
}
```

## License

Code: MIT (see `LICENSE`). The paper text and figures: CC BY 4.0.
