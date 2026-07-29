# unsplittable-flows

Companion repository for the paper *A Planar Lower Bound 1.1397… > 9/8 for
Cost-Preserving Single-Source Unsplittable Flows* by Sergey Nikolenko in
close collaboration with GPT 5.6 Sol, Claude Fable 5, and Claude Opus 5.

The paper is in [`paper/planar-lower-bound-ssuf.pdf`](paper/planar-lower-bound-ssuf.pdf).
The rest of the repository holds the machine-checkable artifacts behind it:
independent verifiers that rebuild every instance from its raw arc list and
re-derive every numerical claim in exact rational or symbolic arithmetic.

## The result

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
constant is `16/15`. The natural quantitative question is how large a universal
`C` must be, if one exists at all.

This paper shows

```
C  ≥  C*  =  (299 − 41·√41) / 32  =  1.139747070789…  >  9/8,
```

with an integer subfamily of exact critical constant `182359/160000 − 1/n` for
every `n ≥ 68`, and proves that `9/8` is the exact supremum of the
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

## Contents

| path | what it is |
|---|---|
| `paper/` | the paper |
| `verify/` | the independent verifiers (see below) |
| `certificates/` | the exact rational certificate consumed by `v8` |
| `VERIFICATION.md` | verification report: what was checked, how, and which claims in earlier drafts had to be corrected |

## Verification

Every numerical claim is certified in exact rational (`fractions.Fraction`) or
symbolic (SymPy over ℚ) arithmetic. Nothing is sampled and nothing is assumed
about the graph: each program rebuilds the instance from its **raw arc list**,
discovers all simple source–terminal paths by depth-first search rather than
receiving the intended paths as input, and enumerates all integral routings.

| script | certifies |
|---|---|
| `v1_envelope.py` | the two-parameter envelope: blind path discovery, the intervals and row supports, feasibility, the six witness identities, exact Bernstein certificates for all 36 dominance gaps over the whole parameter rectangle, full 16-routing exactness, the factorization of `C′`, the location of `t*`, and `C(t*) = (299−41√41)/32` |
| `v2_instances.py` | raw-arc verification of the seven-vertex `16/15` gadget, the `9/8` triangle family (`n = 2, 6, 68, 97, 1000`), the `8867/7800` four-interval point, and the integer subfamily (`n = 67, 68, 100, 1000, 10⁶`), plus planarity, the `K₄` minor, and non-divisibility of the demands |
| `v7_template.py` | the triangle-template proposition symbolically from the raw arcs: the two pair obstructions, the equalizer `b = 1−r`, the completing-square identity, the `η > 1` branch, and the integer members |
| `v3_abstract.py` | complement-closed one-sided discrepancy, the Hadamard bound, the two-scale (chain-firewall) variant, the stable-set rank LPs for odd holes, antiholes, complete graphs and the Chvátal graph, the `k = 7` abstract lower-box failure, and the lower-box deficit of the calibration family by matching primal and dual exact LPs |
| `v4_laminar.py` | whether the canonical laminar realization is a two-terminal series-parallel digraph (it is not), and the laminar abstract ceiling over 22 240 random rational instances plus an exhaustive exact sweep for `k = 3` |
| `v5_obstructions.py`, `v6_repair_search.py` | the ten-arc local-repair instance, and a search for a genuine local-repair obstruction |
| `v8_networkbox.py` | the seven-terminal signed-tree-network lower-box certificate `δ = 31/250`, rebuilt from its rooted tree and re-realized as a raw DAG |

`v1`, `v2` and `v7` certify the main theorem. The others support the additional
structural material discussed in the paper.

### Running them

```bash
python3 -m pip install -r verify/requirements.txt
cd verify && ./run_all.sh
```

Each script prints an `_OK` sentinel on success and raises on any mismatch, so
a non-zero exit means a claim failed to reproduce. `v1` takes a few minutes
(the Bernstein certificates are exact); `v4` and `v6` are search scripts and
are skipped unless you pass `--with-searches`.

Tested with Python 3.10, SymPy 1.14, NumPy 1.26, SciPy 1.13, NetworkX 3.2.

## License

Code: MIT (see `LICENSE`). The paper text and figures: CC BY 4.0.
