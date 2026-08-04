#!/usr/bin/env python3
"""
Exact verifier for the final 17-terminal SSUFP certificate.

Verification checks
-------------------
* No critical check uses Python ``assert``; checks remain active under ``-O``.
* Every one of the 2^17 routings is recomputed from the raw graph from scratch.
* The closed interval formulas are evaluated for every routing and must
  agree route-by-route with the raw-graph calculation.
* All verification arithmetic uses Python integers or Fraction; no floating
  point arithmetic is used.
* Writing a report is optional, so verification does not depend on directory
  write permissions.

Claim
-----
C = 160325086265636045340018727271
    / 125000000000000000000000000000
  = 1.282600690125088362720149818168...
"""

from __future__ import annotations

from collections import deque
from fractions import Fraction
import argparse
import json
from math import comb
from pathlib import Path
import sys

Q = 10**15
Q2 = Q * Q
K = 17
N = 1 << K

INTERVALS = [
    (0,16),(1,19),(2,17),(3,21),(4,18),(5,25),
    (6,20),(7,27),(8,26),(9,22),(10,30),(11,23),
    (12,31),(13,29),(14,24),(15,28),(16,32),
]

D = [
    1000000000000000, 858393867241515, 1000000000000000,
    858393867241515, 1000000000000000, 788277893510092,
    1000000000000000, 843862199642941, 903365952966881,
    1000000000000000, 856688474148761, 999999999999506,
    799963325105069, 917435984766947, 1000000000000000,
    999999999999996, 1000000000000000,
]

F = [
    141606132761254, 505813311259968, 424206822884194,
    505813311262898, 141606132758738, 193268091609523,
    635928929373725, 65868937079789, 193268090530434,
    884911940543826, 66214441719928, 315124734351381,
    219659854159410, 194258970913591, 288280150799138,
    882527340338262, 341642807653942,
]

EXPECTED_RAW = 1282600690125088362720149818168
EXPECTED_VALUE = Fraction(
    160325086265636045340018727271,
    125000000000000000000000000000,
)
EXPECTED_GOOD = 109294
EXPECTED_BAD = 21778
EXPECTED_MINIMIZERS = {
    "01100100010100110",
    "10100100010100011",
    "01100100010100011",
    "10001100010100011",
    "00101100010100011",
    "10000110010100011",
    "00100110010100011",
}


class VerificationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def route_string(mask: int) -> str:
    return "".join("1" if (mask >> i) & 1 else "0" for i in range(K))


def add_arc(adj, reverse, arcs, u: str, v: str) -> None:
    edge = (u, v)
    require(edge not in arcs, f"duplicate arc: {edge}")
    arcs.append(edge)
    adj.setdefault(u, []).append(v)
    adj.setdefault(v, [])
    reverse.setdefault(v, []).append(u)
    reverse.setdefault(u, [])


def build_graph():
    adj, reverse, arcs = {}, {}, []

    for j in range(33):
        add_arc(adj, reverse, arcs, f"v{j}", f"v{j+1}")

    for i, (a, b) in enumerate(INTERVALS):
        add_arc(adj, reverse, arcs, f"v{a}", f"t{i}")
        add_arc(adj, reverse, arcs, f"v{b+1}", f"t{i}")

    require(len(arcs) == 67, f"expected 67 arcs, got {len(arcs)}")
    require(len(adj) == 51, f"expected 51 vertices, got {len(adj)}")
    return adj, reverse, arcs


def verify_acyclic(adj, reverse):
    indegree = {v: len(reverse[v]) for v in adj}
    queue = deque(sorted(v for v in adj if indegree[v] == 0))
    order = []

    while queue:
        u = queue.popleft()
        order.append(u)
        for v in adj[u]:
            indegree[v] -= 1
            if indegree[v] == 0:
                queue.append(v)

    require(len(order) == len(adj), "graph is not acyclic")
    return tuple(order)


def all_simple_paths(adj, target: str):
    paths = []

    def dfs(u: str, seen: frozenset[str], edges: tuple[tuple[str, str], ...]):
        if u == target:
            paths.append(edges)
            return
        for v in adj[u]:
            if v not in seen:
                dfs(v, seen | {v}, edges + ((u, v),))

    dfs("v0", frozenset({"v0"}), tuple())
    return paths


def construct_instance():
    require(len(INTERVALS) == K, "wrong interval count")
    require(len(D) == K, "wrong demand count")
    require(len(F) == K, "wrong share count")
    require(all(0 < d <= Q for d in D), "demand outside (0,Q]")
    require(all(0 < f < Q for f in F), "share outside (0,Q)")
    require(sum(F) == 6 * Q + 1, "share sum is not 6Q+1")
    require(max(D) == Q, "d_max is not 1")

    adj, reverse, arcs = build_graph()
    topo = verify_acyclic(adj, reverse)
    arc_index = {edge: j for j, edge in enumerate(arcs)}

    early_paths = []
    late_paths = []

    for i, (a, b) in enumerate(INTERVALS):
        paths = all_simple_paths(adj, f"t{i}")
        require(len(paths) == 2, f"terminal {i} has {len(paths)} paths")

        early_exit = (f"v{a}", f"t{i}")
        late_exit = (f"v{b+1}", f"t{i}")
        early = [path for path in paths if early_exit in path]
        late = [path for path in paths if late_exit in path]

        require(len(early) == 1, f"terminal {i}: bad early path count")
        require(len(late) == 1, f"terminal {i}: bad late path count")
        require(late_exit not in early[0], f"terminal {i}: mixed early path")
        require(early_exit not in late[0], f"terminal {i}: mixed late path")

        early_paths.append(early[0])
        late_paths.append(late[0])

    # Per-unit costs: only early exits carry cost 1/d_i.
    costs = {edge: Fraction(0) for edge in arcs}
    for i, (a, _) in enumerate(INTERVALS):
        costs[(f"v{a}", f"t{i}")] = Fraction(Q, D[i])

    def path_cost(path):
        return sum((costs[edge] for edge in path), Fraction(0))

    for i in range(K):
        require(
            Fraction(D[i], Q) * path_cost(early_paths[i]) == 1,
            f"terminal {i}: whole early route does not cost 1",
        )
        require(
            path_cost(late_paths[i]) == 0,
            f"terminal {i}: late route does not cost 0",
        )

    # Fractional graph loads as numerators over Q^2.
    x = [0] * len(arcs)
    for i in range(K):
        early_amount = D[i] * (Q - F[i])
        late_amount = D[i] * F[i]
        for edge in early_paths[i]:
            x[arc_index[edge]] += early_amount
        for edge in late_paths[i]:
            x[arc_index[edge]] += late_amount

    fractional_cost = Fraction(K, 1) - Fraction(sum(F), Q)
    require(
        fractional_cost == Fraction(11 * Q - 1, Q),
        "incorrect fractional cost",
    )

    return {
        "adj": adj,
        "arcs": tuple(arcs),
        "arc_index": arc_index,
        "topological_order": topo,
        "early_paths": tuple(early_paths),
        "late_paths": tuple(late_paths),
        "x": tuple(x),
        "fractional_cost": fractional_cost,
    }


def graph_overload_from_scratch(mask: int, instance) -> int:
    """Recompute all 67 integral loads directly from the selected paths."""
    y = [0] * len(instance["arcs"])
    paths0 = instance["early_paths"]
    paths1 = instance["late_paths"]
    arc_index = instance["arc_index"]

    for i in range(K):
        selected_path = paths1[i] if ((mask >> i) & 1) else paths0[i]
        amount = D[i] * Q
        for edge in selected_path:
            y[arc_index[edge]] += amount

    return max(
        0,
        max(y[e] - instance["x"][e] for e in range(len(y))),
    )


def interval_overload(mask: int) -> int:
    """Closed-form load computation used as an internal consistency check."""
    phi = 0

    for t in range(33):
        load = 0
        for i, (a, b) in enumerate(INTERVALS):
            if a <= t <= b:
                zq = Q if ((mask >> i) & 1) else 0
                load += D[i] * (zq - F[i])
        if load > phi:
            phi = load

    for i in range(K):
        zq = Q if ((mask >> i) & 1) else 0
        early_difference = D[i] * (F[i] - zq)
        late_difference = D[i] * (zq - F[i])
        if early_difference > phi:
            phi = early_difference
        if late_difference > phi:
            phi = late_difference

    return phi


def verify_all_routings(instance):
    visited = bytearray(N)
    good = 0
    bad = 0
    best = None
    minimizers = []

    for mask in range(N):
        require(visited[mask] == 0, f"routing {mask} visited twice")
        visited[mask] = 1

        graph_phi = graph_overload_from_scratch(mask, instance)
        formula_phi = interval_overload(mask)
        require(
            graph_phi == formula_phi,
            f"graph/formula mismatch at mask {mask}: "
            f"{graph_phi} != {formula_phi}",
        )

        late_count = mask.bit_count()
        cost_preserving = (
            Fraction(K - late_count, 1) <= instance["fractional_cost"]
        )
        integer_cost_test = (K - late_count) * Q <= 11 * Q - 1
        require(
            cost_preserving == integer_cost_test,
            f"cost-test mismatch at mask {mask}",
        )

        if cost_preserving:
            good += 1
            require(late_count >= 7, "cost-preserving routing has rank < 7")
            if best is None or graph_phi < best:
                best = graph_phi
                minimizers = [mask]
            elif graph_phi == best:
                minimizers.append(mask)
        else:
            bad += 1
            require(late_count <= 6, "cost-bad routing has rank > 6")

    require(all(visited), "not all routings were visited")
    require(good + bad == N, "routing partition does not sum to 2^17")
    require(good == EXPECTED_GOOD, f"good count {good} != {EXPECTED_GOOD}")
    require(bad == EXPECTED_BAD, f"bad count {bad} != {EXPECTED_BAD}")
    require(
        good == sum(comb(K, r) for r in range(7, K + 1)),
        "cost-preserving count disagrees with binomial count",
    )
    require(best == EXPECTED_RAW, f"best {best} != {EXPECTED_RAW}")

    minimizer_strings = {route_string(mask) for mask in minimizers}
    require(
        minimizer_strings == EXPECTED_MINIMIZERS,
        "minimizer set does not match the expected exact set",
    )
    require(len(minimizers) == 7, "expected exactly seven minimizers")

    value = Fraction(best, Q2)
    require(value == EXPECTED_VALUE, "reduced exact fraction mismatch")

    # Check the example used in the paper.
    example_string = "10100100010100011"
    example_mask = sum(
        (1 << i) for i, bit in enumerate(example_string) if bit == "1"
    )
    require(
        interval_overload(example_mask) == EXPECTED_RAW,
        "paper example is not a minimizer",
    )

    # Determine the exact critical arc of the example.
    critical_spines = []
    for t in range(33):
        load = sum(
            D[i] * (
                (Q if ((example_mask >> i) & 1) else 0) - F[i]
            )
            for i, (a, b) in enumerate(INTERVALS)
            if a <= t <= b
        )
        if load == EXPECTED_RAW:
            critical_spines.append(t)
    require(
        critical_spines == [22],
        f"paper example has unexpected critical spines: {critical_spines}",
    )

    return {
        "all_routings": N,
        "all_arcs_per_routing": len(instance["arcs"]),
        "cost_preserving": good,
        "cost_bad": bad,
        "raw_numerator_over_Q_squared": str(best),
        "exact_value": str(value),
        "minimizers": sorted(minimizer_strings),
        "example_critical_arc": "v22->v23",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--report",
        type=Path,
        default=None,
        help="optional JSON output path",
    )
    args = parser.parse_args()

    instance = construct_instance()
    result = verify_all_routings(instance)

    if args.report is not None:
        args.report.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    print("SSUFP_EXACT_VERIFICATION_OK")
    for key, value in result.items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()
