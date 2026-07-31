#!/usr/bin/env bash
# Run every independent verifier for both papers.  Each prints an _OK / _DONE
# sentinel on success and raises on any mismatch, so a non-zero exit means a
# claim failed to reproduce.
set -euo pipefail
cd "$(dirname "$0")"

# Part I: the 1.1397... planar lower bound
CORE=(v1_envelope v2_instances v7_template)
EXTRA=(v3_abstract v5_obstructions v8_networkbox)

# Part II: the 1.28249... record, the deletion-star ceilings, the toolkit
PART2=(p2_v1_k17 p2_v6_pell p2_v5_ladder p2_v3_upper p2_v4_walls
       p2_v12_zerosum p2_v11_vup p2_v10_newmat)

SEARCH=(v4_laminar v6_repair_search p2_v2_star)

run () {
  echo "================================================================"
  echo "== $1"
  echo "================================================================"
  python3 -u "$1.py"
  echo
}

echo "### Part I core: the lower bound"
for s in "${CORE[@]}"; do run "$s"; done

echo "### Part I extended version"
for s in "${EXTRA[@]}"; do run "$s"; done

echo "### Part II: record, ladder, ceilings, upper-bound toolkit"
for s in "${PART2[@]}"; do run "$s"; done

if [[ "${1-}" == "--with-searches" ]]; then
  echo "### search scripts (slow)"
  for s in "${SEARCH[@]}"; do run "$s"; done
else
  echo "(skipping the slow search scripts ${SEARCH[*]}; pass --with-searches to include them)"
fi

echo "ALL VERIFIERS PASSED"
