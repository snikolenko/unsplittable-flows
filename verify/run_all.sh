#!/usr/bin/env bash
# Run every independent verifier.  Each prints an _OK sentinel on success and
# raises on any mismatch, so a non-zero exit means a claim failed to reproduce.
set -euo pipefail
cd "$(dirname "$0")"

CORE=(v1_envelope v2_instances v7_template)
EXTRA=(v3_abstract v5_obstructions v8_networkbox)
SEARCH=(v4_laminar v6_repair_search)

run () {
  echo "================================================================"
  echo "== $1"
  echo "================================================================"
  python3 -u "$1.py"
  echo
}

echo "### core: the lower bound"
for s in "${CORE[@]}"; do run "$s"; done

echo "### extended version"
for s in "${EXTRA[@]}"; do run "$s"; done

if [[ "${1-}" == "--with-searches" ]]; then
  echo "### search scripts (slow)"
  for s in "${SEARCH[@]}"; do run "$s"; done
else
  echo "(skipping the slow search scripts ${SEARCH[*]}; pass --with-searches to include them)"
fi

echo "ALL VERIFIERS PASSED"
