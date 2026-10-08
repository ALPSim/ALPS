#!/usr/bin/env bash
# Judges whatever is in contributions/submission/.
#
# Usage: judge.sh <out_dir>
#   Writes <out_dir>/report.md and <out_dir>/name (the manifest's name), and
#   exits 0 if and only if the overall pass rate meets the threshold.
set -euo pipefail

THRESHOLD=50

here=$(cd "$(dirname "$0")" && pwd)
sub=$(cd "$here/../submission" && pwd)
out=${1:?usage: judge.sh <out_dir>}
mkdir -p "$out"

fail() { echo "::error::$*"; exit 1; }

[ -f "$sub/algorithm.hpp" ] || fail "submission is missing algorithm.hpp"
[ -f "$sub/manifest.yml" ]  || fail "submission is missing manifest.yml"
compgen -G "$sub/*.cpp" > /dev/null || fail "submission has no .cpp file"

name=$(sed -n 's/^name:[[:space:]]*//p' "$sub/manifest.yml" | head -n 1 | tr -d "\"' \r")
[[ "$name" =~ ^[A-Za-z0-9_-]+$ ]] \
    || fail "manifest.yml needs a 'name:' made of letters, digits, '_' or '-' (got '$name')"
echo "$name" > "$out/name"

make -B -C "$here" algorithm_test SUBMISSION="$sub"

# The harness exit code is ignored on purpose: the gate here is the pass rate,
# not whether every case passed. A crash leaves no gate file and is rejected.
"$here/algorithm_test" "$out/report.md" "$out/gate.txt" || true
[ -s "$out/gate.txt" ] || fail "the harness wrote no results (did the algorithm crash?)"

read -r passed judged < "$out/gate.txt"
[ "$judged" -gt 0 ] || fail "no case was judged"

pct=$(( passed * 100 / judged ))
echo "Contributed algorithm '$name' passed $passed of $judged cases (${pct}%)."
[ "$pct" -ge "$THRESHOLD" ] \
    || fail "pass rate ${pct}% is below the ${THRESHOLD}% threshold, so the submission is rejected"
echo "Pass rate ${pct}% meets the ${THRESHOLD}% threshold, so the submission is accepted."
