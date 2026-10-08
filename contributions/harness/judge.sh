#!/usr/bin/env bash
# Judges a submission: contributions/submission/ unless another folder is given.
#
# Usage: judge.sh <out_dir> [submission_dir]
#   Writes <out_dir>/report.md, <out_dir>/name (the manifest's name) and
#   <out_dir>/generated/ (the wrapper built from the manifest), and exits 0 if
#   and only if the overall pass rate meets the threshold.
set -euo pipefail

THRESHOLD=50

here=$(cd "$(dirname "$0")" && pwd)
sub=$(cd "${2:-$here/../submission}" && pwd)
out=${1:?usage: judge.sh <out_dir> [submission_dir]}
mkdir -p "$out"
out=$(cd "$out" && pwd)

fail() { echo "::error::$*"; exit 1; }

[ -f "$sub/manifest.yml" ] || fail "submission is missing manifest.yml"
compgen -G "$sub/*.cpp" > /dev/null || fail "submission has no .cpp file"

# generate.py validates the manifest and writes the wrapper; make reruns it.
rm -rf "$out/generated"
if ! make -B -C "$here" algorithm_test SUBMISSION="$sub" GENERATED="$out/generated"; then
    if [ -f "$out/generated/signature" ]; then
        fail "build failed; contract $(sed -n 's/^contract:[[:space:]]*//p' "$sub/manifest.yml") expects exactly: $(cat "$out/generated/signature")"
    fi
    fail "build failed; see the error above"
fi
name=$(cat "$out/generated/name")
echo "$name" > "$out/name"

# The harness exit code is ignored on purpose: the gate here is the pass rate,
# not whether every case passed. A crash leaves no gate file and is rejected.
"$here/algorithm_test" "$out/report.md" "$out/gate.txt" || true
[ -s "$out/gate.txt" ] || fail "the harness wrote no results (did the algorithm crash?)"

read -r passed judged < "$out/gate.txt"
[ "$judged" -gt 0 ] || fail "no declared case was judged; check models/geometry/lattices in manifest.yml"

# Rounded for display, like the report; the gate itself compares exactly.
pct=$(( (passed * 200 + judged) / (2 * judged) ))
echo "Contributed algorithm '$name' passed $passed of $judged cases (${pct}%)."
[ $(( passed * 100 )) -ge $(( THRESHOLD * judged )) ] \
    || fail "pass rate ${pct}% is below the ${THRESHOLD}% threshold, so the submission is rejected"
echo "Pass rate ${pct}% meets the ${THRESHOLD}% threshold, so the submission is accepted."
