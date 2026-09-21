#!/usr/bin/env bash
# Regenerate every controlled-re-analysis output. Each script prints to stdout and is captured into
# outputs/<name>.txt. Rerun this after editing any script so the tracked outputs/ files stay in lockstep
# with the code and with the numbers quoted in ../paper/manuscript.md (sections 4.1 to 4.4 and 5).
#
# These scripts read the four source-domain repositories by absolute path (see the paths at the top of
# each file): ../../selectivity, ../../psychedelic-selectivity, ../../epigenetic-clock-desiderata and
# ../../methylation-biomarker-agreement. Clone the sibling repositories next to this one before running.
# The two R scripts (smoking_control.R, fasd_auc.R) need R and are run separately; they read
# ../../methylation-biomarker-agreement.
set -euo pipefail
PY=/Users/polina/miniforge3/bin/python
cd "$(dirname "$0")"
mkdir -p outputs

for s in kinase_concentration kinase_concentration2 kinase_pairwise avg_rank \
         candidate_bench candidate_fix candidate_fix2 candidate_fix3 \
         external_anchor transfer repro_borrowed; do
  echo "== $s"
  $PY "$s.py" > "outputs/$s.txt"
done

echo "All controlled re-analyses complete; tracked outputs in $(pwd)/outputs/"
