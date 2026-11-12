#!/usr/bin/env bash
# Regenerate every inflammation output. Each script writes its own outputs/<name>.txt (path defined
# inside the script) and echoes to console. Rerun this after editing any script so the tracked outputs/
# files stay in lockstep with the code and the numbers quoted in ../README.md.
set -euo pipefail
PY=/Users/polina/miniforge3/bin/python
cd "$(dirname "$0")"

# Topic 1 — clinical inflammatory indices (NHANES)
if [[ "${SKIP_BUILD:-0}" != 1 ]]; then $PY build_indices.py; fi
$PY infl_families.py
$PY infl_anchor.py
$PY infl_pairwise.py
$PY infl_reproducibility.py
$PY infl_avg_rank.py            # checked insertion Gibbs on equivalence classes

# Topic 3 — sepsis transcriptomic signatures (GSE65682)
if [[ "${SKIP_BUILD:-0}" != 1 ]]; then $PY build_sepsis_scores.py; fi
$PY sep_families.py
$PY sep_anchor.py
$PY sep_pairwise.py
$PY sep_avg_rank.py             # checked insertion Gibbs
$PY sep_reproducibility.py

# Topic 2 — inflammaging clocks (thin; SImAge cohort)
$PY infl_clocks.py

echo "All inflammation analyses complete; tracked outputs in $(pwd)/outputs/"
