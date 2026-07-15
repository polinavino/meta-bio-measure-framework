#!/usr/bin/env bash
# Regenerate the unified matrix and every analysis output for the TP53 instance.
# Each script writes its own outputs/<name>.txt (path defined inside the script) AND echoes to console;
# this runner just runs them in dependency order. Rerun this whenever any script changes so the tracked
# outputs/ files stay in lockstep with the code and the numbers quoted in ../README.md and paper/main.md.
set -euo pipefail
PY=/Users/polina/miniforge3/bin/python
cd "$(dirname "$0")"
$PY build_matrix.py          # writes ../data/tp53_matrix.csv          + outputs/build_matrix.txt
$PY tp53_families.py         #                                          outputs/families.txt
$PY tp53_anchor.py           #                                          outputs/anchor.txt
$PY tp53_pairwise.py         #                                          outputs/pairwise.txt
$PY tp53_reproducibility.py  #                                          outputs/reproducibility.txt
$PY tp53_circularity.py      #                                          outputs/circularity.txt
$PY tp53_avg_rank.py         # writes ../data/tp53_canonical.csv         + outputs/avg_rank.txt  (slow: Bubley-Dyer ~1-2 min)
echo "All TP53 analyses complete; tracked outputs in $(pwd)/outputs/"
