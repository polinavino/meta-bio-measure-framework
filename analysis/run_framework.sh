#!/usr/bin/env bash
# Full local reproduction, including post-remediation audit. No downloads.
set -euo pipefail
cd "$(dirname "$0")/.."
PY=/Users/polina/miniforge3/bin/python
export MPLCONFIGDIR="${TMPDIR:-/tmp}/meta-bio-framework-matplotlib"

# Build dependencies before cross-domain controls read them.
$PY tp53/analysis/build_matrix.py
$PY inflammation/analysis/build_indices.py
$PY inflammation/analysis/build_sepsis_scores.py
SKIP_BUILD=1 bash tp53/analysis/run_all.sh
SKIP_BUILD=1 bash inflammation/analysis/run_all.sh
$PY cvdrisk/analysis/build_cvd.py
$PY cvdrisk/analysis/scores.py
$PY cvdrisk/analysis/cvd_protocol.py
bash analysis/run_all.sh
$PY paper/make_figures.py > analysis/outputs/figures.txt
$PY -m unittest discover -s analysis -p 'test_*.py'
$PY paper/check_claims.py > analysis/outputs/claim_checks.txt
$PY analysis/provenance.py
$PY analysis/reaudit.py > analysis/outputs/reaudit.txt
echo "Full framework workflow and re-audit complete."
