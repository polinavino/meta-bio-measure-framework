# Controlled framework analyses

Run `bash analysis/run_framework.sh` from any directory for the full local workflow, including TP53,
inflammation, CVD, R controls, figures, mathematical regression checks, provenance and a final audit.
The source datasets and sibling repositories must already be present. Python needs numpy, pandas, scipy
and sklearn; R needs dplyr; `clang++` compiles the native sampler into the system temporary directory.
`run_all.sh` runs only the core/sibling-domain analyses. Domain runners rebuild their derived matrices.

| Files | Purpose |
|---|---|
| `order_utils.py`, `_poset_sampler.cpp` | Strict-reversal/tie semantics, consensus dominance and uniform extension averaging on equivalence classes. Four dispersed insertion-Gibbs chains; split R-hat <1.01 and maximum batch MCSE <0.0025 of class count are required. These are empirical checks, not a finite-time mixing proof. |
| `test_order_utils.py` | Exact small-poset and native-kernel distribution checks; ties, equivalent objects, row and monotone-rescaling invariance. |
| `avg_rank.py` | Tie-preserving gated kinase parameter-median poset and checked aggregate. |
| `kinase_concentration*.py`, `kinase_pairwise.py` | Rank-boundary and separation controls. Active count is a shared-input proxy; pairwise output includes S-score-withheld results. |
| `kinase_controls.py` | Current exponential candidate, conditional construction checks, published versus log-affinity entropy, floor and eligibility sensitivity; conditional rank recovery, retention and unconditional pair recovery on four datasets. |
| `kinase_platform.py` | CID-resolved compounds and shared curated gene panel; direct incomparable-set overlap and matched-compound bootstrap. |
| `serotonin_controls.py` | Observed-only versus floor-imputed scenarios, ratio implementation, fixed coverage counts and identical receptor panels. |
| `protocol_controls.py` | Base-set omission checks, exact-family checks, population denominators and participant-bootstrap correlation-structure replication (sepsis preprocessing and signs are discovery-frozen). |
| `external_anchor.py` | Orientation with declared dependence: active count is shared input, age is a clock training target, smoking is metadata. Historical inverted softplus is a negative control. |
| `transfer.py` | Descriptive consensus-axis curve distances. Directions are declared in advance; strict indicator-pair reversals are used. No held-out predictive-transfer claim. |
| `repro_borrowed.py`, `smoking_control.R`, `fasd_auc.R` | Confirm sibling results and classification/endpoint controls; R stdout is retained in outputs. |
| `provenance.py`, `reaudit.py` | Source/code hashes, sibling revisions and post-remediation checks. |

`candidate_bench.py` and `candidate_fix*.py` are historical softplus/hinge explorations. Their criterion
numbers and one-value checks are not validation of the current candidate. They retain historical
undefined-score conventions; ties are now handled with average ranks. `framework_audit.py` contains the
initial diagnostic reconstruction and a demonstration of the old sampler's problem. Its stored initial
report is historical; `outputs/reaudit.txt` is the current audit.

The current candidate has zero weights at/below the assay floor and exponential weights above it, with
T = tau* - x0 and an explicit activity eligibility gate. Activity eligibility, boundary locality,
weak-addition monotonicity and equal-strength orientation are distinct from finite-panel recovery.
Published classifier coefficients are not reconstructed for the sepsis gene-set proxies. No source is
silently promoted from a proxy to a faithful published implementation.

Numbers in the manuscript are checked by `paper/check_claims.py` for the registered claims. Replication
is conditional stability, not a test that uniquely identifies concept-level incomparability. No continuous
exposure covariate is present for the smoking concentration analysis; that scope limitation remains.
