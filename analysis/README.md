# Control analyses (this work)

Reproduces the controlled empirical findings in `../synthesis/repo-findings.md` §6.2 and
`../paper/main.md` §4–§5. Scripts read data from the sibling repos via absolute paths.

- **Python:** `/Users/polina/miniforge3/bin/python` (numpy/pandas/scipy). **R:** 4.5.2.
- Run from anywhere (paths are absolute).

| Script | What it tests | Key result |
|---|---|---|
| `kinase_concentration.py` | Reliability-gate resolution of the zero-active extreme; interior-concentration vs an independence null (consensus coordinate) | zero-active disagreement 66 vs 23 (gated); observed concentration slope 0.27 inside null band [0.18,0.42] → not structural |
| `kinase_concentration2.py` | Same, with n_active as an INDEPENDENT coordinate (avoids circularity) | slope 0.11 inside null [−0.13,0.14]; disagreement flat across selectivity |
| `kinase_pairwise.py` | Gate-cutoff sensitivity; pairwise discordance vs separation and position | discordance driven by separation (0.76→0.19); position adds ~0 (midpos coef −0.11). **No consequential middle** |
| `avg_rank.py` | Average-rank-over-linear-extensions canonical measure (Bubley–Dyer sampler); benchmark 4 measures | poset 37.2% incomparable; Gini +0.955, entropy +0.950, S-score +0.869, ratio +0.796 vs canonical |
| `candidate_bench.py` | Repo's EXACT softplus-entropy candidate vs the 4 measures; floor-domination test | candidate **anti-correlated** (entropy −0.87, S-score −0.94) — inverted by the softplus floor pedestal (93.6% of entries at floor). **Flag for the kinase paper.** |
| `candidate_fix.py` | Pedestal-subtracted hinge; re-run D1–D4 + panel convergence | un-inverts (entropy +1.00, p*=110, D4 100%) **but D3 fails** (−0.94) |
| `candidate_fix2.py` | T-sweep + scale-free variant to test whether both orientation & D3 achievable | no T gives both; original's D3 pass is a pedestal-domination artifact; scale-free (per-compound min) gets both weakly (+0.51) |
| `candidate_fix3.py` | **Repair:** per-compound-quantile-anchored hinge (baseline-free); test all D1–D4 | q10/q25/median anchor satisfies all four — orientation **+0.996**, D3 +1.00 (shift-invariant), D4 100%, p*=110. A minimal fix for the kinase candidate. |
| `transfer.py` | Cross-domain transfer of the separation→discordance curve (kinase/serotonin/clocks/smoking); sign-aligns measures, overlays decile curves + pairwise RMSE (no fitting) | curve monotone-decreasing in all 4 (Spearman −0.99 to −1.00); transfers across kinase/clocks/smoking (RMSE 0.10–0.14); serotonin outlier (0.40–0.50, tiny 13-receptor panel). Needs `data/smoke*.csv` (exported from methylation repo via R). |
| `repro_borrowed.py` | Confirm source-repo-carried numbers: serotonin S-score vs Gini; clock R²/cosine/cell-type from the clock repo's computed parquets | serotonin −0.682 exact; position-clock R² 0.991–0.995; rate–position signed-log R² 0.93–0.95; cell-type max |r| 0.372; cosine off-diag −0.17 to 0.11 (corrects the cited "0.03–0.11") |
| `fasd_auc.R` | FASD per-signature AUC per cohort (from `methylation-biomarker-agreement/data/fasd_all_scores.rds`) | buccal 0.68–0.79; van der Laan (blood) 0.935 discovery / 0.959 replication — matches cited |
| `external_anchor.py` | Demonstrates two refinements across all 4 domains: **(1) orientation check** (each measure vs external anchor — n_active / chronological age / exposure ordinal) and **(2) anchor-validates-consensus** | (1) all standard measures correctly oriented in every domain; the un-gated candidate control is caught (+0.94 vs n_active). (2) consensus–anchor Spearman kinase +0.94 / serotonin +0.82 / clocks +0.65 / smoking +0.68; monotone across anchor tertiles; agreement higher at extremes (kinase 0.92 vs 0.56). Needs `data/smoke*.csv`. |
| `smoking_control.R` | Per-cohort C3 agreement (reproducibility); threshold robustness; threshold-free rank spread | former lowest in both cohorts; former-smoker divergence is a binary-threshold effect; rank-spread metric is density-confounded (not used as evidence) |

**Candidate benchmark: done** (`candidate_bench.py`) — and it surfaced a likely bug in the kinase
paper's candidate measure (softplus floor-domination inverts it). Flagged to the author; not edited in
the kinase repo. **Blocked:** continuous-exposure concentration in smoking — neither GSE50660 nor
GSE42861 has a continuous exposure variable in GEO.
