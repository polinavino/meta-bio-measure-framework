# Inflammation analysis scripts

Three topics (see `../README.md`). Python `/Users/polina/miniforge3/bin/python`
(numpy/pandas/scipy/sklearn/pyarrow). Deterministic (seeded subsamples + Bubley–Dyer chains). Every
script defines its own output path and tees console output to `outputs/<name>.txt` — those tracked files
are the single source of truth for the numbers in `../README.md`. Rerun `bash run_all.sh` after any edit.

| Script | Topic | Output | Key result |
|---|---|---|---|
| `build_indices.py` | indices | (writes `../data/indices_*.csv`) | 2 NHANES cycles; complete panel 6197 / 5882 |
| `infl_families.py` | indices | `families.txt` | count-ratio family (mean r 0.61) vs CRP/CAR; between 0.19; CRP≈CAR (r=1.00) |
| `infl_anchor.py` | indices | `anchor.txt` | all oriented OK; consensus mortality AUC 0.669; extremes 0.692 vs middle 0.508 |
| `infl_pairwise.py` | indices | `pairwise.txt` | near-tie law (0.987→0.285); flat on independent age axis (R² 0.0015) |
| `infl_avg_rank.py` | indices | `avg_rank.txt` | poset 72.8% incomparable; SII best proxy 0.983, CRP/CAR worst ~0.33 |
| `infl_reproducibility.py` | indices | `reproducibility.txt` | corr structure replicates **Spearman 0.995** across cycles |
| `build_sepsis_scores.py` | sepsis | (writes `../data/sepsis_scores.csv`) | 802 samples; 5 signatures scored (mean-z) |
| `sep_families.py` | sepsis | `sepsis_families.txt` | Hallmark family vs SRS family near-independent (cross r +0.02) |
| `sep_anchor.py` | sepsis | `sepsis_anchor.txt` | individual mortality AUC 0.40–0.50; consensus 0.601; MARS endotype caveat |
| `sep_pairwise.py` | sepsis | `sepsis_pairwise.txt` | near-tie law (0.990→0.386); high overall (families orthogonal) |
| `sep_avg_rank.py` | sepsis | `sepsis_avg_rank.txt` | poset 79.2% incomparable; IFNg best proxy 0.76 |
| `sep_reproducibility.py` | sepsis | `sepsis_reproducibility.txt` | **WEAK/mixed**: corr structure 0.32, MARS8 sign flips |
| `infl_clocks.py` | clocks | `clocks.txt` | SImAge↔ipAGE 0.90; thin (2 public clocks); near-tie 0.96→0.12 |

**Data provenance.** NHANES 2015-16 / 2017-18 (CDC, CBC + hs-CRP + albumin + Linked Mortality).
GSE65682 (GEO, MARS sepsis) + GPL13667 annotation + MSigDB Hallmark + SRS/MARS gene lists. SImAge/ipAGE
tables (Kalyakulina GitHub). Full source list + the measure-selection rule and its completeness caveat
are in `../README.md`.
