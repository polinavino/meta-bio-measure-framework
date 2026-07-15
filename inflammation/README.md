# Inflammation — three candidate instances of the measurement protocol

**Status: DONE (v1).** Three sub-topics were built and run on public data. Every script writes its console
output to `analysis/outputs/<name>.txt` (path defined inside the script); those tracked files are the
single source of truth for the numbers below. Rerun `bash analysis/run_all.sh` after any script edit so
the outputs and the numbers here stay in lockstep.

Quality summary (honest): **clinical indices = strong** (best reproducibility in the whole project);
**sepsis signatures = strong disagreement, weak reproducibility** (reported as a limitation);
**inflammaging clocks = thin** (only 2 public clocks; a qualitative check, not a full instance).

## Why inflammation
"Inflammatory burden" is a latent concept quantified by many competing formulas that visibly disagree —
the setting the meta-paper's protocol targets (`../paper/main.md`). It splits into three sub-topics of
differing fit to the protocol's machinery (which needs scalar measures that induce *orders*):
1. **Clinical inflammatory indices** — NLR, PLR, MLR, SII, SIRI, CRP, CAR.
2. **Inflammaging clocks** — SImAge, ipAGE (+ canonical cytokine markers).
3. **Sepsis transcriptomic signatures / endotypes** — SRS, MARS, Hallmark inflammatory/IFN-γ.

## Provenance and completeness of the candidate measures (read this)
**Where the measures come from.** For each topic the candidates are measures in *established or routine
use* in that field, taken from named sources — not invented here:
- *Indices:* the systemic-inflammation indices in routine clinical/epidemiological use (NLR, PLR, MLR,
  SII, SIRI) plus CRP and the CRP-to-albumin ratio; all computed from standard CBC + CRP + albumin.
- *Clocks:* published inflammatory-age clocks (SImAge — Kalyakulina 2023; ipAGE — Kalyakulina 2022),
  plus canonical single-cytokine inflammaging markers (CXCL9 — the dominant iAge driver, Furman 2021;
  IL-6; TNF).
- *Sepsis:* published blood signatures — Sepsis Response Signature (SRS/SRSq; Davenport 2016,
  Cano-Gamez 2022), the MARS endotype bi-signature (Scicluna 2017), and two MSigDB Hallmark sets
  (inflammatory response, IFN-γ response).

**How do we know we have all of them? We do not — and cannot.** There is no closed enumeration of "all
measures of inflammation." New indices and signatures are published continually. This is the framework's
own **base-set problem**, stated as a first-class limitation in `../paper/main.md` §5: the consensus
order, the canonical average-rank aggregate, and the incomparable set are all defined *relative to the
chosen base set of measures*; adding or dropping measures can change them. The choice is real, not
eliminated — it is pushed up one level (from "which measure?" to "which set of measures?").

Three things keep the base set from being arbitrary, and none is a completeness proof:
1. **Selection rule stated up front** — a measure is included if it is (a) in routine use / named in
   guidelines or widely cited, (b) computable from the available public data, and (c) a distinct formula,
   not a trivial rescaling of another. The rule is applied before looking at results.
2. **Robustness is tested, not assumed** — the family and reproducibility analyses show whether the
   conclusions depend on any single measure. (They mostly do not: e.g. CRP and CAR are rank-identical, so
   dropping one changes nothing; the family split survives leaving any one index out.)
3. **The base set is reported explicitly** so anyone can extend it and rerun. The protocol is stable
   under extension; the specific numbers are not claimed to be.

So the honest claim is *"these are the commonly-used measures we could compute on public data, and here
is what the protocol says about that set,"* never *"this is every possible measure."*

## Protocol applied (per topic, mirroring `../tp53/` and `../analysis/`)
families (correlation structure) · consensus poset + average-rank canonical (Bubley–Dyer) · near-tie law
(discordance vs separation) · orientation + anchor validation · reproducibility across cohorts.

---

## Topic 1 — Clinical inflammatory indices (NHANES). STRONG.
**Data.** NHANES 2015-16 and 2017-18 (two independent survey cycles = two cohorts). Per participant:
CBC (neutrophils, lymphocytes, monocytes, platelets), high-sensitivity CRP, albumin, age, sex, and
all-cause mortality (NHANES Linked Mortality File). Indices computed by `analysis/build_indices.py`.
Complete 7-index panel: 2015-16 n=6,197 (5,327 with mortality, 224 deaths); 2017-18 n=5,882 (5,104;
113 deaths). Anchor = mortality + age. Measures (all higher = more inflammation): NLR, PLR, MLR, SII,
SIRI, CRP, CAR (=CRP/albumin).

**Results.**
- *Families:* two families by measurement basis — leukocyte-count ratios (NLR/PLR/MLR/SII/SIRI, within-*r*
  mean 0.61) vs CRP-protein (CRP/CAR). Between-family mean *r* = 0.19 (near-independent). **CRP and CAR
  are rank-identical (r=1.00)** — one is redundant.
- *Orientation:* all 7 indices correctly oriented (mortality AUC 0.55–0.72, highest MLR 0.72, lowest
  PLR 0.55; all age *r* > 0).
- *Anchor validates consensus:* consensus mortality AUC 0.669; death rate rises across consensus tertiles
  (0.021 → 0.032 → 0.074); it discriminates mortality at the **extremes (AUC 0.692) but not the middle
  (0.508 ≈ chance)**. Mean age also rises across tertiles (44.6 → 52.0) — the anchor is crude and
  age-confounded, reported as such.
- *Near-tie law:* discordance falls steeply with consensus separation (0.987 → 0.285); a middle residual
  on the consensus axis (midpos coef +0.35) but **on the independent age axis discordance is flat**
  (0.75 → 0.70, R² 0.0015) — age is too weak a concept-proxy to serve as the separation axis here.
- *Canonical aggregate:* poset 72.8% incomparable; best proxied by SII (ρ=0.983), worst by CRP/CAR
  (ρ≈0.33). The canonical is not the best single mortality predictor (in the same 600-sample computation,
  MLR alone AUC 0.77 vs canonical 0.66; only 25 deaths there, so noisy) — it is the minimal-commitment
  consensus, not an outcome-optimised score. (Full-sample per-index mortality AUCs are in `anchor.txt`.)
- *Reproducibility (2015 vs 2017):* the inter-index correlation structure replicates at Spearman **0.995**
  (mean |Δr| 0.013); near-tie slopes and anchor orientations match. **The strongest reproducibility in
  the project.**

**Limits.** CRP≈CAR redundancy; age confounds the mortality anchor; PLR is a weak, less stable predictor
(mortality AUC 0.551 → 0.506 across cycles); poset/canonical use a 600-sample subsample for tractability.

---

## Topic 3 — Sepsis transcriptomic signatures (GSE65682). STRONG disagreement, WEAK reproducibility.
**Data.** GSE65682 (MARS consortium, whole-blood, Affymetrix U219), 802 samples; 479 with 28-day
mortality (114 deaths) and a MARS discovery(263)/validation(216) split. Five published signatures scored
per sample as a mean-*z* over their genes by `analysis/build_sepsis_scores.py`: Hallmark Inflammatory
(148/200 genes), Hallmark IFN-γ (177/200), SRS7 (7/7), SRSq19 (19/19), MARS8 (8/8). Anchor = mortality.
**Scoring caveat:** mean-*z* sign is not a-priori meaningful for the mixed-direction SRS/MARS classifier
genes, so each signature is oriented to positively predict mortality; mortality is therefore not re-used
as an independent validator — the discovery/validation split plays that role.

**Results.**
- *Families / headline:* the **Hallmark-inflammatory family (Inflammatory–IFN-γ r=0.72) and the SRS
  family (SRS7–SRSq19 r=0.84) are near-independent (cross-family mean r=+0.02)**. Two published sepsis
  "severity/inflammation" signatures barely agree — a sharp instance of the framework's premise.
- *Anchor:* individual signatures barely track mortality (raw AUC 0.40–0.50); the oriented consensus
  reaches AUC 0.601 with death rate rising across tertiles (0.161 → 0.259 → 0.294).
- *Near-tie law:* discordance falls with separation (0.990 → 0.386); overall rate is high (0.79) because
  the two families are near-orthogonal, so even well-separated pairs often disagree.
- *Canonical aggregate:* poset 79.2% incomparable; chain convergence ρ=0.96; best proxied by IFN-γ
  (ρ=0.76), worst by Inflammatory (ρ=0.47).
- *Endotype (clustering) caveat:* the consensus only partially separates the given MARS endotypes
  (mean consensus-percentile Mars1 0.62 → Mars4 0.27). Comparing two label systems (SRS vs MARS) needs
  cluster agreement (ARI), not the average-rank machinery — an adjacent problem, not folded in.
- *Reproducibility (discovery vs validation): WEAK/mixed.* The near-tie *shape* reproduces and 4/5
  orientation signs match, but the fine correlation structure replicates only weakly (Spearman **0.32**,
  vs 0.99 for the indices) and MARS8's mortality orientation flips between cohorts. Reported as a
  limitation: the qualitative regularity is stable, the detailed disagreement structure is only partly so.

**Limits.** mean-*z* is an approximation of the official PCA-based SRSq and the MARS classifier; weak
individual mortality signal; weak cross-cohort reproducibility of the detailed structure.

---

## Topic 2 — Inflammaging clocks (SImAge cohort). THIN — a qualitative check only.
**Data.** SImAge cohort (Kalyakulina), n=343 (300 controls, 43 ESRD), 46 serum cytokines + chronological
age. Only two composite inflammatory-age clocks are publicly reconstructable on shared samples: **SImAge**
(per-sample values provided) and **ipAGE** (recomputed here as an ElasticNet age-regression on the 46
cytokines, fit on controls; 37/46 nonzero coefficients; ρ=0.744 vs age). iAge and IMM-AGE are not
publicly obtainable. Two measures are degenerate for the machinery, so the concept is widened to
"inflammatory burden" and three canonical single-cytokine markers are added (CXCL9, IL-6, TNF). Anchor =
age. `analysis/infl_clocks.py`.

**Results.** The two clocks agree strongly (ρ=0.90); the single cytokines agree less with the clocks
(mean 0.65) and among themselves (0.51). Each tracks age (SImAge 0.79, ipAGE 0.74, CXCL9 0.69, IL-6 0.51,
TNF 0.36); the consensus rises monotonically with age (tertile mean age 36.5 → 52.0 → 71.9). Near-tie law
holds (discordance 0.96 → 0.12 with separation); poset 54.1% incomparable.

**Limits.** Only 2 true public clocks; the panel is small (n=343), heterogeneous (clocks + raw cytokines),
and mixes "inflammatory age" with "inflammatory burden." Treat as a qualitative check, not a strong
instance. iAge/IMM-AGE would need their gated cohorts / pipelines.

---

## What the three add to the meta-paper
- The **near-tie law** recurs in all three (discordance falls steeply with separation).
- A **family split by measurement basis** recurs in indices (count-ratio vs CRP) and sepsis
  (Hallmark vs SRS) — as in kinase and TP53.
- **Reproducibility varies by domain**, honestly: excellent for population indices (0.99), weak for
  sepsis signatures (0.32). This is itself informative — it shows the constitutive/artifactual test
  doing its job, distinguishing a stable disagreement structure from a partly sample-specific one.
- One genuinely useful practical finding: the leukocyte-ratio indices are near-independent of CRP, so a
  single "inflammation" number hides which axis is elevated.

## Reproducing
The small derived/tidy tables are committed, so the analysis scripts run from a clean clone and
`analysis/outputs/*.txt` regenerate directly. The large raw inputs are git-ignored (NHANES XPT/mortality
`.dat`, the 187 MB `sepsis_expr.parquet`, the GSE65682 matrix, the GPL13667 annotation); to rebuild the
derived tables from scratch, re-download those per the provenance list above, then run `build_indices.py`
and `build_sepsis_scores.py`. `bash analysis/run_all.sh` runs the full pipeline (~2–3 min; Bubley–Dyer
samplers dominate). Deterministic (seeded subsamples and chains).
