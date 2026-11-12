# TP53 variant effect — a fifth domain instance of the measurement protocol

**One line.** Competing measures of *TP53 variant deleteriousness* (8 yeast transactivation readouts,
7 computational pathogenicity predictors, and mammalian proliferation DMS screens) are analysed with the
meta-paper's protocol. They mostly behave as the framework predicts: they agree at the extremes, split
into an experimental family and a computational family, and disagree on near-ties, with the comparability
skeleton reproducing across independent assay platforms. One honest exception: a modest genuine middle
residual that kinase lacked. This folds in as **§4.5** of `../paper/main.md`.

Concept κ = functional impact / deleteriousness of a TP53 missense variant. Objects Ω = TP53 missense
variants on the canonical 393-aa transcript (P04637 / NM_000546). This is the *variant-interpretation*
side of the TP53 field (which variant is damaging?), not the drug side (rezatapopt-style Y220C
reactivators) — the drug problem is single-mutant medicinal chemistry and is out of scope here.

## Data (all public; `data/`, assembled by `analysis/build_matrix.py`)
- **8 Kato/Ishioka yeast transactivation promoters** — WAF1, MDM2, BAX, 14-3-3σ, AIP1, GADD45, NOXA,
  P53R2, % of WT activity (NCI *TP53 Database* r21; = Kato *PNAS* 2003). n=2314 missense variants.
- **Seven computational predictors:** AlphaMissense (Cheng 2023), EVE (Frazer 2021), ESM1b (Brandes
  2023), REVEL (Ioannidis 2016), CADD (Kircher 2014, with the v1.6 scores of Rentzsch 2021), PrimateAI
  (Sundaram 2018), BayesDel (Feng 2017). EVE and ESM1b are native, the other four come via dbNSFP v4.1a
  (Liu 2020). Coverage of the core
  panel: ESM1b/REVEL/CADD/BayesDel 2,314, PrimateAI 2,086, EVE 1,942 (its scored region, residues 45–374).
- **DMS proliferation screens** (ProteinGym): Giacomelli 2018 ×3 conditions (WT-Nutlin, Null-Nutlin,
  Null-Etoposide), Kotler 2018. Mammalian cDNA — an experimentally *independent platform* from yeast.
- **ClinVar** germline clinical significance + review-star rating — the independent external anchor.
- **Core panel** = 2,314 variants with all 8 promoters + AlphaMissense (2,308 also Giacomelli, 973 also
  Kotler). ClinVar anchor in-panel: 227 pathogenic / 142 benign / 988 VUS (≥1★: 220 P / 139 B).

## Protocol results (each maps to a step of `../paper/main.md` §2; scripts in `analysis/`)

**Step 1 — orientation vs external anchor (`tp53_anchor.py`).** All 8 transactivation promoters and
AlphaMissense are correctly oriented against ClinVar (path>benign, AUC 0.96–0.99). **All four DMS
proliferation screens read MIS-ORIENTED (AUC 0.016–0.064)** under ProteinGym's "higher = higher
fitness" convention — i.e. higher proliferative fitness = *more* function = *less* deleterious, so the
screens must be flipped before use. A concrete, real instance of the orientation check catching a
backwards measure (the analog of the un-gated kinase candidate), and the reason the reproducibility
script sign-aligns the DMS screens.

**Step 5 — families (`tp53_families.py`).** The 15 core-panel measures split into two families by
measurement basis. The proliferation screens are **held out** of this analysis, so the family means are
over the 8 transactivation readouts and the 7 predictors only.
Within the 8 experimental transactivation readouts, mean *r* = 0.78 (0.66–0.89); within the 7
computational predictors, mean *r* = 0.76 (0.56–0.93); between the two blocks, mean *r* = 0.50
(0.35–0.64). Both blocks are internally coherent and agree less across than within. AlphaMissense tracks
the experimental block best (0.56), CADD worst (0.43). This is the TP53 analog of kinase's family
structure: experimental function and computational prediction induce different orders of the same concept.

**Steps 4 and 7 — current numerical results.** Strict promoter reversals fall from 0.924 to 0.045
across consensus-separation quintiles, and from 0.606 to 0.256 across AlphaMissense-separation quintiles.
The latter score is outside the promoter panel, but a residual does not identify a biological mechanism.
The nine-measure poset is 56.6% incomparable and contains 2,313 equivalence classes among 2,314 variants.
Four dispersed insertion-Gibbs chains pass split-R-hat and batch-MCSE criteria. Median-of-eight agreement
with the uniform average-class-rank aggregate is 0.955, and computational predictors range 0.572–0.729.
The same-start adjacent-swap estimates are superseded. This is similarity to a chosen aggregate, not
an independent biological-validity benchmark.

**Step 3 — cross-platform replication.** Giacomelli directions are fixed to the declared fitness
convention, not fitted to yeast agreement. The yeast and mammalian scalar consensuses agree overall at
0.49. Their per-tertile structure is asymmetric: functional-end agreement is near chance, while the
deleterious end agrees more strongly. Direction agreement by separation does not directly test whether
the platforms reproduce the same incomparable sets. Replication establishes descriptive stability and
does not distinguish concept ambiguity from shared assay resolution or systematic error.

**Companion axis — evaluation overstatement (`tp53_circularity.py`).** Two findings, one an honest null.
(i) A hotspot-codon holdout (ClinGen PM1 set 175/245/248/249/273/282) shows **no inflation** for any
predictor (AUC drops ≤ 0.01): at the unique-variant level ClinVar is not hotspot-dominated (~10% of
pathogenic variants), so the scaffold-memorization analogy does not bite in this slice — reported as a
null. (ii) The real gap holds **across all seven predictors**: each scores 0.90–0.99 against curated
ClinVar but 0.79–0.89 against the experimentally-independent functional truth (Kato ≤20% rule) on the 988
VUS where prediction is actually used — a gap of +0.09 to +0.12 (e.g. AlphaMissense 0.988→0.889, REVEL
0.971→0.866, CADD 0.928→0.813). The endpoints and variant populations differ, so this gap does not isolate benchmark inflation or training leakage.
This is the TP53 face of the ADMET "evaluations overstate performance" sibling result.

## Honest limits / next steps
- **Base-set choice.** The consensus poset keeps the original 9 measures (8 experimental readouts +
  AlphaMissense); the other six predictors are benchmarked against it externally. Making the consensus
  purely experimental is an equally valid alternative — the base-set choice is real, not eliminated
  (cf. paper §5).
- **Coverage varies across predictors.** EVE scores only residues 45–374 (1,942 of the core panel);
  the family and canonical analyses use pairwise-complete correlations, so this is handled but reduces
  power for EVE. dbNSFP v4.1a lacks MetaRNN and the native AlphaMissense/EVE/ESM1b columns of later
  releases; those three are taken from their own sources instead.
- **Independent-platform breadth.** The 8 promoters are one study (Kato/NCI). True cross-platform
  independence comes from the mammalian screens, and in the reproducibility test that means the 3
  Giacomelli conditions alone (Kotler overlaps the core panel too thinly to serve there, 973 variants).
  **Funk 2025** (endogenous-locus
  CRISPR RFS, 9,225 variants — a third, mechanistically distinct platform), **Boettcher 2019**
  (dominant-negative), and the **Fayer 2021** integrated consensus are the obvious next cohorts.
- **Anchor caveat.** AlphaMissense is trained on ClinVar-adjacent signal, so its agreement with ClinVar
  is partly circular; the consensus-vs-ClinVar 0.997 is descriptive of the selected labelled variants and is not an independent validation of the aggregate.
- **Field-formula provenance.** Kato median-of-8 cutoffs (non-functional ≤20%, functional >75%),
  Giacomelli DNE+LOF (z≥0.61 & z≤−0.21), Kotler RFS≥−1.0, and PM1 hotspots are the ClinGen TP53 VCEP
  specifications (Fortuno et al. 2021); used here for orientation of the field's own consensus rule.

## What this adds to the meta-paper

A variant-level example with experimental and predictive agreement clusters, explicit anchor and
selection caveats, and cross-platform scalar comparisons. The corrected sampler supports a useful
median-of-eight approximation, with lower agreement than the earlier insufficiently mixed estimate.
Current numbers are in `analysis/outputs/`; `paper/manuscript.md` gives the controlled interpretation.
