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
- **Seven computational predictors:** AlphaMissense (Cheng 2023), EVE, ESM1b, REVEL, CADD, PrimateAI,
  BayesDel (EVE + ESM1b native; REVEL/CADD/PrimateAI/BayesDel via dbNSFP v4.1a). Coverage of the core
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

**Step 5 — families (`tp53_families.py`).** The measures split into two families by measurement basis.
Within the 8 experimental transactivation readouts, mean *r* = 0.78 (0.66–0.89); within the 7
computational predictors, mean *r* = 0.76 (0.56–0.93); between the two blocks, mean *r* = 0.50
(0.35–0.64). Both blocks are internally coherent and agree less across than within. AlphaMissense tracks
the experimental block best (0.56), CADD worst (0.43). This is the TP53 analog of kinase's family
structure: experimental function and computational prediction induce different orders of the same concept.

**Step 7 — consensus poset + canonical aggregate (`tp53_avg_rank.py`).** The 9-measure consensus poset
is **56.6% incomparable** — the nine measures jointly order fewer pairs than kinase's four did (37%). The
weight-free average-rank-over-linear-extensions canonical measure (Bubley–Dyer sampler, two chains
ρ=0.998) is best proxied by **the field's own proposed formula — the median of the 8 promoters
(ρ=0.987)** — and worst by **AlphaMissense (ρ=0.597)**, quantifying that the predictor commits most
beyond the consensus. The other computational predictors sit in the same range (ρ 0.46–0.59: EVE 0.59,
REVEL/BayesDel 0.54, PrimateAI 0.52, ESM1b 0.50, CADD 0.46) — the whole computational family diverges
from the experimental consensus, consistent with the family split. Notably the field's Kato/ClinGen
median-of-8 rule *already is* an excellent approximation of the principled minimal-commitment aggregate:
the protocol here vindicates existing practice rather than overturning it. Canonical and median both
discriminate ClinVar at AUC ≈ 0.99.

**Step 4 — locate disagreement, controlled (`tp53_pairwise.py`).** Discordance among the 8 promoters
(2.68M pairs; overall 0.48) is **dominated by separation** (near-tie-ness): 0.93 for the closest
separation quintile → 0.05 for the widest. Unlike kinase, a **modest residual middle effect survives**
(mid-axis near-ties disagree more than end near-ties): midpos coefficient +0.33 on the consensus axis,
shrinking to **+0.10 on the independent AlphaMissense axis**. Honest reading: the near-tie law holds and
dominates, but TP53 retains a small *genuine* middle component — the partial-LOF / separation-of-function
variants — that kinase did not. A real, interpretable cross-domain difference; not the retired P3
"consequential middle."

**Step 3 — constitutive vs artifactual, cross-platform (`tp53_reproducibility.py`).** Yeast
transactivation vs mammalian proliferation — different organism, readout, lab. Consensus orders agree
**ρ=0.49 overall, +0.63 at the extremes, −0.05 in the middle**; cross-platform direction-agreement
rises from **0.53 on near-ties (≈chance) to 0.85 on well-separated pairs**. The comparability skeleton is
reproducible across independent platforms (**constitutive**); both platforms lose resolution on the
*same* near-ties (the concept is genuinely silent there), rather than one assay's noise.

**Companion axis — evaluation overstatement (`tp53_circularity.py`).** Two findings, one an honest null.
(i) A hotspot-codon holdout (ClinGen PM1 set 175/245/248/249/273/282) shows **no inflation** for any
predictor (AUC drops ≤ 0.01): at the unique-variant level ClinVar is not hotspot-dominated (~10% of
pathogenic variants), so the scaffold-memorization analogy does not bite in this slice — reported as a
null. (ii) The real gap holds **across all seven predictors**: each scores 0.90–0.99 against curated
ClinVar but 0.79–0.89 against the experimentally-independent functional truth (Kato ≤20% rule) on the 988
VUS where prediction is actually used — a gap of +0.09 to +0.12 (e.g. AlphaMissense 0.988→0.889, REVEL
0.971→0.866, CADD 0.928→0.813). The curated benchmark overstates skill on the decision-relevant variants.
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
  independence comes only from the Giacomelli/Kotler mammalian screens; **Funk 2025** (endogenous-locus
  CRISPR RFS, 9,225 variants — a third, mechanistically distinct platform), **Boettcher 2019**
  (dominant-negative), and the **Fayer 2021** integrated consensus are the obvious next cohorts.
- **Anchor caveat.** AlphaMissense is trained on ClinVar-adjacent signal, so its agreement with ClinVar
  is partly circular; the consensus-vs-ClinVar 0.997 is safe because it is carried by the 8
  experimentally-independent transactivation measures.
- **Field-formula provenance.** Kato median-of-8 cutoffs (non-functional ≤20%, functional >75%),
  Giacomelli DNE+LOF (z≥0.61 & z≤−0.21), Kotler RFS≥−1.0, and PM1 hotspots are the ClinGen TP53 VCEP
  specifications (Fortuno et al. 2021); used here for orientation of the field's own consensus rule.

## What this adds to the meta-paper
A **fifth, independent, clinically-hot domain** where the protocol runs end to end on public data, with
an unusually strong anchor (ClinVar + experimental function) and a genuine cross-*platform* (not just
cross-cohort) reproducibility test. It reproduces the near-tie law and a two-family structure,
vindicates a field-standard aggregation formula as an approximation of the principled canonical measure,
and — via the VUS-vs-functional-truth gap — reconnects the definitions axis to the ADMET evaluations
axis. It also contributes one honest null (no hotspot inflation) and one honest cross-domain *difference*
(a small genuine middle residual absent in kinase), both of which strengthen rather than dilute the
protocol's credibility.
