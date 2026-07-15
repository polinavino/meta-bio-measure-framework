# TP53 analysis scripts

Fifth domain instance (see `../README.md`). Python `/Users/polina/miniforge3/bin/python`
(numpy/pandas/scipy/sklearn). Deterministic except the Bubley–Dyer sampler (seeded; two-chain
convergence reported). Run order: `build_matrix.py` first, then any analysis.

**Stored outputs = single source of truth for the paper's numbers.** Every script defines its own output
path (`outputs/<name>.txt`) and tees its console output there on each run, so the tracked `outputs/` files
always reflect the current code. Every number quoted in `../README.md`, `../../paper/main.md` §4.5, and
`../../synthesis/repo-findings.md` §6.7 is copied from these files. **Whenever you edit a script, rerun it
(or `bash run_all.sh`) so its `outputs/` file — and any number derived from it — stays in lockstep.**
`run_all.sh` regenerates the matrix + all outputs in dependency order (~2–3 min; the sampler dominates).

| Script | What it tests (protocol step) | Key result |
|---|---|---|
| `build_matrix.py` | Assemble variant × measure matrix from the 5 public sources | 7,486 missense keys; core panel 2,314 (8 promoters + AlphaMissense); ClinVar 227 P / 142 B / 988 VUS |
| `tp53_families.py` | Families by induced order (step 5) | two families: 8 experimental readouts (within-r mean 0.78) vs 7 computational predictors (mean 0.76); **between-block mean 0.50** |
| `tp53_anchor.py` | Orientation vs anchor + anchor-validates-consensus (steps 1, 6) | promoters+AM oriented OK (AUC .96–.99); **all 4 DMS screens MIS-ORIENTED** (.02–.06); consensus AUC vs ClinVar **0.997**; pathogenic frac by consensus tertile 0/5/98% |
| `tp53_avg_rank.py` | Consensus poset + average-rank canonical (step 7) | **56.6% incomparable**; field median-of-8 best proxy **ρ=0.987**; all 7 predictors ρ 0.46–0.59; canonical ClinVar AUC 0.993 |
| `tp53_pairwise.py` | Near-tie law, controlled (step 4) | discordance 0.93→0.05 by separation; **modest residual middle** (midpos coef +0.33 consensus / **+0.10 independent axis**) — unlike kinase |
| `tp53_reproducibility.py` | Constitutive vs artifactual, cross-**platform** (step 3) | yeast vs mammalian consensus ρ 0.49 overall / **+0.63 extremes / −0.05 middle**; direction-agreement 0.53→0.85 |
| `tp53_circularity.py` | Evaluation overstatement (companion axis) | hotspot holdout **null** (all predictors, drop ≤0.01); **all 7 predictors 0.90–0.99 curated ClinVar vs 0.79–0.89 functional-truth-on-VUS** (gap +0.09–0.12) |

**Data provenance.** ProteinGym v1 (Giacomelli 2018 ×3, Kotler 2018); NCI *TP53 Database* r21 =
Kato *PNAS* 2003 (8-promoter transactivation); AlphaMissense (Cheng 2023, Zenodo 8208688); EVE
(evemodel.org, P53_HUMAN); ESM1b (Brandes 2023, HuggingFace `ntranoslab/esm_variants`, P04637);
REVEL/CADD/PrimateAI/BayesDel (dbNSFP v4.1a, Zenodo 4323592, ENST00000269305 slice); ClinVar
(NCBI variant_summary, GRCh38). Field formulas / hotspot codons = ClinGen TP53 VCEP (Fortuno 2021).
Candidate additions not yet folded in: Funk 2025 (endogenous CRISPR RFS), Boettcher 2019 (DN),
Fayer 2021 (integrated consensus).
