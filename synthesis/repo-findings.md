# Reconciled ground truth from the five repos

**Purpose.** The `meta-paper-plan.md` was written from memory. This file is the *verified* version,
extracted directly from the repos (not the ChemRxiv preprint). It exists so no future conversation
has to re-derive the numbers or the exact axiom text. Where the plan and the repos disagree, the
repo wins and the discrepancy is flagged. Cross-domain **tensions** (things that complicate the
clean thesis) are called out explicitly — they are the parts a reviewer will attack.

Local paths (all siblings of this repo under `/Users/polina/Documents/BioInfStuff/`):
`selectivity` (= kinase-selectivity-definitions), `epigenetic-clock-desiderata`,
`psychedelic-selectivity`, `methylation-biomarker-agreement`, `admet-qsar-evaluation`.

> **⚠ READ §6 FIRST (added later in the project).** The paper was re-scoped from a "new-mathematics"
> claim to a **methods/Perspective** paper (a protocol + justification). Three adversarial referees and
> a set of **controlled re-analyses** substantially changed the empirical story — most importantly, the
> "concentration in the consequential middle" (P3) claim **did not survive controls** and has been
> retired. §6 records the referee findings, the controlled numbers, the prior-art map, and the
> decisions. The live paper docs are `../paper/main.md` (manuscript) and `../paper/formal-spine.md`
> (Appendix A). `../paper/framework-core.md` is **superseded** (old new-math architecture).

---

## 1. Completeness / maturity of each instance (what actually exists as text)

| Instance | Canonical writeup | Form | Maturity |
|---|---|---|---|
| Kinase | `selectivity/paper/sections/*.tex` (+ `main.pdf`) | Full LaTeX prose draft, all sections | **Most mature.** Under revision for *Molecular Informatics* (resubmission); ChemRxiv preprint live; has reviewer responses. |
| Clocks | `epigenetic-clock-desiderata/README.md` (250 lines) | README **is** the draft; all `paper/*.tex` are **0-byte stubs** | Results + prose done; **manuscript not written**. 19 figures exist. Also a separate `intervention-framework/` sub-project with real LaTeX + its own *different* desiderata — do not conflate. |
| Serotonin | `psychedelic-selectivity/README.md` (381 lines) | README is the writeup | Complete & self-contained. Does **not** restate or test D1–D4; only *reuses* the four measures. |
| Methylation | `methylation-biomarker-agreement/docs/writeup.md` (334 lines) | Structured manuscript (Abstract→Refs) | Complete draft; terse/telegraphic register (dropped articles — that's the actual style, not truncation). |
| ADMET | `admet-qsar-evaluation/README.md` (205 lines) | README + per-analysis `.txt` reports | Complete & internally consistent. |

Only the kinase instance has a real formal (LaTeX) manuscript. Only the kinase instance states its
desiderata as formal predicates. The meta-paper cannot lean on "four polished papers"; it leans on
**one mature formal instance (kinase) + three complete-but-informal instances + one sibling.**

---

## 2. The four domain axiom-sets, verbatim, and their reconciliation

### 2.1 Kinase — D1–D4 (from `selectivity/paper/sections/desiderata.tex`)

Setup: profile `x = (x_1,…,x_N) ∈ ℝ^N`, `x_i = pK_{d,i}`; measure `f : ℝ^N → ℝ`, higher = more
selective; ranking = order of `f`-scores.

- **D1 Reliability threshold.** ∃ τ*>0 s.t. `f(x)` is undefined / flagged unreliable when `max_i x_i ≤ τ*`.
- **D2 Bounded gap sensitivity.** Small perturbation to the *identity* of the reference off-target ⇒ bounded change in `f`. (Ratio violates when top-1 ≈ top-2.)
- **D3 Distributional consistency.** Ranking stable under small perturbation of the activity baseline β. (Entropy, Gini fail — they're β-parameterized.)
- **D4 Monotonicity under weak off-target addition.** Appending `x_{N+1} ≤ τ*` gives `f(x') ≤ f(x)`. (S-score, Gini violate — you can look more selective by profiling against more inactive kinases.)

Satisfaction (verbatim table; ✓/✕/∼):

| | S-score | Entropy | Gini | Ratio | Candidate |
|---|---|---|---|---|---|
| D1 | ✕ | ✕ | ✕ | ✕ | ✓ |
| D2 | ∼ | ✓ | ✓ | ✕ | ✓ |
| D3 | ∼ | ✕ | ✕ | ✓ | ✓ |
| D4 | ✕ | ✓ | ✕ | ✓ | ✓ |

"No existing definition satisfies all four simultaneously; the candidate does." *(This is the kinase
paper's verbatim claim. Our controlled re-analysis in §6.2 finds the candidate as implemented does
**not** hold — it is inverted, and the fix trades orientation against D3. Read §2.1 as attribution, §6.2
as the finding.)* Section was renamed
"Required Properties for a Well-Formed Selectivity Measure" per reviewer request.

### 2.2 Clocks — D0–D4 (from `epigenetic-clock-desiderata/README.md`)

- **D0 Type declaration.** Every clock declares position (𝒯_τ), deviation (𝒯_δ), or rate (𝒯_τ̇). Cross-type comparisons are not well-formed.
- **D1 Monotonicity.** Expected output increases with chronological age across the population.
- **D2 Directional stability.** Sensitive to motion *along* the canonical trajectory, insensitive to motion *off* it (noise, cell-type shifts).
- **D3 Rank consistency.** Two clocks of the *same declared type* agree on the biological-age ordering of individuals.
- **D4 Intervention type consistency.** Clocks respond to interventions consistently with declared type (position→accumulated change; rate→change in current rate).

### 2.3 Methylation — C1–C4 (from `methylation-biomarker-agreement/docs/writeup.md`)

Signatures `s_1..s_k`, each maps methylation profile M → scalar score.

- **C1 Sensitivity threshold.** Below some exposure level, classification isn't reliable; score distributions of exposed/unexposed must separate at the threshold. (Author maps **C1 = kinase D1**.)
- **C2 Score stability.** Adding a CpG with effect size below δ should change the sample ranking by less than a bound ρ. (Author's README maps **C2 = "D3"**; but see tension 3.1 — the *content* is closest to kinase D4 ∪ D3.)
- **C3 Cross-cohort consistency.** Signatures rank samples consistently across independent cohorts (operationalized as inter-signature agreement rate within each exposure group, not matched-sample correlation).
- **C4 Monotonicity.** Increasing cumulative exposure monotonically increases (or decreases) the score. (Author maps **C4 = "D4"**.)

### 2.4 Serotonin — inherits kinase D1–D4

Does not restate or test them; reuses the four measures {S-score, entropy, Gini, ratio} + a bias
factor (Kenakin operational-model transduction ratio) for biased agonism.

### 2.5 Proposed unified schema (G-axioms) — folded into the paper (§2 protocol step 5 + "Sharpening the desiderata")

The domain lists are reorderings/projections of a small domain-independent set. Proposed mapping
(this is a *synthesis*, to be defended in the paper, not something any single repo states):

| General axiom | Kinase | Clocks | Methylation | Meaning (measurement-theoretic) |
|---|---|---|---|---|
| **G0 Type/scale declaration** | (implicit: "ratio answers a different question") | **D0** | (implicit) | Declare the empirical relational structure + scale type; cross-type comparison ill-formed. |
| **G1 Reliability / domain gate** | **D1** | — | **C1** | The homomorphism is defined only where the relation is observable. *Can be load-bearing for correctness, not just a flag* — kinase: without the gate the measure inverts (§6.2). |
| **G2 Monotonicity in the concept** | **D4** (directional: junk↛selective) | **D1** (age↑→out↑) | **C4** (exposure↑→score↑) | Order-preservation (the homomorphism condition), *verified against an external anchor, not internal consistency* — a measure can be stable + monotone yet oriented backwards (kinase inversion, §6.2). |
| **G3 Ordinal stability under nuisance** | **D2, D3** | **D2** | **C2** | Induced order invariant to *analyst-chosen* nuisance parameters; *apparatus-fixed* ones (detection floor, assay type) are declared context, not varied (kinase floor-vs-baseline, §6.2). |
| **G4 Cross-instance reproducibility** | (panel-size stability) | **D3** | **C3** | The order is a property of the concept, not the sample; replicates across cohorts/panels. |
| **G5 Intervention consistency** | — | **D4** | — | Order preserved under the relevant operations/interventions. |

Note the numbering across repos is **not** harmonized by the author (kinase D1=reliability;
clocks D1=monotonicity; methylation "C2=D3" but content ≈ D4). Harmonizing them into G0–G5 and
showing each list is a projection **is** a genuine deliverable of §3.

---

## 3. Cross-domain TENSIONS (the honest cracks — address, don't hide)

**3.1 The numbering doesn't line up across repos.** Kinase D1 = reliability, clock D1 = monotonicity,
and the methylation README maps C2→"D3" though its content (stability when adding a weak feature) is
closest to kinase **D4** (monotonicity under weak addition) blended with **D3** (distributional
stability). The G-schema above resolves this, but the paper must *show* the resolution, not assert it.

**3.2 "Families" membership is NOT stable across domains — the biggest crack.**
- Kinase (Davis/Klaeger): distribution family {S-score, entropy, Gini} within-r 0.74–0.99 (Davis reaches 0.999); **ratio is the outlier** (cross-family r 0.14–0.62). "Ratio is the consistent outlier."
- Serotonin ChEMBL (n=297): {S-score, Gini, ratio} cluster; **entropy is the outlier** (entropy vs all ≈ 0; S-score–Gini = −0.682). Shulgin/PDSP (n=36): entropy again the outlier.
- So *which* measure is categorically different flips between domains. Contributing factors: (a) **direction-convention inconsistency** — the serotonin scripts compute S-score as fraction *above* vs *below* threshold in different sub-analyses, and measures mix "higher=selective" vs "lower=selective", so raw correlation *signs* are not directly comparable to kinase; (b) **panel size** — 13 serotonin receptors vs 300–433 kinases; entropy's behavior is panel-size-dependent (it only stabilizes above ~110 targets in kinase data), so on 13 receptors it may genuinely measure something different.
- **Safe claim for the paper:** the *phenomenon* "measures partition into a low-disagreement cluster plus at least one categorically-different outlier, by the order they induce" recurs; the *identity* of the outlier is scale- and panel-dependent. Do **not** claim entropy (or ratio) is universally the odd one out. This is more honest and still supports P2.

**3.3 The methylation smoking arm does NOT split into families.** Its four signatures agree broadly
(r 0.81–0.97), unlike kinase's two-family split. Its structured disagreement shows up as the
**intermediate-zone classification divergence** (former smokers 71% vs never 96% / current 82%), not
as families. So P2 (families) and P3 (intermediate-zone concentration) are *separable* phenomena;
not every instance exhibits both. The paper should present P2 and P3 as distinct predictions of the
framework, each instantiated where it occurs, rather than implying all four instances show all five Ps.

**3.4 ADMET is not in the same measurement-theory language at all.** The repo frames itself purely as
a CV-methodology / scaffold-generalization question ("evaluations overstate performance"); it has no
desiderata, no order relation, no families. The "second face of the same problem" reading is an
*interpretive* mapping we would impose. → leans toward **companion citation, not deep fold-in** (§6
decision). See `../paper/main.md` §5 (companion axis) and §6 for the argument.

**3.5 Non-independence caveat (methylation FASD).** The strong blood episignature (van der Laan 2025)
shares lab/tissue/likely samples with the eval cohorts → treat as **upper bound**, not clean
out-of-sample. Does not affect the two robust findings (independent buccal signatures transfer poorly
to blood; between-signature agreement results are invariant to which scores best). Keep this caveat
visible if the FASD arm is cited.

---

## 4. Key numbers per instance (for the evidence subsections)

**Kinase.** Datasets: Davis 68×433 (competition binding), Klaeger 222×343 (chemoproteomics),
Anastassiadis 178×300 (functional %inhibition), Metz 704×172 (pK_i). Two-family clustering
replicates across all four (within-distribution r 0.74–0.99; ratio-vs-distribution 0.14–0.62).
Rank instability concentrated in: zero-active compounds (σ≈74 vs 32), near-tied top targets (top1–top2
gap predicts ratio instability r=−0.34, p<0.001), broad-profile+dominant-target (entropy–ratio
disagreement). Panel size to reach Spearman 0.90 vs full panel: **entropy ~110, S-score ~140,
Gini ~290, ratio only at full 343**. Candidate: entropy + D1 gate + softplus hinge (`w_i =
T·log(1+e^{(x_i−x_0)/T})`, floor x_0, width T≈1); satisfies D1–D4; converges at p*≈170
("half a large panel"); **explicitly non-unique** (T is a design choice; T→0 recovers standard
entropy + its D3 failure).

**Clocks.** GSE40279 (Hannum, n=656) + GSE87571 (n=729) = 1,385 blood 450k. Canonical trajectory =
principal curve (Hastie–Stuetzle) under age-informativeness norm `‖·‖²_* = Σ w_i(·)², w_i = R²_i·σ²_i`,
top 200 CpGs; decompose profile into position τ (arc-length projection) + off-manifold residual r.
Types: position 𝒯_τ (Horvath, Hannum), deviation 𝒯_δ (PhenoAge, GrimAge), rate 𝒯_τ̇ (DunedinPACE).
Results: position clocks **linearly** related R²=0.991–0.995 ("same measurement, different scale");
rate↔position related by **signed log** R²=0.93–0.95 vs 0.75–0.79 linear (**Weber–Fechner**, coeff
a≈0.0135); **all clocks fail D3** cell-type confounding (max |r|: Hannum 0.372); **GrimAge D0
violation** (two-stage composite), pairwise coefficient-vector cosine **near-orthogonal** (reproduced
from the clock repo's `cosine_similarity.parquet`: off-diagonal −0.17 to 0.11, incl. Horvath–DunedinPACE
= 0.000; the earlier "0.03–0.11" figure omitted the −0.171 PhenoAge–DunedinPACE value — all |cos| ≤ 0.17,
so "near-orthogonal" stands; GrimAge excluded from the matrix).
Uses Kendall's κ not Pearson (ordinal-level agreement, deliberate). Closest precursor: Klemera &
Doubal (2006) — same complaint, pre-methylation, linear clinical biomarkers only.

**Serotonin.** ChEMBL 13,584 compounds × 13 human 5-HT receptors; 297 with the required overlap.
Correlation of the four definitions: **entropy ≈ 0 with all; S-score–Gini = −0.682**. Psilocin:
non-selective (highest affinity 5-HT2B pK_i 8.34 & 5-HT1D, not 5-HT2A 6.72) — "inconsistent with its
common description as a '5-HT2A agonist.'" Biased-agonism QSAR: **LOO R²≈0.71–0.78 collapses to LOSO
R²<0** (drop ~0.9) — pure scaffold confounding; this is the case that *motivated* the ADMET repo.
Bias factor = log10[(E_max,βarr/EC50_βarr)/(E_max,Gq/EC50_Gq)].

**Methylation.** Smoking (within-tissue): AHRR(1) / EpiSmoke(4) / Joehanes(9) / EpiTob(12 CpGs),
cohorts GSE50660 n=464 + GSE42861 n=689 (blood 450k). C3 agreement: never **96%** (358/372),
current **82%** (183/222), **former 71%** (350/491) — intermediate-zone divergence, consistent across
cohorts. All pass C4 monotonicity; all fail C1 (formers score intermediate). FASD (cross-tissue):
Portales(657)/Lussier(161,183, buccal) / van der Laan(204, **blood**, 2025). Buccal→blood **AUC
0.68–0.79** (don't track severity) vs blood-native **AUC 0.93–0.96** (but upper bound, tension 3.5).
Cross-tissue convergent validity: van der Laan vs Lussier r 0.16–0.43.

**ADMET (sibling).** ESOL n=1,128 (solubility), Lipophilicity n=4,200 (logD), hERG n=8,257 (pIC50).
Murcko scaffold decomposition; CV schemes random / k-fold(LOO) / LOSO / LOSO* (drops heterogeneous
catch-all scaffold groups). LOO→LOSO* R² drop tracks **fragment-additivity**: non-additive
(ESOL +0.18 Ridge, up to +0.28 RF; hERG +0.10) vs additive (**logD ≈ 0, even −0.02**). "hERG model
reporting R²=0.48 by random split ≈ 0.37 on novel scaffolds — 22% relative overstatement." Note:
**no measurement-theory vocabulary anywhere**; framed as CV methodology.

---

## 5. Existing formal framing already in the repos (what the meta-paper builds on)

Only the **kinase** paper has explicit axiomatic framing, and it is exactly the right scaffold:

- Cites **Shannon's uniqueness theorem** (entropy is the unique function satisfying continuity,
  maximality, expansibility, additivity) as the *model* for what a well-founded measure looks like,
  and notes selectivity entropy "provides no analogous axiomatic justification."
- Cites **Arrow's impossibility** (no aggregator satisfies all reasonable axioms) — flagging that
  axiomatic analysis can yield *impossibility* rather than uniqueness.
- Cites **Singleton & Booth (2022) truth-discovery axiomatization** as "direct structural precedent."
- States the four properties are **necessary, explicitly not proven sufficient or unique**.
- Conclusion names the open goal verbatim: *"Establishing whether these four properties characterize a
  selectivity measure uniquely, as Shannon's axioms do for entropy, is the natural next step."*
- Has the ordinal-invariance premise already: analyses "concern rank correlations … invariant to a
  common monotone rescaling within a dataset."

**No repo uses representational-measurement-theory vocabulary** (scale type, representation/uniqueness
theorem, Krantz/Luce/Suppes) or category theory. The clocks work comes closest *implicitly*
("linearly related = same measurement on different scales"; uses ordinal Kendall κ over Pearson).

→ The measurement-theory layer and a representation/uniqueness-flavored result were *initially* thought
to be the meta-paper's genuine additions. **§6 corrects this:** that layer is established prior art
(Patil & Taillie; RMT; Szpilrajn), so it is *foundation to cite*, not novelty. The genuine additions are
the **protocol**, the **constitutive-vs-artifactual reproducibility test**, and the **controlled
empirical findings**. See `../paper/main.md`.

---

## 6. Controlled re-analyses, prior art, and decisions (added later — supersedes the "new-math" framing)

### 6.1 Why the reframe
Three hostile referees (math rigor / empirical support / novelty) were run against the initial
new-math draft. Consensus verdict: (a) the formalism is **classical and preempted** — it is the
partial-order theory of multi-indicator systems (**Patil & Taillie 2004**, *Environ Ecol Stat*
11:199–228; **Brüggemann & Patil 2011**, *Ranking and Prioritization for Multi-indicator Systems*,
Springer) plus RMT (Krantz–Luce–Suppes–Tversky) plus **Szpilrajn (1930)**; (b) the "open problem" of a
canonical minimal-commitment measure is **already solved** — the **average rank over linear extensions**
(De Loof, De Baets & De Meyer; interval-order variant for uncertainty); (c) the empirical signature
"agree at extremes / diverge in middle" is textbook **convergent validity / method variance** (Campbell
& Fiske 1959) and Bland–Altman range effect. Decision (user): the paper is a **methods/Perspective**
proposing a usable protocol and *showcasing how established math justifies its correctness* — not new
mathematics. Prior art becomes cited foundation. Lead conceptual novelty: separating **constitutive**
disagreement (concept genuinely silent, survives reproducibility) from **artifactual** (near-ties,
detection-floor, threshold boundaries) — a departure from the psychometric "method variance to remove."

### 6.2 Controlled re-analyses run in this work (numbers to cite)
Interpreter: `/Users/polina/miniforge3/bin/python`; R 4.5.2. Scripts in the session scratchpad.

**Kinase (Klaeger, n=222; gated-in = n_active>0 at pKd>6 ⇒ 206):**
- *Gate-cutoff robustness:* zero-active vs active disagreement (std of rank across 4 measures) =
  92 vs 25 (cutoff 5.5), **66 vs 23 (6.0)**, 55 vs 22 (6.5), 40 vs 21 (7.0). The detection-floor
  instability is robust and the high-disagreement compounds ARE the ungated ones (reproduces repo's
  30-config ~74 vs ~32). → G1/D1 reliability gate resolves the "extreme instability," not P3.
- *Concentration null test:* observed interior-concentration slope (Spearman disagreement~middleness)
  = 0.27 (consensus coord) / 0.11 (n_active coord); **independence-null band [0.18, 0.42] / [−0.13,
  0.14]** — observed does **NOT** exceed null. The "middle peak" is rank-boundary-compression artifact.
- *Pairwise discordance (21,115 pairs):* driven by **separation** — Δn_active=0 → 0.757 discordant;
  Δ∈[1,2] → 0.666; [3,5] → 0.467; [6,15] → 0.332; ≥16 → 0.186. Among near-ties (Δ≤2), position is flat
  (ends 0.740, mid 0.658, center 0.705); linear-prob model R² sep-only 0.164 → sep+midpos 0.166,
  midpos coef **−0.106**. → **Disagreement = near-ties; axis position adds nothing. No consequential
  middle.**
- *Average-rank canonical extension* (Bubley–Dyer sampler, 2 chains Spearman 0.997): consensus poset
  **37.2% incomparable pairs**; benchmark vs canonical — **Gini +0.955, entropy +0.950, S-score +0.869,
  ratio +0.796**; mean |rank−canonical| Gini 13.0 / entropy 13.6 / S-score 22.2 / ratio 29.0.
  → Gini/entropy best proxy the weight-free consensus aggregate; ratio commits most beyond consensus
  (quantitative "ratio is the outlier").
- *Candidate (softplus-entropy) — ISSUE FOUND (flag for kinase paper).* Using the repo's **exact**
  `candidate_measure.py`, the candidate is **anti-correlated** with the distribution family it should
  join: Spearman vs entropy −0.87, Gini −0.87, S-score −0.94, consensus-of-4 −0.86; robust across floors
  5.0→0.0. 93.6% of Klaeger entries are at the detection floor and softplus never zeroes them (pedestal
  T·log2≈0.69/kinase), so a truly selective compound (1 spike among 342 pedestals) looks near-uniform →
  high entropy → ranked *non*-selective (schematic: 1-target candidate-H≈8.4 near the panel max, many-target
  H lower; inversion direction confirmed on data — candidate_bench −0.87 vs entropy, un-gated +0.94 vs n_active). Not a sign
  typo. **Fix attempt reveals a fundamental tradeoff:** pedestal-subtraction restores orientation
  (entropy +1.00, S-score +0.94, canonical +0.996, p*=110, D4 100%) **but breaks D3** (worst-case −0.94,
  ≈ hard entropy). T-sweep: original is D3-robust(0.94–1.00) but inverted(−0.87); fixed is oriented(+1.00)
  but D3-fails(−0.94) — never both. So original's D3 pass is an artifact of the inversion-causing
  pedestal domination. **REPAIR FOUND** (`candidate_fix3.py`): anchoring the hinge at a **per-compound low
  quantile** (q10–median) instead of a fixed floor satisfies all four at once — orientation **+0.996**,
  D3 **+1.00** (shift-invariant, no baseline parameter), D4 **100%**, p*=110. (Per-compound *min* anchor
  is the weak +0.51 version; a robust quantile works because ~94% of entries at floor ⇒ q25 robustly
  estimates the inactive bulk.) **Two clean fixes for the kinase paper:** (i) fixed floor → per-compound-
  quantile anchor (small edit); (ii) the baseline-free **average-rank canonical extension**. Both
  weight/baseline-free. Scripts: `candidate_bench.py`, `candidate_fix.py`, `candidate_fix2.py`,
  `candidate_fix3.py`. NOT edited in the kinase repo — flagged to author.

**Smoking (per-cohort, from `data/scores_*.rds`):**
- *Per-cohort C3 agreement (thr = never_mean−2SD):* GSE50660 never .950 / current .909 / former .776;
  GSE42861 never .964 / current .760 / former .662. Former lowest in **both** cohorts (reproducibility,
  directional). Never-smokers have the HIGHEST agreement at every threshold multiplier (1.0–3.0) in both
  cohorts (rising 0.80→0.97 / 0.79→0.98 as the cutoff tightens); the LOWEST-agreement group is former at
  the 2SD headline and low/moderate multipliers, but current overtakes former at strict multipliers
  (≥2.5) in GSE42861 (current 0.62/0.56 vs former 0.77/0.82). That *which* non-never group looks worst
  shifts with the arbitrary cutoff is itself evidence the effect is a binary-threshold boundary
  phenomenon, not a stable concept-level property.
- *Threshold-free rank-spread* gives never > former > current (48.8/32.2/19.5; 69.9/58.6/41.0) — but
  this metric is **density/noise-floor confounded** (never-smokers cluster near the unexposed baseline
  where signatures rank by noise), so it is **not** used as evidence. Conclusion: the former-smoker
  "divergence" is a **binary-threshold boundary** effect, not concept-level middle concentration.
- *Data limitation:* neither GSE50660 nor GSE42861 carries a **continuous** exposure variable
  (pack-years/cigs-day) — only categorical never/former/current. The sharpest location analysis is
  therefore available only in kinase.

### 6.3 What survives as claimable (the paper's empirical core)
1. Disagreement is a **near-tie** phenomenon (separation-driven; position uninformative). Universal.
2. **Detection-floor / low-signal instability** is real and gate-necessary (G1). Robust.
3. **Cross-cohort reproducibility** of the disagreement structure — directionally supported, evidence
   thin (2 cohorts/domain; GSE42861 is an RA cohort → disease confound). Stated as a limitation.
4. A **weight-free canonical aggregate** (average-rank extension) exists, is computable, and ranks the
   existing measures by how far they stray from consensus.
- **Retired:** "concentration in the consequential middle" (P3) as a concept-level law.
- **Kept only as practical caveat:** classification disagreement peaks near decision thresholds
  (a boundary effect, explicitly not concept structure).

### 6.4 Kinase paper — no corrections needed
The live kinase manuscript (`selectivity/paper/sections/`) does **not** claim P3. Its abstract traces
instability to three structural sources ("noise below the detection threshold, near-tied top targets,
distribution-ratio disagreement for promiscuous compounds") — all **confirmed** by the controls above
(the near-tie result generalizes its "near-tied top targets"). The "intermediate promiscuity →
contradictory rankings" paragraph was **already commented out** by the author in
`introduction.tex:44–51` — i.e. the P3 overclaim was deliberately avoided. P3 lived only in the
meta-paper's synthesis and the original plan, and has now been removed there too.

*Update (candidate repaired + validated):* the author fixed the inversion (§6.2) with a **hard gate at
the fixed assay floor** (`selectivity/candidate_measure.py`; `candidate(P) = negH( 1[P>floor]·softplus((P−β)/T) )`)
and added `candidate_validation.py` — validating the candidate across **all four datasets** (orientation
−0.79 to −0.97, D4 100%, T-robust, cross-T agreement 0.999) with a per-dataset floor (incl. a
%-inhibition adaptation for Anastassiadis), plus an explicit inversion demo (un-gated +0.94 vs n_active).
The D3 story was corrected to be honest (emphasis-baseline robust; boundary pinned to the physical floor,
its sensitivity reported). This cross-dataset validation is the **evidence base for the three protocol
refinements** now folded into the meta-paper (external-anchor orientation → G2; apparatus-fixed vs
analyst-chosen params → G3; gate load-bearing for correctness → G1; see §2.5 table and main.md §2).

### 6.5 Cross-domain transfer (the "earns its generality" test)
Script `analysis/transfer.py`. For four domains (different concepts, measures, object types) — kinase
selectivity (4 measures, 222 cpds), serotonin selectivity (4, 297), epigenetic age (5 clocks, 1385
samples), smoking exposure (4 signatures, 1153 samples) — sign-align each domain's measures to a common
direction, then compute inter-measure **discordance** (min(a,k−a)/⌊k/2⌋) vs **normalized separation** on
the consensus axis. Result: the separation→discordance curve is **strongly monotone-decreasing in every
domain** (Spearman −0.99 to −1.00) — near-ties disagree, well-separated objects agree — and it
**transfers** (decile-curve RMSE): kinase↔clocks 0.10, kinase↔smoking 0.10, clocks↔smoking 0.14.
**Serotonin is the outlier** (RMSE 0.40–0.50): discordance stays high even for well-separated compounds
(0.26 at the far decile vs ≈0.00 elsewhere), consistent with its tiny 13-receptor panel + convention
issues (§3.2). Honest reading: the *near-tie law* (disagreement = near-ties) is the transferable,
domain-general regularity; the *magnitude/decay* is domain-specific (panel size, measure design). This
supports methodological generality (one protocol, one universal qualitative law) without claiming a
single quantitative curve fits all domains.

### 6.6 Orientation check + external-anchor validation across domains (backs refinements #1, #4)
Script `analysis/external_anchor.py`. For all four domains, with a crude external anchor per concept
(kinase/serotonin selectivity → n_active; clocks → chronological age; smoking → never<former<current):
- **Orientation check (G2/G1).** Every standard measure is correctly oriented against its anchor — kinase
  S-score −1.00 / entropy −0.94 / Gini −0.94 / ratio −0.32; serotonin all NEG; clocks all POS (+0.33 to
  +0.81 vs age); smoking all POS (+0.60 to +0.69). **Control:** the un-gated kinase candidate is +0.94 vs
  n_active → flagged MIS-ORIENTED. So the check generalizes and catches a real inversion. (N: kinase 222,
  serotonin 949 [≥4 receptors measured], clocks 1385, smoking 1085 [never/former/current only].)
- **Anchor validates the consensus (formal-spine §1.1).** Consensus–anchor Spearman: kinase +0.94,
  serotonin +0.82, clocks +0.65, smoking +0.68; consensus rises monotonically across anchor tertiles;
  agreement is **higher at the extremes than the middle** where computable (kinase 0.92 vs 0.56; clocks
  0.68 vs 0.36; serotonin/smoking middle undefined — discrete/tied anchor). Turns the §1.1 method from
  asserted to demonstrated.
Note (clocks): chronological age is also the training target of the position clocks, so their agreement
is expected — but that is precisely the role of a crude external anchor (the obvious independent proxy),
and the clocks are imperfect predictors (r 0.33–0.81), so there is genuine content.
