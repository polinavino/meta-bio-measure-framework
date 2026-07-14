> **⚠ SUPERSEDED / HISTORICAL (do not treat as current).** This is the original planning brief. The
> project has since been (a) re-scoped from a "new-mathematics" paper to a **methods/Perspective** paper,
> and (b) corrected by controlled re-analyses. In particular, **P3 ("concentration in the consequential
> zone") has been RETIRED** — disagreement is a near-tie phenomenon, not a middle-concentration law (see
> `synthesis/repo-findings.md` §6). The formalism is **prior art to cite**, not a novel contribution.
> Current docs: `README.md`, `paper/main.md`, `paper/formal-spine.md`, `synthesis/repo-findings.md` §6.
> Read those; treat the claims below as the starting hypothesis, much of which did not survive.

---

# Meta-paper plan — "Structured disagreement of biological measurement definitions"  [SUPERSEDED]

**Author:** Polina Vinogradova
**Status:** planning / handoff document (SUPERSEDED — see banner)
**Purpose of this file:** This is a self-contained brief for a *future conversation* (or collaborator)
that will develop a synthesis / meta-paper. It assumes the reader has **not** seen the underlying
repos. Everything needed to pick up cold is here: the thesis, the recurring result, the theoretical
anchor, the four domain instances with links, a proposed paper structure, and — most importantly —
the honest gap that must be closed for this to be a *paper* and not just four case studies stapled
together.

---

## 1. One-line thesis

Across chemistry and biology, whenever multiple formal definitions claim to quantify a single latent
concept (kinase *selectivity*, *biological age*, *exposure*, a *disease signature*), the definitions
**disagree systematically rather than randomly**; the disagreement is **concentrated in the
intermediate/ambiguous region**, can be **organized into families by what the measure actually
tracks**, and can be **disciplined by a small set of desiderata** — of which **no existing measure
satisfies all**, but a constructive candidate can be proposed. This is *representational measurement
theory applied to biological measurement.*

## 2. Why this is one paper, not five

The author has independently produced the **same result** in four domains, plus a fifth,
sibling result about evaluation validity. The repetition is not redundancy — it is evidence of a
**general phenomenon** that has never been stated as such. The meta-paper's job is to state it once,
formally, and present the domain studies as instances.

Positioning note (from the strategy conversation): this paper is the author's **flagship identity
piece**. It should foreground the formal-methods background (category theory / formal verification /
measurement theory), not hide it. Target self-description: *"applied measurement theory for
computational biology."*

---

## 3. The recurring pattern (the general result to formalize)

In each domain the author found some subset of:

- **(P1) Structured, not random, disagreement.** Competing measures correlate highly on extreme
  cases and diverge in the middle. The divergence is reproducible across datasets/cohorts, so it is
  a property of the *definitions*, not noise.
- **(P2) Families.** Measures cluster into a small number of groups by the order relation they
  induce; within a family they nearly agree, between families they systematically disagree.
- **(P3) Concentration in the consequential zone.** Disagreement peaks exactly where classification
  matters most (former smokers; intermediate-selectivity compounds; mild disease; mid-range age).
- **(P4) Desiderata.** A short axiom list (reliability threshold, monotonicity, directional/rank
  stability, type/scale declaration, intervention consistency) that a *well-formed* measure should
  satisfy. **No existing measure satisfies all of them.**
- **(P5) A constructive candidate.** The author proposes a measure that satisfies more of the
  desiderata (e.g. entropy + reliability gate + smooth cutoff in the kinase work), while being
  explicit that *uniqueness* is left open.

## 4. Theoretical anchor (the missing formal spine)

The domain papers rediscover, ad hoc, ideas from **representational measurement theory**
(Krantz, Luce, Suppes, Tversky — *Foundations of Measurement*, vols I–III). The meta-paper should
make this explicit and use it as the formal backbone:

- Each candidate measure presupposes a **scale type** (ordinal vs interval vs ratio). Much of the
  "disagreement" is a **scale-type mismatch** — e.g., a ratio-family metric answers a different
  question than a distribution-family metric (this is literally the kinase finding: "ratio is the
  consistent outlier").
- The desiderata are informal **representation + uniqueness conditions**. The strongest possible
  contribution is a **representation/uniqueness-flavored result**: *given desiderata X on the
  admissible orderings, the class of admissible measures is characterized as Y* (even a partial or
  conjectural characterization is a real formal contribution and plays to the author's strengths).
- The author's **category-theory / formal-methods** background is the natural toolset here and is
  the differentiator. This is where to lean in.

**This section is the intellectual heart and the biggest risk (see §8).**

---

## 5. The domain instances

| # | Domain | Concept being measured | Competing definitions | Desiderata | Preprint | Repo |
|---|--------|------------------------|-----------------------|------------|----------|------|
| 1 | Kinase pharmacology | inhibitor **selectivity** | S-score, selectivity entropy, Gini, ratio | D1–D4 | ChemRxiv | kinase-selectivity-definitions |
| 2 | Epigenetics | **biological age** (clocks) | position / deviation / rate clocks | D0–D4 | bioRxiv (forthcoming) | epigenetic-clock-desiderata |
| 3 | Serotonin pharmacology | 5-HT2A/2B **selectivity** + biased agonism | S-score, entropy, Gini, ratio; bias factor | (shares kinase D1–D4) | — | psychedelic-selectivity |
| 4 | Epigenetic biomarkers | environmental **exposure / disease signature** | AHRR / EpiSmoke / Joehanes / EpiTob (smoking); Portales / Lussier×2 / van der Laan (FASD) | C1–C4 | — | methylation-biomarker-agreement |
| 5 (sibling) | QSAR / ADMET | model **generalization** (eval validity) | random / k-fold / LOSO CV schemes | — | — | admet-qsar-evaluation |

### Repo links
- **kinase-selectivity-definitions** — https://github.com/polinavino/kinase-selectivity-definitions
  (ChemRxiv preprint: https://doi.org/10.26434/chemrxiv.15001618/v1)
- **epigenetic-clock-desiderata** — https://github.com/polinavino/epigenetic-clock-desiderata
- **psychedelic-selectivity** — https://github.com/polinavino/psychedelic-selectivity
- **methylation-biomarker-agreement** — https://github.com/polinavino/methylation-biomarker-agreement
- **admet-qsar-evaluation** — https://github.com/polinavino/admet-qsar-evaluation

### Per-instance details (so the next conversation doesn't have to re-derive them)

**1. Kinase selectivity (D1–D4).** Four metrics compared across four datasets spanning three assay
technologies: Davis (68×433, competition binding), Klaeger (222×343, chemoproteomics),
Anastassiadis (178×300, functional % inhibition), Metz (704×172, pKi). Key results: **ratio is the
consistent outlier** (answers a different question than the distribution-family metrics: S-score /
entropy / Gini); entropy stabilizes above ~110 kinases, ratio needs >320. **Candidate measure:**
entropy with a reliability gate + smooth (not hard) activity cutoff; recovers full-panel ranking
from ~half a panel. Claim: "no existing formula satisfies all four" desiderata. This is the most
mature instance and the template for the formal treatment.

**2. Epigenetic clocks (D0–D4).** Argues different clocks measure fundamentally different things
("odometer vs speedometer"). Models healthy aging as a **canonical trajectory** (principal curve,
Hastie–Stuetzle) under an age-informativeness-weighted norm; decomposes each sample into position
(τ) + off-manifold residual. Clock **types**: position, deviation, rate. Five desiderata: type
declaration, monotonicity, directional stability, rank consistency, intervention type consistency.
Datasets: GSE40279 (Hannum, n=656), GSE87571 (n=729), 1,385 blood samples total. Results: position
clocks linearly related (R²=0.991–0.995); rate clocks relate to position clocks via a **logarithm**
(Weber–Fechner-like); **all clocks fail D3** (cell-type confounding, |r| up to 0.372); GrimAge is a
two-stage composite, structurally incomparable (near-orthogonal, cosine 0.03–0.11). Cites Klemera &
Doubal (2006) as closest precursor. "No existing test satisfies all five."

**3. Psychedelic / serotonin selectivity.** Applies the kinase selectivity metrics to 5-HT2A vs
5-HT2B (cardiac-safety motivation: 5-HT2B agonism → valvulopathy). ChEMBL, 13,584 compounds across
13 serotonin receptors. **S-score and Gini strongly anti-correlated (r = -0.682); entropy
uncorrelated with the others** — direct evidence for the "families" claim (P2). Also: LOO-CV
misleadingly optimistic vs LOSO (ties to instance 5). Finding: psilocin is non-selective,
"inconsistent with its common description as a '5-HT2A agonist.'" Generalizes the kinase entropy
result into a second pharmacology domain.

**4. Methylation biomarker agreement (C1–C4).** THIS repo. Framework: C1 sensitivity threshold, C2
score stability, C3 cross-cohort/cross-signature consistency, C4 monotonicity — explicitly built as
analogs of the kinase desiderata. **Smoking arm** (within-tissue): four signatures, two blood 450k
cohorts (GSE50660 n=464, GSE42861 n=689). Signatures agree on never-smokers (96%) and current (82%)
but diverge on **former smokers (71%)** — the intermediate zone (P3), consistent across cohorts.
**FASD arm** (cross-tissue reliability audit): dominant failure mode is **poor cross-tissue transfer**
(buccal signatures → blood AUC 0.68–0.79, don't track severity; blood-native episignature AUC
0.93–0.96). Notable rigor catch: the strong episignature (van der Laan 2025) is **non-independent**
(same lab / tissue / likely sample overlap with the eval cohorts) → treated as an upper bound.
Full writeup lives in `docs/writeup.md` of this repo.

**5. QSAR/ADMET evaluation validity (sibling theme).** Not about *definitions* disagreeing but about
*evaluations* lying: standard CV inflates apparent performance because models memorize scaffolds.
Murcko scaffold decomposition; four CV schemes (random, k-fold/LOO, LOSO, LOSO*). Key result: the
performance drop under honest splits depends on whether the endpoint is **fragment-additive** —
meaningful drop for non-additive (solubility, hERG: 0.10–0.28 R² units), ~0 for additive (logD).
Datasets: ESOL (1,128), Lipophilicity (4,200), hERG (8,257). This is a **second axis** of the same
program: "measures don't agree" (P1–P5) and "evaluations overstate truth" are two faces of *how
computational biology deceives itself about measurement.* Decide in §6 whether to fold this in or
cite it as companion.

---

## 6. Proposed paper structure

1. **Introduction** — the shared phenomenon; why "which definition is correct?" is the wrong
   question; roadmap.
2. **A framework for biological measurement definitions** — formalize P1–P5. Define: latent concept,
   family, structured disagreement, the desiderata schema. Connect each desideratum to a
   representation/uniqueness condition (§4).
3. **General desiderata** — the domain-independent axiom set that the domain-specific lists
   (D1–D4, D0–D4, C1–C4) are instances of. Show the specialization explicitly.
4. **Instances** — four subsections (kinase, clocks, serotonin, methylation), each ~2 pages: the
   concept, the competing measures, which desiderata hold/fail, the families, the concentration
   result. Keep these tight — they are evidence, the framework is the contribution.
5. **(Optional) The evaluation-validity corollary** — fold in the QSAR/ADMET result as a second
   axis, or cite as companion. Decision needed.
6. **Toward canonical measures** — the constructive side (P5): what a candidate that satisfies the
   desiderata looks like; the uniqueness question stated honestly as open.
7. **Discussion** — implications for how biomarkers/metrics should be reported (declare scale type;
   report the ambiguous zone; validate in-tissue); relation to *Foundations of Measurement*.

---

## 7. What is genuinely new (the contribution claim)

Not "we noticed the same thing four times." The defensible novelty is:

- **A domain-independent desiderata schema** for biological measurement definitions, with the four
  domain axiom-sets shown to be specializations of it.
- **The families / scale-type-mismatch account** of *why* the disagreement is structured (P2 + §4),
  rather than merely reporting that it is.
- **A partial representation/uniqueness result** disciplining the space of admissible measures
  (even conjectural — this is the formal-methods contribution and the reason it's *this* author's
  paper).

## 8. Honest gaps / risks (read before starting)

- **The unifying formalism must be real.** The single biggest failure mode is a paper that is four
  case studies with a shared vibe and no theorem. §4 (measurement theory anchor) and §7 (a
  representation/uniqueness-flavored result) are non-negotiable for this to be more than a review.
  If the general result can't be made crisp, consider instead a **Perspective/Opinion** piece
  (still valuable, lower bar) rather than a research article.
- **Venue fit is awkward** — it spans chem + epigenetics + measurement theory. Candidates:
  *Patterns* (Cell Press), *PLOS Computational Biology* (Perspective), *Briefings in
  Bioinformatics*, *GigaScience*, or a measurement-theory-friendly venue. bioRxiv/arXiv first
  regardless.
- **Two preprints are prerequisites.** The kinase (ChemRxiv, live) and clock (bioRxiv,
  "forthcoming") papers should be **posted and citable** before or with the meta-paper — it leans
  on them.
- **Scope discipline.** Four instances + optional fifth is already a lot. Resist adding more domains;
  depth on the formalism beats breadth of examples.
- **Keep it constructive.** Every "X is broken" must be paired with a criterion + candidate (P4+P5),
  or the paper reads as sniping. The author already does this — preserve it.

## 9. Concrete next steps for the picking-up conversation

1. Draft §2–§3 first (the framework + general desiderata) — this is where the paper lives or dies.
   Try to state one representation/uniqueness-style claim, however partial.
2. Confirm the measurement-theory mapping (§4): for each desideratum, name the corresponding
   scale-type / representation condition.
3. Decide the QSAR/ADMET fold-in (§5 vs companion citation).
4. Only then compress the four instances into evidence subsections.
5. Pick a target venue and matching format (research article vs Perspective) before polishing.

## 10. Author context (for the picking-up conversation)

- Formal-methods / category-theory / formal-verification background (Coq, Agda, Lean; EUTxO/Cardano
  formal work). This is the differentiator — foreground it.
- Goal is **impact + a distinctive public research identity**, remote, not a conventional job hunt.
  This paper is the flagship of that identity.
- Related repos not central to this paper but part of the same person's portfolio: rnaseq-pipeline,
  neoantigen-prediction-pipeline, esol-solubility-gnn (applied comp-bio; show range, not thesis).
