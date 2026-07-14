> **⚠ SUPERSEDED (kept for history).** This doc reflects the abandoned "new-mathematics" framing, in
> which the linear-extension/scale-type material was pitched as the meta-paper's novel contribution.
> Three referees + controlled re-analyses forced a re-scope to a **methods/Perspective** paper. Do not
> build on this file. Live docs: **`main.md`** (manuscript), **`formal-spine.md`** (Appendix A, corrected
> — the scale-type "ladder", the log⇒different-family argument, and the max-entropy "open problem" here
> are all WRONG and were removed), and **`../synthesis/repo-findings.md` §6** (controlled results,
> prior art, decisions). This file's still-valid residue: the G0–G5 desiderata mapping (now in `main.md`
> §3 / `formal-spine.md`) and the ADMET companion-not-fold-in decision (kept).

---

# Framework core — the formal spine (draft of §2, §3, §4, §7)  [SUPERSEDED — see banner above]

**What this is.** A first, real attempt at the part of the meta-paper that "lives or dies" (plan §8):
a formalism crisp enough that the paper has a *theorem*, not a shared vibe. Grounded in the verified
repo findings (`../synthesis/repo-findings.md`). Written to be argued with — every claim is flagged
as **[solid]**, **[defensible]**, or **[risky/open]**.

The one-sentence bet: **a biological measure is an order-preserving map (a homomorphism) from an
observed empirical structure into the reals; competing measures disagree exactly and only where that
empirical structure is a genuine partial order rather than a total one; the "families" are its
equivalence classes under the scale-type transformation group.** Everything the four instances found
is a corollary of that.

---

## §2 — The framework

### 2.1 The object of study

Fix a **latent concept** κ (selectivity, biological age, exposure, disease signature) and a
population/domain of **objects** Ω (compounds, individuals, samples). We do **not** assume κ is a
number. We assume only that the concept licenses **empirical comparative judgments**: for some pairs
of objects, domain knowledge (or an idealized perfect experiment) says *a is at least as κ as b*.
Write this relation `a ≽ b`.

**Definition (empirical structure).** The empirical structure for κ on Ω is the pair `(Ω, ≽)` where
`≽` is a **preorder** (reflexive, transitive) — crucially **not assumed total**. Some pairs are
`a ≽ b` and `b ≽ a` (equivalent, `a ~ b`); some are strictly ordered; and some are **incomparable**
(`a ∥ b`: neither `a ≽ b` nor `b ≽ a` is licensed). Incomparability is the load-bearing part.

**Definition (a measure).** A candidate measure is a function `f : Ω → ℝ`. It **represents** `(Ω, ≽)`
if it is order-preserving on the comparable pairs:
> `a ≽ b  ⟹  f(a) ≥ f(b)`   for all comparable `a, b`.  **(REP)**

This is a *representation condition* in the exact sense of representational measurement theory
(Krantz–Luce–Suppes–Tversky): `f` is a homomorphism from the empirical relational structure into
`(ℝ, ≥)`. **[solid — this is textbook; the contribution is applying it here, not inventing it.]**

### 2.2 Where the concept actually pins objects down: the comparability skeleton

The key structural object is **not** `f` but the preorder `≽` itself, specifically its comparable
part. Two facts about biological measurement concepts:

1. **At the extremes, `≽` is a total order.** Everyone agrees a compound that binds one kinase at
   1 nM and nothing else is more selective than one that binds 300 kinases uniformly; a 90-year-old
   is biologically older than a 20-year-old; a lifelong pack-a-day smoker is more exposed than a
   never-smoker. The comparative judgment is *forced*.
2. **In the middle, `≽` has incomparabilities.** Is a compound with two tight targets more or less
   selective than one with five moderate targets? Is a former smoker of 20 years more or less
   "exposed" (methylation-wise) than a light current smoker? The concept **does not decide**; only a
   *choice of measure* decides.

**Definition (a measure is a linear extension).** A total-order measure `f` satisfying (REP) induces
a total order that **extends** `≽`: it agrees with `≽` wherever `≽` is defined, and *breaks* the
incomparabilities `a ∥ b` one way or the other. Distinct valid measures are distinct **linear
extensions of the same partial order.** **[defensible — this is the central modeling claim; a
reviewer's fair objection is "is `≽` well-defined independent of the measures?" Answer in 2.3.]**

### 2.3 The three headline phenomena as theorems about linear extensions

Now P1–P3 from the plan stop being observations and become consequences.

**Theorem-shaped claim A (structured, not random — P1).** Any two measures `f, g` both satisfying
(REP) **agree on every comparable pair** of `≽` and can differ only on incomparable pairs. Hence their
disagreement is not noise: it is confined to, and diagnostic of, the incomparable set `I(≽) = {{a,b} :
a ∥ b}`. Reproducibility of the disagreement across cohorts (observed everywhere) is then *expected*,
because `I(≽)` is a property of the concept, not the sample. **[solid given the model.]**

**Theorem-shaped claim B (concentration in the consequential zone — P3).** The comparable pairs
concentrate at the extremes (2.2 fact 1) and the incomparable pairs concentrate in the middle
(2.2 fact 2). Therefore inter-measure disagreement is **necessarily** concentrated in the
intermediate region. This is exactly: former smokers (71% vs 96%/82%), intermediate-selectivity
compounds, mild disease, mid-range age. **The "consequential zone" is the image of `I(≽)`.**
**[defensible — the strongest single result in the paper; it *predicts* the location of disagreement
from the order structure alone. Needs one clean worked example, kinase is the candidate.]**

**Claim C (families — P2).** Group measures by the total order (equivalently, the linear extension)
they induce: `f ≈ g` iff they induce the same order on Ω. This is an equivalence relation; its classes
are the **families**. Within a family, measures agree on *all* pairs (comparable and incomparable) and
so differ only by a monotone reparameterization (§4). Between families, they disagree on incomparable
pairs — they chose different extensions. **[defensible, with the tension below.]**

> **⚠ Honest crack (repo-findings §3.2).** Empirically the *number* and *membership* of families is
> **not** domain-invariant: ratio is the kinase outlier, entropy the serotonin outlier; the smoking
> signatures form a *single* near-family. So Claim C must be stated as: *the partition into families is
> well-defined per (domain, dataset); that a small number of families exists — typically a
> low-disagreement bulk plus ≥1 categorically-different outlier — recurs; which measure is the outlier
> is scale- and panel-dependent.* Overclaiming a universal family structure is the fastest way to lose
> a referee. **[risky if overstated; safe if stated as above.]**

---

## §3 — General desiderata (G0–G5) as conditions on the representation

The domain axiom-lists (kinase D1–D4, clock D0–D4, methylation C1–C4) are **projections** of a
domain-independent schema. Each Gᵢ is a condition on the map `f` and/or the structure `(Ω,≽)`.
(Full cross-map table in `../synthesis/repo-findings.md` §2.5.)

- **G0 — Type/scale declaration.** Declare which empirical structure `f` represents and its scale type
  (the admissible-transformation group, §4). *Comparisons across declared types are meaningless* — a
  ratio-scale answer and an ordinal answer to different `≽`'s are not rival estimates of one number.
  → clock **D0** (explicit); kinase "ratio answers a different question" (implicit); the deepest axiom.
- **G1 — Reliability / domain gate.** `f` is defined only on the sub-domain where `≽` is observable
  (binding above assay floor; exposure above a sensitivity threshold). Off-domain, no representation
  is claimed. → kinase **D1**, methylation **C1**.
- **G2 — Monotonicity in the concept.** (REP) itself: `f` preserves `≽`. Domain forms: age↑⇒output↑
  (clock D1), exposure↑⇒score↑ (methylation C4), adding inactive target ⇏ more selective (kinase D4,
  a *directional* monotonicity). → the homomorphism condition.
- **G3 — Ordinal stability under nuisance.** The induced order is invariant to implementation-arbitrary
  parameters (baseline β, which off-target is "the" reference, off-manifold noise, adding a weak
  feature). Formally: if `θ` is a nuisance parameter, the order `≼_{f(·;θ)}` is locally constant in θ.
  → kinase **D2, D3**; clock **D2**; methylation **C2**.
- **G4 — Cross-instance reproducibility.** The induced order replicates across independent
  cohorts/panels/subsamples (because `≽` is a property of κ, per Claim A). Operationally: rank
  agreement across datasets; panel-size convergence. → clock **D3**, methylation **C3**, kinase
  panel-size stability.
- **G5 — Intervention consistency.** The represented order is preserved under the relevant
  operations/interventions, consistently with the declared type. → clock **D4**. (This is where the
  ADMET sibling *almost* connects — §6.)

**Deliverable for §3 of the paper:** the table above, proved-up as "each domain list is the image of
G0–G5 under forgetting the axioms that domain didn't test." E.g. kinase never states G4 as an axiom
(it tests it empirically via panel size); serotonin tests none of them (it inherits). **[solid — this
is concrete, doable, and is the "specialization shown explicitly" the plan asks for.]**

---

## §4 — Measurement-theory mapping: families = scale-type equivalence classes

This is the section that names *why* the disagreement is structured, via the **uniqueness** side of
representational measurement theory (the partner of the representation side in §2).

**Scale type = the admissible-transformation group.** In RMT, a representation `f` is determined only
up to a group `G` of transformations that preserve the empirical structure:
- **Ordinal** scale: `G` = all strictly increasing `φ: ℝ→ℝ`. Only the order is meaningful.
- **Interval** scale: `G` = affine `x ↦ ax+b`, `a>0`. Differences are meaningful.
- **Ratio** scale: `G` = linear `x ↦ ax`, `a>0`. Ratios are meaningful; there's a true zero.

**Mapping claim.** *Within a family* (Claim C: same induced order), two measures represent the **same**
empirical structure and therefore differ by an element of `G`. The **observed functional form of the
disagreement reveals the shared scale type:**
- Position clocks are related by `f ≈ 0.95·g + b`, **R²=0.991–0.995** — an **affine** map ⇒ they are
  the *same interval-scale measurement*, different units. The repo's own words: "essentially the same
  measurement on different scales." **[solid — this is a direct, quantitative confirmation of the RMT
  reading, and it's already in the data.]**
- Rate vs position clocks are related by a **signed logarithm** (R²=0.93 vs 0.75 linear;
  Weber–Fechner). A log is *not* in the affine group ⇒ rate and position are **different scale types /
  different `≽`'s**, not rival estimates. This is a *between-type* relation, formalizing D0. **[solid.]**
- Kinase "ratio is the consistent outlier" = ratio induces a different order (not `G`-related to the
  distribution family for any `G`) ⇒ it represents a **different empirical structure** — it answers
  "how many fold over the second-best target," a genuinely different `≽`. **[defensible.]**

**So the meta-claim of §4 is:** disagreement decomposes into (i) **within-family**, which is a *choice
of scale representative* and is harmless for ordinal claims but fatal if someone interprets interval/
ratio structure that isn't licensed; and (ii) **between-family**, which is a **scale-type / structure
mismatch** — the measures are homomorphisms of different empirical structures and *no* transformation
reconciles them. Calling them "competing definitions of one quantity" is the category error the whole
program diagnoses. **[this is the intellectual heart; it is defensible and it is genuinely new as a
unifying statement — no instance says it.]**

Category-theory framing (optional, plays to author's background, §7-adjacent): a measure is a functor
from the thin category `(Ω,≽)` (a preorder) to `(ℝ,≤)`; representation = functoriality (REP); families
= naturally-isomorphic functors; scale type = the automorphism group of the target object. This is a
*presentation* choice — the content stands without it — but it is the differentiator and can be the
paper's voice. **[defensible; use lightly, don't let it become the paper.]**

---

## §7 — The representation/uniqueness-flavored result (stated honestly, uniqueness open)

The plan (§7, §8) says a partial representation/uniqueness result is *non-negotiable* for this to be a
research article. Here is the one I think is true and provable-at-the-right-altitude.

> **Proposition (characterization up to scale type — target statement).**
> Let `(Ω,≽)` be the empirical structure for κ, with comparable part `C(≽)` and incomparable part
> `I(≽)`. Let 𝓕 be the class of measures satisfying G0–G4. Then:
> 1. **(Representation)** Every `f ∈ 𝓕` restricted to `C(≽)` is the same up to the scale-type group `G`
>    declared under G0; i.e. on comparable pairs, 𝓕 is a *single* equivalence class. Disagreement
>    within 𝓕 lives entirely in `I(≽)`.
> 2. **(Non-uniqueness / degrees of freedom)** The residual freedom is exactly the set of linear
>    extensions of `≽`; `|𝓕/≈|` equals the number of order-distinct extensions realizable by G0–G4
>    measures. G0–G4 therefore **characterize the admissible class, not a unique element.**
> 3. **(What uniqueness would need)** A unique `f` (Shannon-style) requires an *additional* structural
>    axiom that selects one extension — e.g. an aggregation/independence axiom (Shannon additivity +
>    expansibility for entropy; a concatenation operation for extensive measurement). Absent such an
>    axiom, uniqueness is *false*, and the honest output is a **canonical construction** (an existence
>    witness), not a uniqueness theorem.

**Why this is the right result:**
- It is **provable** at the stated altitude (parts 1–2 are essentially the definition of linear
  extension + the RMT uniqueness theorem; the work is stating them cleanly for this setting).
- It **explains every instance**: kinase's candidate is explicitly "a candidate, not unique" (part 3 —
  the softplus width T is a free design choice; T→0 collapses it) → it is precisely an *existence
  witness*. Clocks' "no test satisfies all five," methylation's "no signature satisfies all four" →
  the class 𝓕 is nonempty-but-not-singleton exactly as parts 1–2 predict.
- It is **honest** (uniqueness open, flagged as needing an extra axiom) — matching plan §8's demand to
  not overclaim, and matching the kinase paper's own stance verbatim.
- It **discharges the "review vs research article" risk**: parts 1–2 are a theorem; part 3 tells you
  what the constructive candidates are witnessing. That's a formal contribution, not a survey.

**Open problem to state as such (invites follow-up, signals rigor):** *For a given κ, what natural
aggregation axiom selects a unique linear extension of `≽`?* The kinase softplus-entropy candidate
conjecturally corresponds to a maximum-entropy / minimal-commitment extension; proving that would be a
genuine uniqueness theorem for one domain and is the obvious next paper. **[risky/open — do NOT claim
it's solved; claiming it is the failure mode. Stating it as the open problem is a strength.]**

---

## §6 decision — QSAR/ADMET: companion, not fold-in **[my recommendation]**

ADMET has no `≽`, no desiderata, no families; it's about *evaluation* schemes overstating truth, in CV
vocabulary (repo-findings §3.4). It *is* the same person's program and the serotonin LOSO collapse
literally motivated it — but folding it into the G0–G5 / linear-extension machinery would require
re-casting "an evaluation" as "a measure of generalization," which is a second formalism the paper
doesn't have room to build rigorously. **Recommendation: one tight subsection or a boxed "second axis"
+ companion citation**, framed as "measures disagree (this paper) and evaluations overstate (companion)
are two faces of how comp-bio deceives itself about measurement." Keeps scope discipline (plan §8).
The one real bridge worth stating: fragment-*additivity* in ADMET is itself a *scale-structure*
property (additive = extensive/concatenation structure = ratio-scale-friendly), which is the same RMT
vocabulary — so the companion citation can be made to rhyme with §4 in one sentence.

---

## Immediate open decisions for Polina (bring to the conversation)
1. **Is the linear-extension-of-a-partial-order the spine you want?** It's the load-bearing bet. It
   makes P1/P3 theorems and P2 a (carefully hedged) corollary, and it's category-theory-native.
   Alternative spines: pure Shannon-style axiomatic-uniqueness-per-domain (narrower, what kinase already
   does), or a looser Perspective piece (plan §8 fallback). I recommend the linear-extension spine.
2. **How hard to formalize?** Full theorem-with-proof (research article, aim *Patterns*/*Brief.
   Bioinform.*) vs. stated-precisely-but-semi-formal (Perspective, *PLOS Comp Biol*). Gated by whether
   Prop parts 1–2 can be made airtight — I believe they can.
3. **Does the family-instability crack (3.2) worry you?** I want to state it openly as a finding
   ("outlier identity is scale/panel dependent"), not bury it. Confirm you're comfortable with that.
4. **ADMET: companion citation** (my rec) vs a real §5 section?
5. **Prereq preprints:** kinase is live (ChemRxiv) but under revision for Mol. Inf.; clocks bioRxiv is
   "forthcoming" but the manuscript is **0-byte stubs** — it does not exist yet. The meta-paper leans
   on both. Is writing the clock manuscript a blocker we schedule, or do we cite the repo/README?
