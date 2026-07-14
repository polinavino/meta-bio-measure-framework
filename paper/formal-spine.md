# Appendix A — Formal justification of the protocol

**Role of this appendix (revised).** This is *not* a claim of new mathematics. The protocol in the
main text rests on established results — representational measurement theory (Krantz, Luce, Suppes &
Tversky), Szpilrajn's extension theorem, and the partial-order theory of multi-indicator systems
(Patil & Taillie 2004; Brüggemann & Patil 2011), including the average-rank canonical extension
(De Loof, De Baets & De Meyer). This appendix states, in one place, *which established result
justifies which step of the protocol*, and fixes the object precisely so the steps are well-defined.
Every result here is classical or definitional; the paper's contribution is the protocol and the
controlled empirical findings (main text §4), not the theorems.

Correctness note: an earlier draft of this appendix contained (i) a scale-type "ladder" claiming
families are interval/ratio classes, (ii) an argument that a logarithmic relation implies different
families, (iii) a "maximum-entropy extension" open problem, and (iv) a consensus-bounds statement in
the wrong direction. All four were wrong and are corrected/removed below; the reasons are recorded so
the errors are not reintroduced.

---

## 1. The object

Fix a latent concept κ (selectivity, biological age, exposure, disease signature) and a finite set of
objects Ω. Domain knowledge licenses **comparative judgments** `a ≽ b` ("a is at least as κ as b").

**Definition 1.** The empirical structure is `(Ω, ≽)` with `≽` a **preorder** (reflexive, transitive),
**not** assumed total. Derived: `a ~ b` (both ways), `a ≻ b` (strict), `a ∥ b` (incomparable — neither
direction licensed). `C(≽)` = comparable pairs, `I(≽)` = incomparable pairs.

**Definition 2.** A *measure* is `f : Ω → ℝ`. It is **admissible** if `a ≽ b ⟹ f(a) ≥ f(b)` and
`a ≻ b ⟹ f(a) > f(b)` (an order-preserving homomorphism into `(ℝ,≥)` — the representation condition of
representational measurement theory). `𝓕(≽)` = admissible measures.

**Definition 3.** `f` induces the total preorder `a ⪰_f b :⟺ f(a) ≥ f(b)`. `f ≈ g :⟺ ⪰_f = ⪰_g`;
classes are **families**. A total preorder `⪰` is a **linear extension** of `≽` if `a ≽ b ⟹ a ⪰ b`.
On finite Ω, `f ↦ ⪰_f` maps `𝓕(≽)` onto the linear extensions of `≽`, factoring through `≈`.

This is exactly the "multiple indicators → partial order → linear extensions" setup of Patil & Taillie
(2004); we adopt it, we do not claim it.

### 1.1 How `≽` is estimated in practice (consensus), stated in the correct direction

`≽` is not observed. Given the set 𝓜 of measures domain experts actually use, define the **consensus
order** `a ≽_𝓜 b :⟺ f(a) ≥ f(b) ∀ f ∈ 𝓜`. Under the assumption that every reasonable measure is
admissible for the true `≽` (`a ≽ b ⟹ ∀f, f(a) ≥ f(b)`), we get **`≽ ⊆ ≽_𝓜`**. Consequences (corrected
from the earlier draft, which stated these backwards):

- `C(≽) ⊆ C(≽_𝓜)`: consensus **over-declares** comparability. `C(≽_𝓜)` is an **upper** bound on the
  truly comparable pairs; `I(≽_𝓜)` a **lower** bound on the truly incomparable pairs.
- Failure mode: measures may agree on a pair *by coincidence* though `κ` leaves it incomparable, so
  `≽_𝓜` can wrongly certify comparability.
- Guard: **reproducibility across independent cohorts** (protocol step 3). A comparability that
  replicates across cohorts is unlikely to be coincidental; a disagreement that replicates is a
  property of κ, not the sample. Reproducibility, not the consensus construction, is what licenses
  reading `I(≽_𝓜)` as concept-level incomparability.

---

## 2. Lemma (justifies the consensus/disagreement step) — analytic, not a discovery

**Lemma 1.** If `f, g ∈ 𝓕(≽)` rank a pair in opposite strict directions (`f(a)>f(b)`, `g(a)<g(b)`),
then `a ∥ b`. Equivalently, admissible measures agree on the strict order of every comparable pair;
all inter-measure ranking reversals lie in `I(≽)`.

*Proof.* If `a ≽ b`, (REP) on `g` gives `g(a) ≥ g(b)`, contradiction. If `b ≽ a`, (REP) on `f` gives
`f(b) ≥ f(a)`, contradiction. So `a ∥ b`. ∎

**Status (honest).** This is the contrapositive of admissibility — analytic, with no empirical content
on its own; for the *estimated* order `≽_𝓜` it is true by construction (the reversal set *is* the
disagreement set). Its role is only to justify **protocol step 2**: the pairs on which reasonable
measures disagree are exactly the candidate incomparabilities to examine. The empirical weight sits in
steps 3–4 (reproducibility and the controlled location analysis), not here.

---

## 3. Families are ORDINAL classes (corrected — no interval/ratio ladder)

**Lemma 2.** If `f, g ∈ 𝓕(≽)` and `f ≈ g`, there is a strictly increasing `φ` with `g = φ∘f`;
conversely any strictly increasing `φ` gives `φ∘f ≈ f`. So a family is an equivalence class under
strictly monotone reparameterization — an **ordinal** class (Stevens; RMT).

*Proof.* Define `φ(f(a)) := g(a)`; well-defined and strictly increasing because `f ≈ g`. Converse
immediate. ∎

**Corrections recorded (do not reintroduce):**
- The earlier "interval/ratio ladder" (claiming affine-relatedness proves an interval scale) is
  **removed.** Interval/ratio structure requires a difference/concatenation empirical structure with
  its own axioms; none is established for these concepts, and Def 1 supplies only order. Two admissible
  measures being affine-related on one dataset is fully explained by two *ordinal* measures that happen
  to be linearly related there. Empirically this is expected when both are regressions on a shared
  external target (e.g. two epigenetic clocks both trained on chronological age), so it is **not**
  evidence of a shared latent interval structure.
- The claim "a logarithmic relation ⟹ different family" is **wrong and removed**: a signed logarithm is
  strictly increasing, so by Lemma 2 it keeps two measures in the **same** family. Family membership is
  diagnosed by **rank reversals** (whether `⪰_f = ⪰_g`), not by the regression form relating two
  measures. Descriptively, a non-affine monotone relation (e.g. rate vs position clocks) is worth
  reporting, but as an observation about two measures, not a scale-type proof.

Consequence for the protocol: step 5 clusters measures into families by their induced order (rank
agreement), and reports the relation between families descriptively — it does **not** assign scale
types.

---

## 4. Non-uniqueness (justifies "map the class, don't crown a winner") — classical

**Proposition (Szpilrajn + folklore).** On finite `(Ω,≽)`: `𝓕(≽) ≠ ∅`; `𝓕(≽)/≈` is in bijection with
the linear extensions of `≽`; and `|𝓕(≽)/≈| = 1` **iff `≽` is total**. If `I(≽) ≠ ∅` there are ≥2
order-distinct admissible measures, and no requirement that only entails order-preservation (the
desiderata G1–G4) can single one out.

*Proof.* Szpilrajn's extension theorem gives a linear extension; on a finite set it has a real
representation, which is admissible. A finite poset has a unique linear extension iff it is a chain;
if `a ∥ b`, both `a≻b` and `b≻a` extensions exist. ∎

**Status (honest).** Entirely classical (Szpilrajn 1930; standard order theory). We use it, by name, to
justify **protocol step 6**: because no admissible measure is canonical whenever incomparabilities
exist, "benchmark to crown the single correct measure" is mis-specified; the correct output is the
comparability skeleton plus the ambiguous (incomparable) set. This is the same conclusion reached in
social choice (Arrow) and truth-discovery axiomatics (Singleton & Booth 2022) — cited, not rederived.

---

## 5. The canonical aggregate: average rank over linear extensions (justifies step 7)

When a single summary is unavoidable, the protocol recommends the **average rank of each object over
all linear extensions** of the consensus poset (Patil & Taillie 2004; De Loof, De Baets & De Meyer
2008/2011), optionally the interval-order variant when scores carry uncertainty. It is **weight-free**
and commits no further than the consensus order forces — the principled resolution of the non-uniqueness
above. (The earlier draft posed a "maximum-entropy extension" as an open problem; this was ignorance of
the literature — average-rank already *is* the minimal-commitment canonical choice. Removed.)

**(i) any admissible closed-form vs (ii) the class centroid.** By the Theorem, the orderings induced by
the admissible measures `𝓕(≽)` are *exactly* `Lin(≽)`. So:
- **(i)** a single closed-form `f ∈ 𝓕(≽)` realizes *one* element of `Lin(≽)` — it satisfies the
  properties but makes an arbitrary commitment on each incomparable pair (which extension to pick).
- **(ii)** the average rank `r̄(x) = 𝔼_{⪰ ∈ Lin(≽)}[rank_⪰(x)]` is the **centroid** of `Lin(≽)`: it
  agrees with every `f ∈ 𝓕(≽)` on the comparable pairs (Prop 1) and, on each incomparable pair `{a,b}`,
  returns the fraction of extensions ranking `a` above `b` (the mutual-ranking probability) rather than
  forcing a side. It privileges no tie-breaking and has no free parameter.

Hence **(ii) dominates (i)** in the minimal-commitment sense *whenever `Lin(≽)` can be explored* — which,
since exact enumeration is #P-complete (Brightwell & Winkler 1991), means *approximately sampled* (Bubley
& Dyer 1999). This is a conditional dominance: the capability to sample the admissible class is exactly
what converts "pick one arbitrary formula" into "report the consensus of all of them." Two honest limits:
(ii) is canonical *relative to* `≽` (i.e. to the chosen base set `𝓜`), and average-rank is the scalar
summary of the fuller mutual-ranking-probability matrix. (i) remains preferable for per-object scoring in
isolation or when an interpretable closed form is needed.

**Demonstration (kinase, Klaeger, this work).** Consensus poset over the 206 reliability-gated
compounds (a ≽ b iff all four measures agree): **37.2% of pairs incomparable** — substantially partial.
Average rank estimated by a Bubley–Dyer uniform linear-extension sampler (two independent chains,
Spearman 0.997). Benchmarking each existing measure's ranking against the canonical average-rank:

| measure | Spearman vs canonical | mean |rank − canonical| |
|---|---|---|
| Gini | **+0.955** | 13.0 |
| entropy | **+0.950** | 13.6 |
| S-score | +0.869 | 22.2 |
| ratio | +0.796 | 29.0 |

So Gini/entropy are the best single proxies for the weight-free consensus aggregate; **ratio commits
most beyond consensus** — a quantitative version of "ratio is the outlier."

**⚠ Candidate-measure issue found (flag for the kinase paper).** Benchmarking the repo's *exact*
`candidate_measure.py` softplus-entropy candidate on Klaeger, it is **anti-correlated** with the
distribution family it is meant to join: Spearman(candidate, entropy) = −0.87, vs Gini −0.87, S-score
−0.94, consensus-of-4 −0.86; robust across floors 5.0→0.0. Cause: 93.6% of Klaeger entries sit at the
detection floor and softplus never zeroes them (pedestal T·log2 ≈ 0.69 per kinase), so a genuinely
selective compound (one spike among 342 pedestals) looks maximally uniform → high entropy → the
candidate calls it *non*-selective (worked example: 1-target H≈8.39 > 50-target H≈7.92, i.e. it ranks
the promiscuous compound as more selective). Not a sign typo — the pedestal dominates. **Unaffected:**
the candidate's panel-size convergence (p*≈110–170) and D4 monotonicity, which concern self-consistency,
not agreement direction.

**Deeper finding (a tradeoff, not a one-line fix).** Attempting the obvious fix — subtract the floor
pedestal so inactive kinases contribute ≈0 — *does* restore correct orientation (Spearman with entropy
+1.00, S-score +0.94, canonical +0.996; p*=110; D4 100%) **but breaks D3**: its baseline-robustness
worst-case Spearman falls to −0.94, because pedestal subtraction makes it ≈ the hard-cutoff entropy,
which itself fails D3. A T-sweep confirms the tradeoff is fundamental: across every smoothing width the
original is D3-robust (0.94–1.00) but inverted (−0.87), and the fixed version is correctly oriented
(+1.00) but D3-failing (−0.94) — **never both.** So the original candidate's D3 pass is an *artifact of
the same pedestal-domination that inverts it.*

**Repair found (satisfies all four).** The failure is the *arbitrary fixed floor*. Anchoring the softplus
hinge at a **per-compound low quantile** of the profile (q10–median of the compound's own kinase values)
instead of a fixed floor removes the baseline parameter entirely, and delivers all four desiderata at
once (`analysis/candidate_fix3.py`): orientation vs the consensus selectivity order **+0.996** (vs the
original −0.87), **D3 satisfied by construction** (shift-invariant, +1.00 — there is no baseline
parameter to vary), **D4 100%** monotone, and **p*=110** (fast, matching entropy). Intuition: since ~94%
of Klaeger entries sit at the floor, a per-compound q25/median robustly *estimates* the inactive bulk, so
the measure recovers correct hard-entropy orientation while being baseline-free. (Anchoring at the
per-compound *min* is the weak +0.51 version; a robust quantile is what works.) **Two clean resolutions
therefore exist:** (i) this minimally-modified closed-form candidate (fixed floor → per-compound-quantile
anchor), and (ii) the **baseline-free average-rank canonical extension** (§5). Both are correct and
weight/baseline-free; the average-rank is more principled (no anchor choice at all), the anchored
closed-form is a smaller edit to the existing kinase candidate. Scripts: `analysis/candidate_bench.py`,
`candidate_fix.py`, `candidate_fix2.py`, `candidate_fix3.py`. Not edited in the kinase repo — flagged to
the author to verify.

---

## 6. What the empirical work actually establishes (main text §4) — the real contribution

The theorems above are scaffolding. The controlled findings are the contribution; each is stated with
its honest scope.

1. **Disagreement is a near-tie phenomenon.** Kinase pairwise discordance is driven by *separation* on
   the selectivity axis: equal-n_active pairs are 76% discordant, well-separated pairs (Δn_active ≥ 16)
   only 19%. Axis **position adds nothing** beyond separation (added R² ≈ 0.002; middle-position
   coefficient −0.11). → Objects close in κ get ordered inconsistently; where near-ties fall is a
   density property, **not** a concept-level "consequential middle."
2. **No "consequential middle" concentration (P3 retired).** A naive partial Spearman of disagreement
   on middle-position *is* significant (0.315, p=4e-6, controlling n_active/max_act) — which is exactly
   the trap. But that control does not remove rank-boundary compression; the proper boundary-compression
   **independence null** does, and against it the apparent interior peak vanishes: observed slope
   (0.11–0.27) lies inside the null band (up to 0.42), in fact at/below its mean. So the interior peak is
   a mechanical artifact, not concept structure. In smoking, the former-smoker divergence (per-cohort
   classification agreement 78% / 66%; the source repo's pooled figure is 71%) is a
   **binary-threshold** effect (robust to the cutoff, but a boundary phenomenon); a threshold-free
   metric does **not** put the peak at the intermediate group, and is itself density-confounded — so no
   directional concentration law is claimed.
3. **Reliability-gate / detection-floor instability is real and robust.** Kinase zero-active compounds
   carry the large instability (disagreement ~40–92 vs ~21–25; ≥~2× the active level), stable across every pKd active cutoff
   (5.5–7.0). These are exactly what the gate (G1) excludes; excluding them is necessary, not optional.
4. **Reproducibility across cohorts (G4) is directionally supported.** In smoking the classification-
   agreement ordering (never > current > former) reproduces in both cohorts separately (GSE50660
   .95/.91/.78; GSE42861 .96/.76/.66). Magnitudes differ; one cohort is a rheumatoid-arthritis study
   (disease confound) — stated as a limitation.

**Relation to prior framings (cite, delimit):** the situation is *convergent validity* (Campbell &
Fiske 1959) — multiple operationalizations of one construct that disagree. Psychometrics treats mid-
range divergence as **method variance to remove**; our contribution is to separate divergence that is
**constitutive** (concept-level incomparability, surviving reproducibility) from divergence that is
**artifactual** (near-ties, detection-floor noise, threshold boundaries). The "agree-at-extremes /
diverge-in-middle" pattern alone is also the textbook Bland–Altman proportional-bias/range effect;
what is *not* generic is the reproducibility test that promotes a disagreement to a property of κ.

---

## 7. Proved vs assumed vs open

- **Classical/definitional (used, not claimed):** Lemma 1, Lemma 2, the Szpilrajn non-uniqueness
  proposition, the average-rank canonical extension.
- **Assumed:** `≽` exists as a property of κ and reasonable measures respect it (operationalized as the
  reproducible consensus order). One substantive assumption; argued, not proved.
- **Empirical (this work's contribution):** §6.1–§6.4, with the stated scope and limitations.
- **Deferred, not open problems:** benchmarking the softplus-entropy candidate against the canonical
  extension (needs the repo implementation); continuous-exposure concentration analysis in smoking
  (the public cohorts lack a continuous exposure variable — a data limitation, flagged).
