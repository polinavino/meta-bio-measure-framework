# Appendix A — Formal justification of the protocol

**This appendix claims no new mathematics.** The protocol in the main text rests on established
results: representational measurement theory (Krantz, Luce, Suppes & Tversky), Szpilrajn's extension
theorem, and the partial-order theory of multi-indicator systems (Patil & Taillie 2004, Brüggemann &
Patil 2011), including the average-rank canonical extension (De Loof, De Baets & De Meyer). This
appendix states, in one place, *which established result justifies which step of the protocol*, and
fixes the object precisely so the steps are well-defined. Every result here is classical or
definitional. The paper's contribution is the protocol and the controlled empirical findings (main
text §4–§5), not the theorems.

Correctness note: an earlier draft of this appendix contained (i) a scale-type "ladder" claiming
families are interval/ratio classes, (ii) an argument that a logarithmic relation implies different
families, (iii) a "maximum-entropy extension" open problem, and (iv) a consensus-bounds statement in
the wrong direction. All four were wrong and are corrected or removed below, with the reasons recorded
so the errors are not reintroduced. A later pass corrected four more: the average rank called a
"centroid of the admissible class", the average rank credited with returning the mutual-rank
probability, #P-completeness predicated of exact average ranks rather than of counting, and Szpilrajn
credited with the whole Proposition rather than with non-emptiness alone.

---

## 1. The object

Fix a latent concept κ (selectivity, biological age, exposure, variant impact, inflammatory burden)
and a finite set of objects Ω. Domain knowledge licenses **comparative judgments** `a ≽ b` ("a is at
least as κ as b").

**Definition 1.** The empirical structure is `(Ω, ≽)` with `≽` a **preorder** (reflexive, transitive),
**not** assumed total. Derived: `a ~ b` (both ways), `a ≻ b` (strict), `a ∥ b` (incomparable — neither
direction licensed). `C(≽)` = comparable pairs, `I(≽)` = incomparable pairs.

**Definition 2.** A *measure* is `f : Ω → ℝ`. It is **admissible** if **both** `a ≽ b ⟹ f(a) ≥ f(b)`
**and** `a ≻ b ⟹ f(a) > f(b)`, an order-preserving homomorphism into `(ℝ,≥)`, which is the
representation condition of representational measurement theory. `𝓕(≽)` = admissible measures. Both
clauses are required. The weak clause alone permits a constant map, which preserves no strict
comparison and must not count as a measure of κ.

**Definition 3.** `f` induces the total preorder `a ⪰_f b :⟺ f(a) ≥ f(b)`, and `f ≈ g :⟺ ⪰_f = ⪰_g`,
whose classes are the **families**. A total preorder `⪰` is a **linear extension** of `≽` if
`a ≽ b ⟹ a ⪰ b` and `a ≻ b ⟹ a ≻_⪰ b`, that is, if it preserves strict comparisons as well as weak
ones. On finite Ω, `f ↦ ⪰_f` maps `𝓕(≽)` onto the linear extensions of `≽`, factoring through `≈`.
The strict clause is what makes "onto" true. Without it, a total preorder that merely weakened some
`a ≻ b` to a tie would count as a linear extension while being induced by no admissible `f`, since
Definition 2 forces `f(a) > f(b)`. Note also that linear extensions here are tie-allowing total
preorders, not the antisymmetric linear orders of the Brightwell–Winkler and Bubley–Dyer results in §5.

This is exactly the "multiple indicators → partial order → linear extensions" setup of Patil & Taillie
(2004). We adopt it, we do not claim it.

### 1.1 How `≽` is estimated in practice (consensus), stated in the correct direction

`≽` is not observed. Given the set 𝓜 of measures domain experts actually use, define the **consensus
order** `a ≽_𝓜 b :⟺ f(a) ≥ f(b) ∀ f ∈ 𝓜`. Under the assumption that every reasonable measure is
admissible for the true `≽` (`a ≽ b ⟹ ∀f, f(a) ≥ f(b)`), we get **`≽ ⊆ ≽_𝓜`**. Consequences (corrected
from the earlier draft, which stated these backwards):

- `C(≽) ⊆ C(≽_𝓜)`: consensus **over-declares** comparability. `C(≽_𝓜)` is an **upper** bound on the
  truly comparable pairs, and `I(≽_𝓜)` a **lower** bound on the truly incomparable pairs.
- Failure mode: measures may agree on a pair *by coincidence* though `κ` leaves it incomparable, so
  `≽_𝓜` can wrongly certify comparability.
- Guard: **reproducibility across independent cohorts** (protocol step 3). A comparability that
  replicates across cohorts is unlikely to be coincidental, and a disagreement that replicates is a
  property of κ, not the sample. Reproducibility, not the consensus construction, is what licenses
  reading `I(≽_𝓜)` as concept-level incomparability.
- **(External anchor as a second, independent check.)** Reproducibility guards against sampling
  coincidence but is still *internal* to the measures `𝓜`. A crude, theory-light **external anchor** `A`
  for the concept (number of active targets for selectivity, chronological age for clocks, self-reported
  exposure for smoking, clinical labels for variant effect, mortality for inflammation) gives an
  *independent* check. The consensus order `≽_𝓜` should agree with the order induced by `A` at the
  extremes, where both are unambiguous. Agreement there is external-validity evidence that `≽_𝓜` tracks
  κ rather than a shared artifact of the measures, and disagreement there is a red flag that the
  measures are collectively mis-oriented (the failure mode the step-1 orientation check targets, main
  text §2, §4.1). The anchor is too crude to *define* `≽`, being one more coarse measure, but it is
  exactly the right tool to *falsify* a mis-oriented consensus. **Demonstrated in the four
  anchor-scored domains** (`analysis/external_anchor.py`): every standard measure is correctly oriented
  against its anchor (n_active / chronological age / exposure ordinal), the mis-oriented un-gated
  candidate is caught (+0.94 vs n_active, wrong sign), and the consensus agrees with the anchor
  (Spearman +0.65 to +0.94), more strongly at the extremes than in the middle where the anchor is
  continuous enough to split (kinase 0.92 vs 0.56, clocks 0.68 vs 0.36). The check also has teeth in
  the later domains: it flagged all four TP53 proliferation screens as mis-oriented under the naive
  fitness convention (main text §4.5). **Caveats that keep this honest:** the extremes-versus-middle
  split is undefined for the discrete anchors (serotonin, smoking), so that leg rests on kinase and
  clocks. Two anchors are not independent of the measures they check. For the *position clocks* the
  anchor `A` = chronological age is their **training target**, and for the sepsis signatures the anchor
  is the mortality outcome the signatures were built against, so both of those checks are confirmatory
  rather than independent. The genuinely independent anchors are kinase, serotonin, smoking and the
  TP53 clinical labels, and the "independent check" framing applies fully only to those.

---

## 2. Disagreement between admissible measures lies in the incomparable set (justifies step 2)

**Lemma 1.** If `f, g ∈ 𝓕(≽)` rank a pair in opposite strict directions (`f(a)>f(b)`, `g(a)<g(b)`),
then `a ∥ b`. Equivalently, admissible measures agree on the strict order of every comparable pair, and
all inter-measure ranking reversals lie in `I(≽)`.

*Proof.* If `a ≽ b`, admissibility of `g` gives `g(a) ≥ g(b)`, a contradiction. If `b ≽ a`,
admissibility of `f` gives `f(b) ≥ f(a)`, a contradiction. So `a ∥ b`. ∎

**Status (honest).** This is the contrapositive of admissibility, analytic and with no empirical
content on its own. For the *estimated* order `≽_𝓜` it is true by construction (the reversal set *is*
the disagreement set). Its role is only to justify **protocol step 2**: the pairs on which reasonable
measures disagree are exactly the candidate incomparabilities to examine. The empirical weight sits in
steps 3–4 (reproducibility and the controlled location analysis), not here.

---

## 3. Families are ordinal classes and nothing stronger (corrected, no interval/ratio ladder)

**Lemma 2.** If `f, g ∈ 𝓕(≽)` and `f ≈ g`, there is a **strictly increasing** `φ` with `g = φ∘f`, and
conversely any **strictly increasing** `φ` gives `φ∘f ≈ f`. So a family is an equivalence class under
strictly monotone reparameterization, an **ordinal** class (Stevens, RMT). The hypothesis is strict
increase, not mere monotonicity: a weakly increasing `φ` can collapse distinct values of `f` into one
value of `φ∘f`, which destroys a strict comparison and leaves the family.

*Proof.* Define `φ(f(a)) := g(a)`. It is well-defined and strictly increasing because `f ≈ g`. The
converse is immediate. ∎

**Corrections recorded (do not reintroduce):**
- The earlier "interval/ratio ladder" (claiming affine-relatedness proves an interval scale) is
  **removed.** Interval/ratio structure requires a difference/concatenation empirical structure with
  its own axioms. None is established for these concepts, and Def 1 supplies only order. Two admissible
  measures being affine-related on one dataset is fully explained by two *ordinal* measures that happen
  to be linearly related there. Empirically this is expected when both are regressions on a shared
  external target (e.g. two epigenetic clocks both trained on chronological age), so it is **not**
  evidence of a shared latent interval structure.
- The claim "a logarithmic relation ⟹ different family" is **wrong and removed**. A signed logarithm is
  strictly increasing, so by Lemma 2 it keeps two measures in the **same** family. Family membership is
  diagnosed by **rank reversals** (whether `⪰_f = ⪰_g`), not by the regression form relating two
  measures. Descriptively, a non-affine but strictly increasing relation (e.g. rate versus position
  clocks) is worth reporting, but as an observation about two measures, not a scale-type proof.

Consequence for the protocol: step 5 clusters measures into families by their induced order (rank
agreement), and reports the relation between families descriptively — it does **not** assign scale
types.

---

## 4. No measure is canonical once incomparabilities exist

**Proposition (non-uniqueness).** On finite `(Ω,≽)`: `𝓕(≽) ≠ ∅`, `𝓕(≽)/≈` is in bijection with
the linear extensions of `≽`, and `|𝓕(≽)/≈| = 1` **iff `≽` is total**. If `I(≽) ≠ ∅` there are ≥2
order-distinct admissible measures, and no requirement that only entails order-preservation (the
desiderata G1–G4) can single one out.

*Proof.* Szpilrajn's extension theorem supplies a linear extension, which on a finite set has a real
representation, so `𝓕(≽) ≠ ∅`. A finite poset has a unique linear extension iff it is a chain, and
if `a ∥ b` then both an `a≻b` and a `b≻a` extension exist. ∎

**Status (honest).** Entirely classical. **Attribution, stated precisely:** Szpilrajn (1930) gives the
**non-emptiness** of `𝓕(≽)` and nothing more. The bijection between families and linear extensions,
and the "unique iff total" criterion, are standard order theory rather than Szpilrajn's theorem, so we
do not cite him for them. We use the Proposition, by name, to justify **protocol step 6**: because no
admissible measure is canonical whenever incomparabilities exist, "benchmark to crown the single
correct measure" is mis-specified, and the correct output is the comparability skeleton plus the
ambiguous (incomparable) set. This is the same conclusion reached in social choice (Arrow) and
truth-discovery axiomatics (Singleton & Booth 2022), cited, not rederived.

---

## 5. The average rank over linear extensions is the minimal-commitment aggregate (justifies step 7)

When a single summary is unavoidable, the protocol recommends the **average rank of each object over
all linear extensions** of the consensus poset (Patil & Taillie 2004, De Loof, De Baets & De Meyer
2008, 2011), optionally the interval-order variant when scores carry uncertainty. It is **weight-free**
and commits no further than the consensus order forces, which makes it the principled resolution of the
non-uniqueness above. (The earlier draft posed a "maximum-entropy extension" as an open problem. That
was ignorance of the literature, since average-rank already *is* the minimal-commitment canonical
choice. Removed.)

**(i) any admissible closed form versus (ii) the barycenter of the rank vectors.** By the Proposition,
the orderings induced by the admissible measures `𝓕(≽)` are *exactly* `Lin(≽)`. So:
- **(i)** a single closed-form `f ∈ 𝓕(≽)` realizes *one* element of `Lin(≽)` — it satisfies the
  properties but makes an arbitrary commitment on each incomparable pair (which extension to pick).
- **(ii)** the average rank `r̄(x) = 𝔼_{⪰ ∈ Lin(≽)}[rank_⪰(x)]` is the **barycenter of the rank vectors**
  `{rank_⪰ : ⪰ ∈ Lin(≽)}` under the **uniform distribution on `Lin(≽)`**. It is not "the centroid of the
  admissible class", which is a set of orderings and has no centroid. It agrees with every `f ∈ 𝓕(≽)` on
  the comparable pairs (by admissibility, Lemma 1) and, on each incomparable pair `{a,b}`, it **settles
  the pair by mean rank** rather than by privileging either object. It commits to no tie-breaking rule
  and has no free parameter.

**What the average rank does not return (corrected).** It does **not** return the fraction of extensions
that rank `a` above `b`. That quantity is the **mutual-rank-probability matrix**, a strictly finer
summary, and it is the one that forces no side on any pair. The average rank is the scalar reduction of
that matrix, so it does order every incomparable pair, just without a free tie-breaking parameter.
Report the matrix when the ambiguity on a pair is itself the object of interest.

Hence **(ii) dominates (i)** in the minimal-commitment sense *whenever `Lin(≽)` can be explored* — which,
since **counting** the linear extensions of a finite poset is **#P-complete** (Brightwell & Winkler 1991)
and computing **exact average ranks** is consequently **#P-hard**, means *approximately sampled* (Bubley &
Dyer 1999). The two statements are not interchangeable. #P-completeness is a property of the counting
problem, and #P-hardness is what it implies for the average ranks, so do not predicate #P-completeness of
"exact enumeration" of average ranks. This is a conditional dominance: the capability to sample the
admissible class is exactly what converts "pick one arbitrary formula" into "report the consensus of all
of them." Two honest limits: (ii) is canonical *relative to* `≽`, that is, relative to the chosen base set
`𝓜`, and average-rank is the scalar summary of the fuller mutual-rank-probability matrix. (i) remains
preferable for per-object scoring in isolation or when an interpretable closed form is needed.

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

So Gini and entropy are the best single proxies for the weight-free consensus aggregate, and **ratio
commits most beyond consensus**, a quantitative version of "ratio is the outlier."

**⚠ Candidate-measure issue found (flag for the kinase paper).** Benchmarking the repo's *exact*
`candidate_measure.py` softplus-entropy candidate on Klaeger, it is **anti-correlated** with the
distribution family it is meant to join: Spearman(candidate, entropy) = −0.87, vs Gini −0.87, S-score
−0.94, consensus-of-4 −0.86, robust across floors 5.0→0.0. Cause: 93.6% of Klaeger entries sit at the
detection floor and softplus never zeroes them (pedestal T·log2 ≈ 0.69 per kinase), so a genuinely
selective compound (one spike among 342 pedestals) looks maximally uniform → high entropy → the
candidate calls it *non*-selective (schematic worked example: a single-target profile gives candidate-H
≈8.4, near the 343-kinase maximum, while a many-target profile gives a *lower* H — so the promiscuous
compound is ranked as more selective, and the inversion direction is confirmed on data by candidate_bench.py,
−0.87 vs entropy, and the un-gated +0.94 vs n_active). Not a sign typo — the pedestal dominates. **Unaffected:**
the candidate's panel-size convergence (p*≈110–170) and D4 monotonicity, which concern self-consistency,
not agreement direction.

**Deeper finding (a tradeoff, not a one-line fix).** Attempting the obvious fix — subtract the floor
pedestal so inactive kinases contribute ≈0 — *does* restore correct orientation (Spearman with entropy
+1.00, S-score +0.94, canonical +0.996, with p*=110 and D4 100%) **but breaks D3**: its baseline-robustness
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
per-compound *min* is the weak +0.51 version, and a robust quantile is what works.) **Two clean resolutions
therefore exist:** (i) this minimally-modified closed-form candidate (fixed floor → per-compound-quantile
anchor), and (ii) the **baseline-free average-rank canonical extension** (§5). Both are correct and
weight-free and baseline-free. The average-rank is more principled (no anchor choice at all), the anchored
closed-form is a smaller edit to the existing kinase candidate. Scripts: `analysis/candidate_bench.py`,
`candidate_fix.py`, `candidate_fix2.py`, `candidate_fix3.py`. Not edited in the kinase repo — flagged to
the author to verify.

---

## 6. The empirical work, not the theorems, is the contribution (main text §4–§5)

The theorems above are scaffolding. The controlled findings are the contribution, and each is stated
with its honest scope.

1. **Disagreement is a near-tie phenomenon.** Kinase pairwise discordance is driven by *separation* on
   the selectivity axis: equal-n_active pairs are 76% discordant, well-separated pairs (Δn_active ≥ 16)
   only 19%. Axis **position adds nothing** beyond separation (added R² ≈ 0.002, middle-position
   coefficient −0.11). → Objects close in κ get ordered inconsistently, and where near-ties fall is a
   density property, **not** a concept-level "consequential middle."
2. **No "consequential middle" concentration (P3 retired).** A naive partial Spearman of disagreement
   on middle-position *is* significant (0.315, p=4e-6, controlling n_active/max_act), which is exactly
   the trap. But that control does not remove rank-boundary compression. The proper boundary-compression
   **independence null** does, and against it the apparent interior peak vanishes: observed slope
   (0.11–0.27) lies inside the null band (up to 0.42), at or below the null mean for the consensus
   coordinate (0.265 vs 0.299). So the interior peak is
   a mechanical artifact, not concept structure. In smoking, the former-smoker divergence (per-cohort
   classification agreement 78% and 66%, with the source repo's pooled figure at 71%) is a
   **binary-threshold** effect, robust to the cutoff but still a boundary phenomenon. A threshold-free
   metric does **not** put the peak at the intermediate group, and is itself density-confounded, so no
   directional concentration law is claimed.
3. **Reliability-gate / detection-floor instability is real and robust.** Kinase zero-active compounds
   carry the large instability (disagreement ~40–92 vs ~21–25, roughly 1.9–3.7× the active level),
   stable across every pKd active cutoff (5.5–7.0). The main text reports the same effect as a rank
   standard deviation of 72.9 for zero-active compounds versus 31.3 for active ones. These are exactly
   what the gate (G1) excludes, and excluding them is necessary, not optional.
4. **Reproducibility across cohorts (G4) is directionally supported.** In smoking the classification-
   agreement ordering (never > current > former) reproduces in both cohorts separately (GSE50660
   .95/.91/.78, GSE42861 .96/.76/.66). Magnitudes differ, and one cohort is a rheumatoid-arthritis
   study (disease confound), which is stated as a limitation.

**Relation to prior framings (cite, delimit):** the situation is *convergent validity* (Campbell &
Fiske 1959), multiple operationalizations of one construct that disagree. Psychometrics treats mid-
range divergence as **method variance to remove**. Our contribution is to separate divergence that is
**constitutive** (concept-level incomparability, surviving reproducibility) from divergence that is
**artifactual** (near-ties, detection-floor noise, threshold boundaries). The "agree-at-extremes /
diverge-in-middle" pattern alone is also the textbook Bland–Altman proportional-bias/range effect. What
is *not* generic is the reproducibility test that promotes a disagreement to a property of κ.

---

## 7. Proved vs assumed vs open

- **Classical/definitional (used, not claimed):** Lemma 1, Lemma 2, the non-uniqueness Proposition
  (Szpilrajn 1930 for non-emptiness, standard order theory for the bijection and the "unique iff total"
  criterion), the average-rank canonical extension (Patil & Taillie 2004, De Loof et al. 2008, 2011),
  and the sampling results it depends on (Brightwell & Winkler 1991 for #P-completeness of counting
  linear extensions, Bubley & Dyer 1999 for near-uniform sampling).
- **Assumed:** `≽` exists as a property of κ and reasonable measures respect it (operationalized as the
  reproducible consensus order). One substantive assumption, argued rather than proved.
- **Empirical (this work's contribution):** §6.1–§6.4, with the stated scope and limitations.
- **Deferred, not open problems:** benchmarking the softplus-entropy candidate against the canonical
  extension (needs the repo implementation), and continuous-exposure concentration analysis in smoking
  (the public cohorts lack a continuous exposure variable, a data limitation, flagged).
