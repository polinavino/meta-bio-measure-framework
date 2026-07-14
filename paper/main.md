# Which definition? A protocol for comparing competing measures of a latent biological concept

*Polina Vinogradova*

**Draft — main text (methods / Perspective).** The paper proposes a *technique* bioinformaticians can
run when several formally distinct measures each claim to quantify one latent concept and disagree; it
justifies the technique with established measurement theory (Appendix A) and demonstrates it across four
domains. It does **not** claim new mathematics. Register: comp-bio audience; formal statements in
Appendix A (`formal-spine.md`). All numbers are sourced in `../synthesis/repo-findings.md` and the
controlled re-analyses in this work.

---

## Abstract

Many biological concepts are quantified by several competing formulas at once — kinase inhibitor
*selectivity* (S-score, selectivity entropy, Gini, target-ratio), *biological age* (position, deviation,
and rate epigenetic clocks), environmental *exposure* or a *disease signature* (competing DNA-methylation
panels). On the same data these formulas rank the same objects differently, and the field's usual
response is to benchmark them and pick a winner. We argue this is mis-specified and give a better
procedure. Treating each measure as an order-preserving representation of the concept's comparative
structure — the standard view of representational measurement theory, and of the partial-order theory of
multi-indicator systems (Patil & Taillie 2004) — makes the disagreement interpretable: measures must
agree wherever the concept's order is determined and can differ only where objects are genuinely
incomparable, so no measure is canonical once incomparabilities exist (Szpilrajn). The situation is a
*convergent-validity* problem (Campbell & Fiske 1959); our contribution is a protocol that separates
disagreement that is **constitutive** (a property of the concept, surviving cross-cohort reproducibility)
from disagreement that is **artifactual** (near-ties, detection-floor noise, decision-threshold
boundaries). Applying the protocol with proper controls overturns a natural but incorrect belief — that
disagreement concentrates in a "consequential middle." In the kinase domain, inter-measure discordance is
driven by how close objects are on the selectivity axis (near-ties), not by where they sit (position adds
no predictive value beyond separation); the apparent middle-peak is a rank-boundary artifact. What is
robust and reproducible is (i) large instability at the detection floor, which the reliability gate
correctly excludes, and (ii) the ordering of cohort-level agreement. When a single summary is
unavoidable, we recommend the weight-free *average-rank-over-linear-extensions* canonical measure (De
Loof, De Baets & De Meyer), which commits no further than the consensus order forces; in kinase it is
best approximated by Gini and selectivity entropy and worst by the ratio measure, quantifying the sense
in which the ratio is an outlier. We provide a domain-independent desiderata checklist (reliability gate,
monotonicity, ordinal stability, cross-cohort reproducibility, intervention consistency) of which the
published per-domain lists are specializations, and worked demonstrations in kinase pharmacology,
epigenetic ageing, serotonin pharmacology, and methylation biomarkers.

---

## 1. Introduction

For many biological concepts, measurement is definitional before it is empirical: to say compound A is
more selective than B, or that an intervention slowed someone's ageing, one must first choose a *number*
to stand for the concept — and for most interesting concepts several defensible numbers are in routine
use. These are not implementations of one agreed formula; they are different formulas, and on the same
data they rank the same objects differently. Faced with this, the literature typically benchmarks the
candidates and anoints a winner.

We argue the benchmark question — *which measure is correct?* — is usually the wrong one, and we give a
procedure that asks better questions. The starting point is old and well established: a measure is an
order-preserving map from a concept's comparative structure into the numbers (representational
measurement theory: Krantz, Luce, Suppes & Tversky), and when several indicators target one concept the
admissible aggregates are the *linear extensions* of the partial order they jointly induce (Patil &
Taillie 2004; Brüggemann & Patil 2011). We do not extend this theory; we bring it to bioinformatics as a
usable protocol and show what it reveals — and, importantly, what it *dissolves*.

Two consequences of the view are immediate and classical. First, any two admissible measures must agree
wherever the concept's order is determined and can disagree only on pairs the concept leaves
**incomparable** (Appendix A, Lemma 1). Second, whenever incomparabilities exist, the desiderata one
would impose on a "good" measure pin down a whole *class* of admissible measures — the linear extensions
— and cannot single out a unique one (Szpilrajn; Appendix A §4). So "no measure satisfies every
requirement, and none is canonical" is not a failure of the field's effort; it is a theorem. The right
output is therefore not a winner but a *map*: where does the concept decide, where is it silent, and —
if a single number is unavoidable — what is the principled minimal-commitment summary.

**What is genuinely contributed here** is not the formalism (it is textbook and we cite it) but:
1. **A protocol** (§2) that operationalizes this view for a bioinformatician with several measures and
   some cohorts, including the controls needed to avoid over-interpretation.
2. **A separation of constitutive from artifactual disagreement.** The situation is a convergent-validity
   problem (Campbell & Fiske 1959); psychometrics treats mid-range divergence as *method variance to be
   removed*. We show part of it is *constitutive* — the concept genuinely declines to order some objects —
   and give the reproducibility test that tells the two apart (§2 step 3, §4).
3. **Controlled empirical findings** (§5) that correct a natural error: disagreement does **not**
   concentrate in a "consequential middle." With boundary and separation controls, kinase discordance is
   a near-tie phenomenon whose location carries no extra information; the robust, reproducible signals are
   detection-floor instability (which the reliability gate excludes) and cohort-level agreement ordering.
4. **A recommended canonical measure** for when a single summary is required: the weight-free average-rank
   over linear extensions (§2 step 7, §5), new to bioinformatics practice though standard in the
   multi-indicator literature.

§2 states the protocol; §3 places it against prior work and states what each step borrows; §4 works
through four domains; §5 reports what the demonstrations establish, with honest scope; §6 discusses
reporting practice and limits.

---

## 2. The protocol

Input: a latent concept κ; a set Ω of objects; several candidate measures f₁…f_k each claiming to
quantify κ; ideally ≥2 independent cohorts. Each step names the established result that justifies it
(Appendix A).

1. **Declare type and scale.** For each measure, state what comparative structure it represents (what it
   would mean for one object to exceed another) and its scale of measurement. Measures of different
   declared types are not rival estimates of one number and must not be benchmarked against each other.
   *(RMT; Appendix A Def 1–2.)*
2. **Compute the consensus order and the disagreement set.** Form the partial order in which a ≽ b iff
   *every* measure ranks a ≥ b. The pairs on which the measures disagree are the candidate
   incomparabilities. *(Patil & Taillie; Appendix A Lemma 1.)*
3. **Test reproducibility (the constitutive/artifactual split).** Recompute the consensus order and the
   disagreement set in each cohort separately. Disagreement that reproduces across independent cohorts is
   evidence of concept-level incomparability (constitutive); disagreement that does not is sampling noise
   or method variance (artifactual). *(This is the step that distinguishes us from a convergent-validity
   audit; Appendix A §1.1.)*
4. **Locate disagreement — with controls.** Before claiming disagreement concentrates anywhere, control
   for (a) object separation on the concept axis (near-ties trivially disagree), (b) local object density,
   (c) distance to any decision threshold, and (d) rank-boundary compression (via an independence null).
   Report only what survives these controls. *(§5 shows why: uncontrolled, one "finds" a spurious middle.)*
5. **Cluster measures into families and score them against the desiderata.** Group measures by the
   ordering they induce (rank agreement); within a family they agree up to a monotone reparameterization
   (ordinal equivalence — *not* a shared interval/ratio scale; Appendix A Lemma 2). Then fill the
   desiderata checklist (§3) — reliability gate, monotonicity, ordinal stability under nuisance,
   cross-cohort reproducibility, intervention consistency — recording which measures satisfy which.
6. **Report a map, not a winner.** Because no measure is canonical when incomparabilities exist
   (Appendix A §4), report the comparability skeleton (where the concept decides), the reproducible
   incomparable set (where it does not), the families, and the desiderata table — rather than a single
   "best" measure.
7. **If a single summary is unavoidable, use the minimal-commitment canonical aggregate.** Compute the
   average rank of each object over the linear extensions of the consensus poset (interval-order variant
   if scores carry uncertainty). It is weight-free and commits no further than the consensus order forces.
   Benchmark existing measures against it to see which is the best single proxy and which commit beyond
   the consensus. *(De Loof, De Baets & De Meyer; Appendix A §5.)*

The protocol's value is concentrated in steps 3–4: without them, naive disagreement metrics manufacture
structure that proper controls dissolve (§5).

---

## 3. Foundations and relation to prior work

We build on, and delimit against, established work (verified references and positioning in
`related-work.md`). The formalism is entirely prior art; we cite it as foundation.

- **Representational measurement theory** (Krantz, Luce, Suppes & Tversky, *Foundations of Measurement*,
  1971/1989/1990). Measurement as homomorphism; scale types as admissible-transformation groups. We use
  its representation view and its ordinal-level caution (we do *not* claim interval/ratio structure we
  cannot establish). Its clearest prior application in biology is **Houle et al. (2011), "Measurement and
  Meaning in Biology"** — which stays at the level of single measurements and scale types, not competing
  measures of one concept as extensions of a shared partial order.
- **Partial-order theory of multi-indicator systems** (Patil & Taillie 2004; Brüggemann & Patil 2011;
  in chemistry back to **Halfon & Reggiani 1986**; the social-indicator analogue **Fattore 2016**).
  Multiple indicators of one concept → a dominance partial order → linear extensions → no unique ranking.
  This is our formal spine; our contribution is bringing it to biological measurement and pairing it with
  a reproducibility-based constitutive/artifactual test.
- **The average-rank canonical extension** (De Loof, De Baets & De Meyer 2008/2011; Lerche et al. 2003),
  estimated via near-uniform linear-extension sampling (**Bubley & Dyer 1999** — the sampler we use;
  counting is #P-complete, **Brightwell & Winkler 1991**). The weight-free minimal-commitment summary; we
  recommend it as the protocol's default (step 7) and demonstrate it (§5). An earlier version of this work
  posed finding such a measure as an open problem — it is already solved here; we adopt their answer.
- **Convergent validity / MTMM** (Campbell & Fiske 1959). "Multiple operationalizations of one construct
  disagree" is a 60-year-old named problem, with method variance as the standard explanation for
  mid-range divergence. Our departure: divergence is not always a defect to scrub — part is *constitutive*
  (the concept is silent), and reproducibility separates the two.
- **Social choice / rank aggregation** (Arrow; Singleton & Booth 2022; Dwork et al. 2001 on NP-hard
  Kemeny consensus; Kolde et al. 2012 robust rank aggregation in bioinformatics). No-canonical-aggregator
  results parallel our non-uniqueness; the aggregation literature *forces* a consensus, whereas we
  *preserve* incomparability and only summarize when required.
- **Domain precedents that already show the phenomenon:** competing kinase-selectivity metrics disagree
  (Uitdehaag & Zaman 2011; Bosc et al. 2017); gene-expression breast-cancer predictors share few genes
  yet give concordant patient predictions (**Fan et al. 2006**) — an empirical instance of "agree where
  the order is determined"; epigenetic clocks disagree and differ in reliability (Belsky et al. 2022).

We also note the raw pattern "agree at extremes, diverge in mid-range" is, by itself, the textbook
Bland–Altman proportional-bias/range effect (Bland & Altman 1995); §5 shows the interesting content is
*not* that pattern but the reproducibility test and the controls that determine what is real.

---

## 4. Worked demonstrations

Each subsection applies the protocol. **Provenance:** the kinase near-tie/gate/canonical results (§4.1),
the smoking per-cohort agreement (§4.4), and the cross-domain transfer (§5) are computed here with
controls (`analysis/`). The remaining figures — clock affinities/cosine/cell-type (§4.2), FASD AUCs
(§4.4), serotonin S-score–Gini −0.682 (§4.3) — originate in the source-repo analyses and were
**reproduced/confirmed locally this session** (`analysis/repro_borrowed.py`, `analysis/fasd_auc.R`):
serotonin −0.682 exact; position-clock R² 0.991–0.995, rate–position signed-log R² 0.93–0.95,
cell-type max |r| 0.372, all confirmed; FASD AUCs buccal 0.68–0.79, van der Laan 0.935/0.959 confirmed.
One correction surfaced: the clock coefficient-cosine range is −0.17 to 0.11 (not "0.03–0.11"; the
"near-orthogonal" conclusion is unchanged).

### 4.1 Kinase selectivity (full worked example)
Concept: inhibitor selectivity. Measures: S-score, selectivity entropy, Gini, target-ratio, across four
datasets/three assay technologies (Davis 68×433, Klaeger 222×343, Anastassiadis 178×300, Metz 704×172).
- **Families (step 5):** two families by induced order — distribution {S-score, entropy, Gini}
  (within-r 0.74–0.99) and the target-ratio (cross-family r 0.14–0.62). Ordinal classes only.
- **Reliability gate (step 5):** zero-active compounds carry large instability (disagreement ~40–92 vs
  ~21–25 for active; ≥~2× the active level at every pKd active cutoff 5.5–7.0). These are exactly what the gate
  excludes; the "extreme instability" is a domain-of-definition issue, not concept structure.
- **Locate disagreement, controlled (step 4):** pairwise discordance is governed by separation on the
  selectivity axis — equal-n_active pairs 76% discordant, well-separated pairs (Δ≥16) 19% — and axis
  **position adds essentially nothing** (added R² ≈ 0.002; middle coefficient −0.11). The apparent
  interior peak is inside the rank-boundary independence null. **No "consequential middle."**
- **Canonical aggregate (step 7):** consensus poset over 206 gated-in compounds is 37.2% incomparable;
  the average-rank extension is best proxied by Gini (+0.955) and entropy (+0.950), worst by ratio
  (+0.796) — a quantitative "ratio is the outlier."
- **Desiderata:** no existing measure satisfies all (D1–D4); consistent with non-uniqueness.

### 4.2 Epigenetic clocks (families and a scale caution)
Concept: biological age. Five clocks, three *declared* types (position, deviation, rate) on 1,385 blood
samples (GSE40279 n=656, GSE87571 n=729). Position clocks are affine-related (R² 0.991–0.995) — reported
**descriptively**, with the caution that both are regressions on chronological age, so affinity is
expected and is *not* evidence of a shared interval-scale latent structure (Appendix A §3). Rate vs
position relate non-affinely (signed-log R² 0.93–0.95 vs linear 0.75–0.79); reported as a between-measure
observation, not a scale-type proof. Cell-type reproducibility fails for all clocks (max |r| with cell
fractions up to 0.372, Hannum), and the age-trained clocks are near-orthogonal in coefficient space
(|cosine| ≤ 0.17) — a G1/G4/D0 flag. *(These figures confirmed locally from the clock repo's computed
outputs; `analysis/repro_borrowed.py`.)* This domain best illustrates step 1 (declaring type) and the
ordinal-only caution in step 5.

### 4.3 Serotonin selectivity (families replication — and a convention caveat)
The kinase measures applied to 5-HT receptor pharmacology (ChEMBL, 13,584 compounds × 13 receptors; 297
with the required overlap). A family structure recurs, but **which** measure is the outlier differs
(entropy here, not ratio), and the reported anti-correlation (S-score vs Gini −0.68) is entangled with
inconsistent higher-vs-lower-is-selective conventions across sub-analyses. Demonstrates step 1's
importance (declare orientation/scale first) and the honest limit: family *membership* is not
domain-invariant; only the existence of families and an outlier is.

### 4.4 Methylation biomarkers (a domain-gate demonstration, and a threshold artifact)
Smoking arm (within-tissue): four signatures; two cohorts (GSE50660 n=464, GSE42861 n=689). The
classification-agreement ordering never > current > former reproduces per cohort (GSE50660 .95/.91/.78;
GSE42861 .96/.76/.66; step 3), but the former-smoker "divergence" is a **binary-threshold boundary**
effect: never-smokers agree most at every cutoff, but *which* non-never group agrees least shifts with
the (arbitrary) threshold — former at the 2SD headline, current at strict cutoffs in the RA cohort — a
decision-boundary phenomenon, not concept structure (step 4).
FASD arm (cross-tissue): the dominant failure is a **reliability/domain-gate (G1) violation** — buccal
signatures scored on blood transfer poorly (AUC 0.68–0.79); the strong blood-native episignature
(0.93–0.96) is **non-independent** (its derivation cohort overlaps the evaluation samples — confirmed by
the depositing author) and is an upper bound only. Shows the protocol classifying *which* desideratum
fails, and step 3's independence caveat.

---

## 5. What the demonstrations establish

With the protocol's controls applied, the honest, reproducible findings are:

1. **Disagreement is a near-tie phenomenon (not a consequential middle).** Discordance tracks how close
   objects are on the concept axis, not where they sit; controlling separation, position is
   uninformative (kinase §4.1). The intuitive "measures split on the interesting middle cases" is, under
   control, a rank-boundary and near-tie artifact. We retire it as a general claim.
2. **Detection-floor instability is real and gate-necessary.** The largest disagreement occurs where the
   signal is weakest (zero-active compounds; unexposed samples) — precisely the region the reliability
   gate (G1) excludes. Excluding it is necessary, not optional, and resolves apparent "extreme
   instability."
3. **Reproducibility is what makes disagreement meaningful.** The one leg that promotes "these measures
   disagree here" to "the concept is silent here" is cross-cohort reproducibility (step 3). Its evidence
   is currently modest (two cohorts per domain; one with a disease confound) and we state it as such.
4. **A single formula is valid, but the class centroid is better when the class can be explored.** The
   required properties do not pick a unique measure — they characterize a *class*, precisely the linear
   extensions of the concept's comparable order (Appendix A). Two consequences:
   - **(i) any closed-form measure satisfying the properties is a legitimate representative** — usable,
     and per-object deployable. But it makes an *unavoidable arbitrary commitment on every incomparable
     pair*: it silently picks one linear extension, i.e. one way of breaking exactly the ties the concept
     leaves open. The kinase candidate makes this concrete — pinning down one "well-formed" formula
     required arbitrary choices (smoothing width, a fixed floor vs a per-compound anchor, a hard gate)
     and hit a genuine orientation-vs-baseline-robustness tradeoff (§4.1; Appendix A §5). That difficulty
     is not a failure; it is the non-uniqueness theorem made tangible.
   - **(ii) when the admissible class can be explored — its linear extensions (approximately) sampled —
     a minimal-commitment summary is available and preferable: the average rank over the whole class**
     (the average-rank canonical extension). It is parameter-free, privileges no tie-breaking, and on
     each undecided pair reports how the admissible measures divide rather than forcing a side. In this
     precise sense (ii) *dominates* (i): it is the centroid of exactly the space the properties define.
     We compute it here (kinase: 37% of pairs incomparable; average rank estimated by near-uniform
     linear-extension sampling; Gini/entropy proxy it best at ρ≈0.95, the ratio worst at 0.80 — a
     quantitative "ratio commits most beyond the consensus").
   Caveats: (ii) is canonical *relative to* the comparable skeleton `≽` (the reproducible consensus of the
   reasonable base measures), not in a vacuum — the choice of base set is real, pushed up one level, not
   eliminated. "Explore the class" means *sample* it (exact counting is #P-complete; we use Bubley–Dyer),
   and (i) remains preferable when an object must be scored in isolation or an interpretable closed form
   is required. So: use a single formula if you must, but if you can sample the class, report its centroid.
5. **The near-tie law is domain-general (the generality is earned, qualitatively).** Overlaying the
   separation→discordance decile curves across domains and measuring pairwise agreement
   (`analysis/transfer.py`; kinase, serotonin, clocks, smoking — different concepts, measures, object
   types; no curve is fit), the curve is strongly monotone-decreasing in every domain (Spearman −0.99 to −1.00) and
   transfers across kinase/clocks/smoking (decile-curve RMSE 0.10–0.14). Serotonin is the exception
   (RMSE 0.40–0.50): its measures stay discordant even for well-separated compounds, consistent with its
   tiny 13-receptor panel and mixed conventions (§4.3). Honest reading: the *qualitative* law (near-ties
   disagree, well-separated objects agree) is domain-general and is the transferable content; the
   *quantitative* decay rate is domain-specific. This supports methodological — not numerical — generality.

These are less sweeping than "the same five-part phenomenon recurs in four domains," and deliberately so:
they are what survives controls. The protocol's contribution is precisely the machinery that separates
these from the artifacts.

---

## 6. Discussion

- **Reporting practice.** Declare scale type and orientation (step 1); report the reproducible
  incomparable set and the desiderata table rather than a single ranking (step 6); validate measures in
  their native domain/tissue/panel (G1, §4.4); when a summary is needed, report the average-rank
  canonical aggregate and each measure's distance from it (step 7).
- **The conceptual claim.** Some disagreement among biological measures is not a defect to be
  benchmarked away but a faithful reflection of a concept that does not totally order its objects
  (contra the method-variance reading of convergent validity). The protocol's job is to find the part
  that is real and reproducible, and to stop the field from adjudicating the part that is not there.
- **Limitations.** One substantive assumption (a concept-level order exists and reasonable measures
  respect it). Reproducibility evidence is thin (few cohorts). Family membership is not domain-invariant.
  The public smoking cohorts lack a continuous exposure variable, so the sharpest location analysis is
  available only in kinase. The formalism is entirely borrowed; the contribution is protocol + controls +
  cross-domain demonstration.
- **Scope.** A companion axis — that *evaluations* (not definitions) overstate performance via scaffold
  memorization in QSAR/ADMET — is cited as a second face of the same self-deception about measurement,
  not folded in.

---

## Appendix A — formal justification of the protocol
See `formal-spine.md`: definitions, Lemmas 1–2 (analytic), the Szpilrajn non-uniqueness proposition, the
average-rank canonical extension, and §6 there (the controlled empirical findings with scope). To be
inlined as the technical appendix.
