# A measurement-theoretic protocol for competing quantifications of latent biological concepts

*Polina Vinogradova*
Independent Researcher.
Correspondence: polina.vino@gmail.com. ORCID: 0000-0003-3271-3841.

> **Target:** *Journal of Theoretical Biology* (hybrid, traditional/subscription route, no APC).
> Elsevier "Your Paper Your Way" allows flexible formatting at initial submission. A formatting pass to
> JTB style (reference format, figure files) happens after the scientific draft is complete. This file
> is the canonical submission draft. `main.md` and `formal-spine.md` are the working notes it draws on.
> **All quantitative claims are copied from the stored analysis outputs** (`tp53/analysis/outputs/`,
> `inflammation/analysis/outputs/`, `analysis/`), which are the single source of truth.

---

## Highlights
- Competing biological measures of one concept are linear extensions of a shared partial order.
- So disagreement is expected, and the right output is a comparability map, not a single winner.
- A protocol separates concept-level (constitutive) from artefactual disagreement.
- Demonstrated across six domains: selectivity, ageing, exposure, TP53, inflammation, sepsis.
- Cross-cohort reproducibility is the test that promotes disagreement to a concept property.

## Abstract
Many biological concepts are quantified by several competing formulas at once. This includes 
kinase-inhibitor
selectivity, biological age, the functional impact of a TP53 variant, and inflammatory burden, each 
of which have
multiple defensible measures resulting in distinct rankings. The field's usual
response is to benchmark them and choose a winner. We argue this is mis-specified, and, instead,
suggest treating each measure
as an order-preserving representation of the concept's comparative structure as data.
Admissible measures must agree wherever the concept's order is determined, and can differ
only on objects the concept leaves incomparable, so once incomparabilities exist, no measure is canonical.
The right output is a map of where the concept decides and where it is silent, not a winner. Our
contribution is a reusable protocol that operationalizes this view, with the controls that separate
constitutive disagreement (a concept property that survives cross-cohort reproducibility) from artefactual
disagreement (near-ties, detection-floor noise, decision-threshold boundaries). When a single summary is
unavoidable, we recommend the weight-free average rank over linear extensions, which commits no further
than the consensus order forces. We demonstrate the protocol on public data across six domains, reporting
honestly where it is strong, where a domain differs, and where reproducibility is weak.

**Keywords:** measurement theory, partially ordered sets, convergent validity, rank aggregation,
biomarkers, variant effect prediction

---

## 1. Introduction

For many biological concepts, measurement is definitional before it is empirical. To say that compound A
is more selective than B, that an intervention slowed someone's ageing, or that a genetic variant is more
damaging than another, one must first choose a number to stand for the concept. For most interesting
concepts several defensible numbers are in routine use. These are not implementations of one agreed
formula. They are different formulas, and on the same data they rank the same objects differently. Faced
with this, the literature typically benchmarks the candidates against some reference and anoints a winner.

We argue that the benchmark question, *which measure is correct?*, is usually the wrong one, and that a
better question follows from an old and well-established view of what a measurement is. In representational
measurement theory (Krantz et al., 1971, Suppes et al., 1989, Luce et al., 1990), a measure is an
order-preserving map from a concept's comparative structure into the real numbers. When several indicators
target one concept, the aggregates consistent with all of them are exactly the *linear extensions* of the
partial order the indicators jointly induce (Patil and Taillie, 2004, Brüggemann and Patil, 2011). We do not
extend this theory. We bring it to biological measurement as a usable protocol, and show what it reveals,
and what it dissolves.

Two consequences of the view are immediate and classical. First, any two admissible measures must agree
wherever the concept's order is determined, and can disagree only on pairs the concept leaves
**incomparable** (Section 2, Lemma 1). Second, whenever incomparabilities exist, the properties one would
demand of a "good" measure pin down a whole *class* of admissible measures, the linear extensions, and
cannot single out a unique one (Szpilrajn, 1930, Section 2). So "no measure satisfies every requirement,
and none is canonical" is not a failure of the field's effort. It is a theorem. The right output is
therefore not a winner but a *map*. Where does the concept decide, where is it silent, and, if a single
number is unavoidable, what is the principled minimal-commitment summary?

This reframing matters because the alternative is actively misleading. Benchmarking competing measures to
a winner treats concept-level ambiguity as estimation error to be resolved, when in fact it may be a
faithful reflection of a concept that does not totally order its objects. A TP53 variant with partial,
promoter-selective loss of function, a former smoker between clearly-exposed and clearly-unexposed, a
compound of intermediate selectivity: these are objects on which the concept itself may decline to
adjudicate. Reporting a single ranking hides exactly the cases where honesty about ambiguity matters most.

**What is genuinely contributed here** is not the formalism, which is textbook and which we cite as
foundation, but:

1. **A protocol** (Section 3) that operationalizes this view for a working analyst who has several
   measures and some cohorts, including the controls needed to avoid over-interpretation.
2. **A separation of constitutive from artefactual disagreement.** The situation is a *convergent-validity*
   problem (Campbell and Fiske, 1959). Psychometrics treats mid-range divergence as method variance to be
   removed. We show that part of it is *constitutive*, meaning the concept genuinely declines to order some
   objects, and we give the reproducibility test that tells the two apart (Section 3, step 3).
3. **A cross-domain demonstration** on public data (Section 4) spanning kinase pharmacology, epigenetic
   ageing, serotonin pharmacology, methylation biomarkers, TP53 variant effect, and inflammatory burden.
   We report where the protocol is strong, where a domain behaves differently, and where the evidence is
   weak, including two honest null and negative results.

The theoretical content, in the sense a reader of this journal will care about, is a statement about the
*structure of biological concepts*. Several of them do not impose a total order on their objects, and the
persistent disagreement among their measures is the observable signature of that partial structure. The
protocol is the apparatus that reads the signature.

Section 2 states the measurement-theoretic foundations and what each licenses. Section 3 gives the
protocol. Section 4 works through the six domains. Section 5 reports what the demonstrations establish,
with scope. Section 6 discusses reporting practice and limits. Appendix A collects the formal statements.

---

## 2. Foundations

We fix the object and record which established result justifies which step of the protocol. Everything in
this section is classical or definitional. The contribution is the protocol (Section 3) and the empirical
work (Sections 4 and 5), not the theory. Fuller statements and proofs are in Appendix A.

**The object.** Fix a latent concept κ (selectivity, biological age, variant impact, inflammatory burden)
and a finite set of objects Ω. Domain knowledge licenses *comparative judgments* a ≽ b ("a is at least as
κ as b"). The empirical structure is (Ω, ≽) with ≽ a preorder that is **not** assumed total. Write a ∥ b
when neither direction is licensed, meaning a and b are *incomparable*. A *measure* f : Ω → ℝ is
**admissible** if a ≽ b implies f(a) ≥ f(b), an order-preserving homomorphism, which is the representation
condition of measurement theory. Two measures belong to the same **family** if they induce the same total
preorder. A family is an equivalence class under strictly monotone reparameterization, that is, an
*ordinal* class (Lemma 2). This is exactly the multiple-indicators-to-partial-order-to-linear-extensions
setup of Patil and Taillie (2004), which we adopt.

**Disagreement lies in the incomparable set (Lemma 1).** If two admissible measures rank a pair in
opposite strict directions, that pair is incomparable under κ. Equivalently, admissible measures agree on
the strict order of every comparable pair, and all inter-measure reversals lie among the incomparabilities.
This justifies protocol step 2. The pairs on which reasonable measures disagree are exactly the candidate
incomparabilities.

**The order is estimated by consensus, and the estimate over-declares comparability.** The true ≽ is not
observed. Given the set 𝓜 of measures experts actually use, define the consensus order a ≽𝓜 b iff
f(a) ≥ f(b) for every f in 𝓜. Under the assumption that every reasonable measure is admissible, ≽ ⊆ ≽𝓜,
so the consensus can certify a pair as comparable when the measures agree by coincidence though κ leaves
it incomparable. Two guards make the estimate trustworthy. Cross-cohort reproducibility promotes a
disagreement that replicates across independent cohorts to a property of κ rather than of the sample
(protocol step 3). A crude external anchor for the concept should agree with the consensus at the
extremes, and it flags a collectively mis-oriented set of measures when it does not (step 1).

**No measure is canonical once incomparabilities exist.** On a finite (Ω, ≽), the admissible measures, up
to family, are in bijection with the linear extensions of ≽, and there is a unique one iff ≽ is total
(Szpilrajn, 1930). If any pair is incomparable there are at least two order-distinct admissible measures,
and no requirement that only entails order-preservation can single one out. This justifies step 6.
"Benchmark to crown the correct measure" is mis-specified, and the correct output is the comparability
skeleton plus the reproducible incomparable set. The same conclusion appears in social choice (Arrow, 1963)
and truth-discovery axiomatics (Singleton and Booth, 2022).

**The minimal-commitment summary.** When a single number is unavoidable, the average rank of each object
over the linear extensions of the consensus poset (Patil and Taillie, 2004, De Loof et al., 2008, 2011) is
the centroid of the admissible class. It agrees with every admissible measure on the comparable pairs, and
on each incomparable pair it reports the fraction of extensions ranking one object above the other rather
than forcing a side. It is weight-free and parameter-free. Exact enumeration is #P-complete (Brightwell
and Winkler 1991), so it is estimated by near-uniform sampling (Bubley and Dyer, 1999). This is the
principled resolution of the non-uniqueness above, and it justifies step 7. It is canonical *relative to*
the chosen base set 𝓜, so the choice of which measures to admit is real and only lifted one level, not
eliminated (Section 6).

**A scale-type caution.** Families are ordinal classes only (Lemma 2). Two admissible measures being
linearly related on one dataset does not establish a shared interval scale, which is expected, for
instance, when both are regressions on a shared external target. A non-affine monotone relation between
two measures, such as rate versus position clocks, is a descriptive observation rather than a scale-type
proof. The protocol reports relations between families descriptively and does not assign interval or ratio
structure it cannot establish.

---

## 3. The protocol

Input: a latent concept κ, a set Ω of objects, candidate measures f₁ … f_k each claiming to quantify κ,
and ideally at least two independent cohorts. Each step names the established result that justifies it
(Section 2, Appendix A).

1. **Declare type and scale, and verify orientation against an external anchor.** State, for each measure,
   what comparative structure it represents and its scale of measurement. Measures of different declared
   types are not rival estimates of one number. Then check that each measure moves with a crude,
   theory-light external proxy for the concept (number of active targets for selectivity, chronological age
   for clocks, clinical labels for variant effect, mortality for inflammation). A measure can pass every
   internal stability check yet point backwards, and only an external anchor catches a collectively
   mis-oriented set.
2. **Compute the consensus order and the disagreement set.** Form the partial order in which a ≽ b iff
   every measure ranks a ≥ b. The pairs on which the measures disagree are the candidate incomparabilities
   (Lemma 1).
3. **Test reproducibility, the constitutive/artefactual split.** Recompute the consensus order and the
   disagreement set in each cohort separately. Disagreement that reproduces across independent cohorts is
   evidence of concept-level incomparability (constitutive). Disagreement that does not is sampling noise
   or method variance (artefactual). This is the step that distinguishes the protocol from a
   convergent-validity audit.
4. **Locate disagreement, with controls.** Before claiming disagreement concentrates anywhere, control for
   object separation on the concept axis (near-ties trivially disagree), local density, distance to any
   decision threshold, and rank-boundary compression. Report only what survives.
5. **Cluster measures into families and score them against desiderata.** Group measures by the order they
   induce. Then fill a domain-independent checklist, recording which measures satisfy which: a
   reliability/domain gate, monotonicity in the concept against the anchor, ordinal stability under
   nuisance, cross-cohort reproducibility, and intervention consistency.
6. **Report a map, not a winner.** Because no measure is canonical when incomparabilities exist, report the
   comparability skeleton, the reproducible incomparable set, the families, and the desiderata table.
7. **If a single summary is unavoidable, use the minimal-commitment canonical aggregate.** Compute the
   average rank over the linear extensions of the consensus poset. Benchmark existing measures against it
   to see which is the best single proxy and which commit beyond the consensus.

The protocol's value is concentrated in steps 3 and 4. Without them, naive disagreement metrics
manufacture structure that proper controls dissolve (Section 5).

---

## 4. Worked demonstrations

We apply the protocol in six domains, all on public data. The domains differ in maturity and in how well
the protocol's assumptions hold, and we say so. Kinase selectivity is the running example and the most
developed. TP53 variant effect and inflammatory burden are the newest and most data-rich. Every number
below is copied from a stored analysis output in the project repository.

### 4.1 Kinase-inhibitor selectivity
Concept: how concentrated an inhibitor's activity is across the kinome. Measures: the S-score, selectivity
entropy, the Gini coefficient, and target-ratio, applied across four profiling datasets spanning three
assay technologies (Davis et al., 2011, 68 compounds, Klaeger et al., 2017, 222, Anastassiadis et al.,
2011, 178, Metz et al., 2011, 704, with prior comparisons of these metrics by Uitdehaag and Zaman, 2011,
and Bosc et al., 2017). The measures split into
two families that measure different things. The S-score, entropy, and Gini quantify concentration across
the panel and agree closely (within-family Spearman 0.74 to 0.99), while the target-ratio quantifies the
gap to one off-target and is the consistent outlier (cross-family 0.14 to 0.62). On the reliability-gated
Klaeger compounds the consensus poset leaves **37.2% of pairs incomparable**, and against the weight-free
average-rank aggregate the Gini (0.955) and entropy (0.950) are the best single proxies while the ratio
(0.796) commits most beyond the consensus. Disagreement is a near-tie phenomenon. Pairwise discordance is
driven by separation on the selectivity axis (equal-active-count pairs are 76% discordant, well-separated
pairs 19%), and axis position adds essentially nothing once separation is controlled (added R² about
0.002). Compounds with no activity above the detection floor carry the largest instability (rank standard
deviation about 74 versus 32 for active compounds), which is why the reliability gate is the first
desideratum. No existing measure satisfies the four desiderata we abstract from these failures, while a
gated smooth-hinge effective-target number does.

### 4.2 Epigenetic clocks (biological age)
Concept: biological age. Measures: DNA-methylation clocks of three declared *types*, which the protocol's
step 1 keeps separate rather than pooling. Position clocks report where a sample sits on the ageing
trajectory, deviation clocks report departure from it, and rate clocks report speed along it (1,385 blood
samples, GSE40279 and GSE87571). Rate clocks such as DunedinPACE (Belsky et al., 2022) are one declared
type, and independent comparisons confirm that clocks disagree and differ in predictive value
(Mavrommatis et al., 2025). Position clocks are affine-related (R² 0.991 to 0.995), and rate clocks
relate to position through a signed logarithm (R² 0.93 to 0.95), a Weber-Fechner-like relation that is a
descriptive observation, not a scale-type proof (Section 2). All clocks fail an ordinal-stability check
against cell-type composition (absolute correlation up to 0.372), and GrimAge, a two-stage composite, is
structurally near-orthogonal to the others (cosine of coefficient vectors between −0.17 and 0.11). The
external anchor here, chronological age, is the training target of the position clocks, so for them the
orientation check is confirmatory rather than independent. This is stated as a limitation.

The clock idea is not specific to blood or to methylation. On a rat retinal bulk RNA-seq aging series
(GSE314970, 85 samples over 6 to 27 months, Shavlakadze et al., 2026), five competing transcriptomic
aging measures disagree sharply on the aging axis: the SenMayo and Fridman senescence signatures (Saul et
al., 2022, Fridman and Tainsky, 2008), a Hallmark inflammatory inflammaging proxy, and a fitted
elastic-net clock. Against true age they range from Spearman 0.98 for the fitted clock down to 0.13 for
the Fridman senescence-up signature. The curated senescence signatures are weak individual age-trackers
while the fitted clock is near-perfect. They split into a senescence and inflammaging family plus a
distinct fitted-clock direction, the mean-rank consensus tracks age (0.755) and rises monotonically across
age tertiles, and discordance falls with separation as elsewhere. This is a single-cohort, rat,
single-sex extension with the clock trained in-distribution, so it is qualitative, on par with the
inflammaging clocks. It shows the clock family structure and the near-tie pattern recur in a new tissue
and a new measurement modality. Full analysis in the clock instance repository (`eye_aging/`).

### 4.3 Serotonin-receptor selectivity
Concept: 5-HT2A versus 5-HT2B selectivity. The kinase selectivity measures transfer directly (ChEMBL
n=297, Shulgin/PDSP n=36). The family structure recurs, but the *identity* of the outlier does not. Here
entropy is the outlier (its correlation with the others is near zero, while the S-score and Gini correlate
−0.682 after orientation), whereas in kinase the ratio was the outlier. Across domains the consensus order
still agrees with the external anchor (number of active targets, Spearman about 0.82). This domain is also
where cross-domain transfer is weakest. Discordance stays high even for well-separated compounds, so the
near-tie account explains less here than elsewhere, which we report rather than smooth over.

### 4.4 Methylation exposure biomarkers
Concept: environmental exposure or disease signature. Measures: competing smoking-exposure signatures
(AHRR, EpiSmoke, Joehanes, EpiTob) across two cohorts (GSE50660 n=464, GSE42861 n=689). The signatures
agree almost perfectly on never-smokers (96%) and well on current smokers (82%), and diverge on former
smokers (71% pooled, 78% and 66% per cohort). The divergence reproduces in both cohorts separately, which
makes it constitutive rather than sampling noise, but it is a binary-classification-threshold effect at a
decision boundary, not evidence of a concept-level "consequential middle" (Section 5). A cross-tissue arm
(fetal alcohol spectrum disorder) shows the reliability gate at work. Buccal-derived signatures scored on
blood transfer poorly (AUC 0.68 to 0.79), and the one strong blood-native signature (0.93 to 0.96) is
non-independent of its evaluation cohort (shared laboratory and tissue, confirmed by the depositing
author), so it is an upper bound only.

### 4.5 TP53 variant effect
Concept: the functional impact of a TP53 missense variant. Measures: eight Kato/Ishioka yeast
transactivation readouts (Kato et al., 2003), seven computational pathogenicity predictors (AlphaMissense,
Cheng et al., 2023, EVE, Frazer et al., 2021, ESM1b, Brandes et al., 2023, plus REVEL, CADD, PrimateAI and
BayesDel), and mammalian proliferation screens (Giacomelli et al., 2018, Kotler et al., 2018, via ProteinGym,
Notin et al., 2023), with ClinVar (Landrum et al., 2018) as the external clinical anchor and the ClinGen
TP53 specifications (Fortuno et al., 2021) for the field's functional rules (core panel 2,314 variants). Step 1 earns its place here. The eight
transactivation readouts and AlphaMissense are correctly oriented against ClinVar (AUC 0.96 to 0.99),
while all four proliferation screens read as mis-oriented (AUC 0.016 to 0.064) under the naive "higher
fitness is worse" convention, so they must be flipped before use. The measures form two families by
measurement basis, experimental (within-family mean 0.78) and computational (mean 0.76), agreeing less
across the two blocks (mean 0.50) than within either. The consensus poset is 56.6% incomparable, and the
field's own median-of-eight-promoters rule is the best proxy for the average-rank aggregate (Spearman
0.987) while all seven predictors sit well below it (0.46 to 0.60). The consensus discriminates ClinVar
labels at AUC 0.997. Disagreement is again separation-driven (0.93 for near-ties to 0.05 for well-separated
pairs), with one honest difference from kinase: a modest genuine middle residual survives control
(mid-axis coefficient +0.33 on the consensus axis, shrinking to +0.10 on an independent AlphaMissense
axis), corresponding to partial-loss-of-function variants. Two experimental *platforms*, yeast
transactivation and mammalian proliferation, agree strongly at the extremes (+0.63) and not at all in the
middle (−0.05), so the comparability skeleton is constitutive across platforms while the middle is
genuinely undecided. Finally, all seven predictors score 0.90 to 0.99 against curated ClinVar but only 0.79
to 0.89 against the independent functional truth on the variants of uncertain significance where prediction
is actually used, the variant-interpretation form of an evaluation-overstatement effect. A hotspot-codon
holdout showed no inflation for any predictor and is reported as an honest null.

### 4.6 Inflammatory burden
Concept: systemic inflammatory burden, quantified three ways of differing strength. **Clinical indices
(strong):** seven routine indices (NLR, PLR, MLR, SII, SIRI, CRP, CRP-to-albumin) on two independent
NHANES cycles, anchor all-cause mortality. Two families emerge by measurement basis, leukocyte-count
ratios and CRP-protein, nearly independent of each other (between-family mean 0.19), and CRP and its
albumin ratio are rank-identical (1.00), so one is redundant. The consensus discriminates mortality at the
extremes (AUC 0.69) but not the middle (0.51), and its inter-index correlation structure replicates across
the two cycles at Spearman **0.995**, the strongest reproducibility in this work. **Sepsis signatures
(strong disagreement, weak reproducibility):** five published blood signatures (the Hallmark inflammatory
and interferon-γ sets, Liberzon et al., 2015, the SRS7 and SRSq19 Sepsis Response Signatures, Davenport et
al., 2016, and the MARS8 endotype signature, Scicluna et al., 2017) on GSE65682 (802 samples), anchor
28-day mortality. The
Hallmark and SRS families are near-independent (cross-family mean +0.02), so two published sepsis-severity
signatures barely agree, but the detailed structure replicates only weakly across the discovery and
validation cohorts (Spearman 0.32), which we report as a limitation rather than a success. **Inflammaging
clocks (thin):** only two composite clocks are publicly reconstructable on shared samples, SImAge
(Kalyakulina et al., 2023) and a recomputed ipAGE (Yusipov et al., 2022). Two further inflammaging clocks,
iAge (Sayed et al., 2021) and IMM-AGE, are not reconstructable here, so this is a qualitative check, not a
full instance. Sepsis *endotype labels* (SRS
versus MARS) are partitions rather than scalar orders, a related but different object for which cluster
agreement, not average rank, is the right tool. Only the signature scores are analysed with the protocol.

## 5. What the demonstrations establish

With the protocol's controls applied, the honest and reproducible findings are the following. They are
less sweeping than "one phenomenon recurs identically in six domains," and deliberately so. They are what
survives the controls.

1. **Disagreement is a near-tie phenomenon, not a consequential middle.** In every domain where it can be
   tested, inter-measure discordance falls steeply as objects separate on the concept axis (kinase 76% to
   19%, TP53 0.93 to 0.05, indices 0.987 to 0.285). Once separation is controlled, axis position adds
   little or nothing. The intuitive claim that measures split on the "interesting middle" is, under
   control, largely a near-tie and rank-boundary artifact, and we retire it as a general law. TP53 is the
   honest exception, retaining a small genuine middle residual (partial-loss-of-function variants) that
   shrinks but does not vanish on an independent axis.

2. **Detection-floor and out-of-domain instability is real and gate-necessary.** The largest disagreement
   sits where the signal is weakest, at the kinase detection floor (rank standard deviation about 74
   versus 32) and in cross-tissue methylation transfer (AUC 0.68 to 0.79). This is exactly the region the
   reliability and domain gate excludes. Excluding it is necessary, not optional.

3. **Families recur, but the identity of the outlier does not.** Measures partition into a low-disagreement
   cluster plus at least one categorically different measure in kinase, serotonin, TP53, and the
   inflammatory indices. Which measure is the outlier is domain- and panel-dependent (the ratio in kinase,
   entropy in serotonin, the computational block in TP53, the CRP-protein pair in the indices). The
   existence of families is the invariant, not their membership.

4. **Field heuristics often approximate the minimal-commitment aggregate.** Where a field already uses an
   ad hoc summary, it frequently sits close to the weight-free average-rank aggregate. The TP53 community's
   median-of-eight-promoters rule proxies it at Spearman 0.987, and in kinase the Gini and entropy proxy it
   at about 0.95. The protocol here vindicates existing practice rather than overturning it, while also
   quantifying which measures stray furthest (the ratio in kinase at 0.796, the predictors in TP53 at 0.46
   to 0.60).

5. **Reproducibility is domain-dependent, and the test shows it.** The constitutive/artefactual test gives
   sharply different verdicts across domains. The disagreement structure of the inflammatory indices
   replicates across independent NHANES cycles at Spearman 0.995, while sepsis-signature structure
   replicates only weakly (0.32) and one signature's orientation flips between cohorts. This is the test
   doing its job, separating a stable concept-level disagreement from a partly sample-specific one, and it
   cautions against assuming any single instance's reproducibility generalises.

6. **External anchors confirm orientation and catch mis-orientation.** Across domains the consensus order
   agrees with a crude external anchor (Spearman about 0.65 to 0.94), more strongly at the extremes than
   in the middle where the anchor is continuous enough to split. The check has teeth. It flagged an
   un-gated kinase candidate that pointed backwards, and it flagged all four TP53 proliferation screens as
   mis-oriented under the naive fitness convention. Two anchors are training targets or coincide with the
   outcome used to orient (position clocks against age, sepsis signatures against mortality), so those legs
   are confirmatory rather than independent, and we say so.

These findings sit alongside a companion observation, established most cleanly in TP53 and in the QSAR
setting, that *evaluations* overstate performance when the benchmark overlaps the training signal. That is
a different failure from *definitions* disagreeing, and we treat it as a second face of the same
self-deception about measurement rather than folding it in.

The framing is not new in outline. That multiple operationalisations of one construct disagree is
convergent validity (Campbell and Fiske, 1959), and the raw "agree at the extremes, diverge in the middle"
pattern is the textbook Bland-Altman range effect (Bland and Altman, 1995). What is not generic is the
reproducibility test that promotes a disagreement from method variance to a property of the concept, and
the controls that decide which disagreement is which.

## 6. Discussion

**Reporting practice.** The protocol implies concrete changes to how competing measures are reported.
Declare each measure's type and scale, and verify its orientation against an external anchor, before any
comparison (step 1). Report the reproducible incomparable set and a desiderata table rather than a single
ranking (step 6). Validate a measure in its native domain, tissue, or panel rather than assuming transfer
(the reliability gate, Sections 4.4 and 4.6). When a single number is unavoidable, report the average-rank
canonical aggregate and each measure's distance from it, rather than an arbitrary favourite (step 7).

**The conceptual claim.** Some disagreement among biological measures is not a defect to be benchmarked
away. It is a faithful reflection of a concept that does not totally order its objects. This is the
departure from the psychometric reading of convergent validity, which treats mid-range divergence as
method variance to remove. The protocol's job is to find the part of the disagreement that is real and
reproducible, and to stop the field from adjudicating the part that is not there. Stated as a claim about
biology, several important concepts, including selectivity, biological age, variant impact, and
inflammatory burden, are only partially ordered by their own objects, and the persistent disagreement of
their measures is the observable trace of that partial order.

**Limitations.** Four are substantive. First, the base set is a choice. The consensus order, the
incomparable set, and the canonical aggregate are all defined relative to the set of measures admitted,
and there is no closed enumeration of "all measures" of a biological concept. We mitigate this with a
stated selection rule (measures in routine use, computable on the data, distinct in the order they
induce), with robustness checks, and by reporting the base set explicitly, but the choice is real and only
lifted one level, not removed. Second, the reproducibility evidence is uneven. It is strong for population
indices and thin or weak elsewhere, and one domain (sepsis) replicates only weakly, which we report rather
than hide. Third, two external anchors are not independent of the measures they check (chronological age
is the training target of position clocks, mortality is used to orient the sepsis signatures), so those
orientation legs are confirmatory. Fourth, the framework assumes a concept-level order exists and that
reasonable measures respect it. This is one substantive assumption, argued rather than proved. Separately,
categorical *endotype* systems are partitions rather than scalar orders, an adjacent problem for which
cluster-agreement methods, not the average-rank machinery, apply.

**Scope and relation to prior work.** The formalism is entirely borrowed. It is representational
measurement theory (Krantz et al., 1971), the partial-order theory of multi-indicator systems (Patil and
Taillie 2004, Brüggemann and Patil, 2011), Szpilrajn's extension theorem, and the average-rank canonical
extension (De Loof et al., 2008). The clearest prior application of measurement theory in biology (Houle et
al., 2011) treats single measurements and scale types rather than competing measures of one concept as
extensions of a shared partial order. The partial-order approach has precedent in environmental chemistry
(Halfon and Reggiani, 1986) and in social indicators (Fattore, 2016), and it contrasts with
rank-aggregation methods that force a consensus (Dwork et al., 2001, Kolde et al., 2012) rather than
preserve incomparability. An early biological instance of measures that share few components yet agree
where the concept is determined is the concordance among breast-cancer gene-expression predictors (Fan et
al., 2006). Our contribution is the protocol, the constitutive/artefactual
reproducibility test, and the cross-domain demonstration on public data. The value, for a working analyst,
is a reusable and honest procedure for the recurring situation of several defensible measures of one
latent concept that refuse to agree.

---

## Declarations

**Declaration of competing interest.** The author declares no competing financial or personal interests.

**Funding.** This research received no external funding.

**Declaration of generative AI in scientific writing.** During the preparation of this work the author
used Claude (Anthropic) to assist with data acquisition, analysis scripting, and drafting. The author
reviewed and edited all output and takes full responsibility for the accuracy and integrity of the work.

**Data and code availability.** All analyses use publicly available data and are fully reproducible from
the project repository, which contains, for each domain, a README documenting data provenance and the
measure-selection rule, the analysis scripts, and their stored outputs. *Repository URL and an archived
DOI (Zenodo) to be added on submission.*

**CRediT author contributions.** Polina Vinogradova: conceptualization, methodology, software, formal
analysis, investigation, writing – original draft, writing – review and editing.

---

## Appendix A — formal justification

This appendix states which established result justifies which step of the protocol, and fixes the object
so the steps are well-defined. Every result is classical or definitional. The contribution is the protocol
and the empirical work, not the theorems.

**Definition 1 (the object).** Fix a latent concept κ and a finite object set Ω. The empirical structure
is (Ω, ≽), where ≽ is a preorder (reflexive and transitive) that is not assumed total. Write a ~ b when
both a ≽ b and b ≽ a, write a ≻ b for strict order, and write a ∥ b when neither a ≽ b nor b ≽ a
(a and b are incomparable).

**Definition 2 (admissible measure).** A measure is a map f : Ω → ℝ. It is admissible if a ≽ b implies
f(a) ≥ f(b) and a ≻ b implies f(a) > f(b). This is the order-preserving representation condition of
representational measurement theory. Write 𝓕(≽) for the admissible measures.

**Definition 3 (families and linear extensions).** A measure f induces the total preorder a ⪰_f b iff
f(a) ≥ f(b). Two measures are in the same *family* when they induce the same total preorder. A total
preorder ⪰ is a *linear extension* of ≽ if a ≽ b implies a ⪰ b. On finite Ω, f ↦ ⪰_f maps 𝓕(≽) onto the
linear extensions of ≽. This is the multiple-indicators-to-partial-order setup of Patil and Taillie (2004).

**Consensus (how ≽ is estimated).** Given the set 𝓜 of measures in use, define a ≽𝓜 b iff f(a) ≥ f(b) for
every f in 𝓜. If every reasonable measure is admissible for the true ≽, then ≽ ⊆ ≽𝓜, so the consensus
over-declares comparability. The comparable pairs of ≽𝓜 are an upper bound on the truly comparable pairs,
and the incomparable pairs of ≽𝓜 are a lower bound on the truly incomparable pairs. Reproducibility across
independent cohorts (protocol step 3) is what licenses reading a reproducible incomparability as a property
of κ rather than a coincidence of one sample. An external anchor (step 1) is a second, independent check.

**Lemma 1 (disagreement lies in the incomparable set).** If f, g in 𝓕(≽) strictly disagree on a pair,
meaning f(a) > f(b) while g(a) < g(b), then a ∥ b. *Proof.* If a ≽ b then admissibility of g gives
g(a) ≥ g(b), a contradiction. If b ≽ a then admissibility of f gives f(b) ≥ f(a), a contradiction. So
a ∥ b. This is the contrapositive of admissibility. For the estimated order ≽𝓜 it holds by construction,
and its role is only to justify step 2.

**Lemma 2 (families are ordinal classes).** If f, g in 𝓕(≽) are in the same family, there is a strictly
increasing φ with g = φ ∘ f, and conversely any strictly increasing φ gives a same-family measure.
*Proof.* Define φ(f(a)) := g(a). It is well-defined and strictly increasing because f and g induce the
same order. So a family is an equivalence class under strictly monotone reparameterisation, an ordinal
class. Linear-relatedness on a dataset therefore does not establish a shared interval scale, and a
non-affine monotone relation between two measures keeps them in the same family. Family membership is
diagnosed by rank agreement, not by the regression form relating two measures.

**Proposition (non-uniqueness).** On finite (Ω, ≽): 𝓕(≽) is non-empty, the families are in bijection with
the linear extensions of ≽, and there is a unique family iff ≽ is total. *Proof.* Szpilrajn's extension
theorem gives a linear extension, which on a finite set has a real representation and is admissible. A
finite poset has a unique linear extension iff it is a chain. If a ∥ b, both a ≻ b and b ≻ a extensions
exist. So when incomparabilities exist there are at least two order-distinct admissible measures, and no
requirement that only entails order-preservation can single one out. This justifies step 6.

**The canonical aggregate.** When a single summary is required, the average rank of each object over the
linear extensions of ≽𝓜 (Patil and Taillie, 2004, De Loof et al., 2008) is the centroid of the admissible
class. It agrees with every admissible measure on the comparable pairs, and on each incomparable pair it
returns the fraction of extensions ranking one object above the other. It is weight-free and
parameter-free. Exact enumeration is #P-complete (Brightwell and Winkler, 1991), so it is estimated by
near-uniform sampling (Bubley and Dyer, 1999). It is canonical relative to the chosen base set 𝓜. This
justifies step 7.

## References

*Elsevier name-date (Harvard) style; DOIs verified against primary/publisher records. `[UNVERIFIED]` tags mark residual gaps to confirm at proof. Journal house style (including semicolon-separated in-text citation lists) is applied by the journal at proof.*

Anastassiadis, T., Deacon, S.W., Devarajan, K., Ma, H., Peterson, J.R., 2011. Comprehensive assay of kinase catalytic activity reveals features of kinase inhibitor selectivity. Nature Biotechnology 29, 1039–1045. https://doi.org/10.1038/nbt.2017

Arrow, K.J., 1963. Social Choice and Individual Values, 2nd ed. Wiley, New York.

Belsky, D.W., Caspi, A., Corcoran, D.L., Sugden, K., Poulton, R., Arseneault, L., Baccarelli, A., Chamarti, K., Gao, X., Hannon, E., Harrington, H.L., Houts, R., Kothari, M., Kwon, D., Mill, J., Schwartz, J., Vokonas, P., Wang, C., Williams, B.S., Moffitt, T.E., 2022. DunedinPACE, a DNA methylation biomarker of the pace of aging. eLife 11, e73420. https://doi.org/10.7554/eLife.73420

Bland, J.M., Altman, D.G., 1995. Comparing methods of measurement: why plotting difference against standard method is misleading. Lancet 346, 1085–1087. https://doi.org/10.1016/S0140-6736(95)91748-9

Bosc, N., Meyer, C., Bonnet, P., 2017. The use of novel selectivity metrics in kinase research. BMC Bioinformatics 18, 17. https://doi.org/10.1186/s12859-016-1413-y

Brandes, N., Goldman, G., Wang, C.H., Ye, C.J., Ntranos, V., 2023. Genome-wide prediction of disease variant effects with a deep protein language model. Nature Genetics 55, 1512–1522. https://doi.org/10.1038/s41588-023-01465-0

Brightwell, G., Winkler, P., 1991. Counting linear extensions. Order 8, 225–242. https://doi.org/10.1007/BF00383444

Brüggemann, R., Patil, G.P., 2011. Ranking and Prioritization for Multi-indicator Systems: Introduction to Partial Order Applications. Springer, New York.

Bubley, R., Dyer, M., 1999. Faster random generation of linear extensions. Discrete Mathematics 201, 81–88. https://doi.org/10.1016/S0012-365X(98)00333-1

Campbell, D.T., Fiske, D.W., 1959. Convergent and discriminant validation by the multitrait-multimethod matrix. Psychological Bulletin 56, 81–105. https://doi.org/10.1037/h0046016

Cheng, J., Novati, G., Pan, J., Bycroft, C., Žemgulytė, A., Applebaum, T., Pritzel, A., Wong, L.H., Zielinski, M., Sargeant, T., Schneider, R.G., Senior, A.W., Jumper, J., Hassabis, D., Kohli, P., Avsec, Ž., 2023. Accurate proteome-wide missense variant effect prediction with AlphaMissense. Science 381, eadg7492. https://doi.org/10.1126/science.adg7492

Davenport, E.E., Burnham, K.L., Radhakrishnan, J., Humburg, P., Hutton, P., Mills, T.C., Rautanen, A., Gordon, A.C., Garrard, C., Hill, A.V.S., Hinds, C.J., Knight, J.C., 2016. Genomic landscape of the individual host response and outcomes in sepsis: a prospective cohort study. Lancet Respiratory Medicine 4, 259–271. https://doi.org/10.1016/S2213-2600(16)00046-1

Davis, M.I., Hunt, J.P., Herrgard, S., Ciceri, P., Wodicka, L.M., Pallares, G., Hocker, M., Treiber, D.K., Zarrinkar, P.P., 2011. Comprehensive analysis of kinase inhibitor selectivity. Nature Biotechnology 29, 1046–1051. https://doi.org/10.1038/nbt.1990

Dwork, C., Kumar, R., Naor, M., Sivakumar, D., 2001. Rank aggregation methods for the web. In: Proceedings of the 10th International Conference on World Wide Web (WWW '01). ACM, New York, pp. 613–622. https://doi.org/10.1145/371920.372165

Fan, C., Oh, D.S., Wessels, L., Weigelt, B., Nuyten, D.S.A., Nobel, A.B., van't Veer, L.J., Perou, C.M., 2006. Concordance among gene-expression-based predictors for breast cancer. New England Journal of Medicine 355, 560–569. https://doi.org/10.1056/NEJMoa052933

Fattore, M., 2016. Partially ordered sets and the measurement of multidimensional ordinal deprivation. Social Indicators Research 128, 835–858. https://doi.org/10.1007/s11205-015-1059-6

Fortuno, C., Lee, K., Olivier, M., Pesaran, T., Mai, P.L., de Andrade, K.C., Attardi, L.D., Crowley, S., Evans, D.G., Feng, B.J., Foreman, A.K.M., Frone, M.N., Huether, R., James, P.A., McGoldrick, K., Mester, J., Seifert, B.A., Slavin, T.P., Witkowski, L., Zhang, L., Plon, S.E., Spurdle, A.B., Savage, S.A., 2021. Specifications of the ACMG/AMP variant interpretation guidelines for germline TP53 variants. Human Mutation 42, 223–236. https://doi.org/10.1002/humu.24152

Frazer, J., Notin, P., Dias, M., Gomez, A., Min, J.K., Brock, K., Gal, Y., Marks, D.S., 2021. Disease variant prediction with deep generative models of evolutionary data. Nature 599, 91–95. https://doi.org/10.1038/s41586-021-04043-8

Fridman, A.L., Tainsky, M.A., 2008. Critical pathways in cellular senescence and immortalization revealed by gene expression profiling. Oncogene 27, 5975–5987. https://doi.org/10.1038/onc.2008.213

Giacomelli, A.O., Yang, X., Lintner, R.E., McFarland, J.M., Duby, M., Kim, J., Howard, T.P., Takeda, D.Y., Ly, S.H., Kim, E., Gannon, H.S., Hurhula, B., Sharpe, T., Goodale, A., Fritchman, B., Steelman, S., Vazquez, F., Tsherniak, A., Aguirre, A.J., Doench, J.G., Piccioni, F., Roberts, C.W.M., Meyerson, M., Getz, G., Johannessen, C.M., Root, D.E., Hahn, W.C., 2018. Mutational processes shape the landscape of TP53 mutations in human cancer. Nature Genetics 50, 1381–1387. https://doi.org/10.1038/s41588-018-0204-y

Halfon, E., Reggiani, M.G., 1986. On ranking chemicals for environmental hazard. Environmental Science & Technology 20, 1173–1179. https://doi.org/10.1021/es00154a010

Houle, D., Pélabon, C., Wagner, G.P., Hansen, T.F., 2011. Measurement and meaning in biology. Quarterly Review of Biology 86, 3–34. https://doi.org/10.1086/658408

Kalyakulina, A., Yusipov, I., Kondakova, E., Bacalini, M.G., Franceschi, C., Vedunova, M., Ivanchenko, M., 2023. Small immunological clocks identified by deep learning and gradient boosting. Frontiers in Immunology 14, 1177611. https://doi.org/10.3389/fimmu.2023.1177611

Kato, S., Han, S.Y., Liu, W., Otsuka, K., Shibata, H., Kanamaru, R., Ishioka, C., 2003. Understanding the function-structure and function-mutation relationships of p53 tumor suppressor protein by high-resolution missense mutation analysis. Proceedings of the National Academy of Sciences of the United States of America 100, 8424–8429. https://doi.org/10.1073/pnas.1431692100

Klaeger, S., Heinzlmeir, S., Wilhelm, M., Polzer, H., Vick, B., Koenig, P.A., Reinecke, M., Ruprecht, B., Petzoldt, S., Meng, C., Zecha, J., Reiter, K., Qiao, H., Helm, D., Koch, H., Schoof, M., Canevari, G., Casale, E., Re Depaolini, S., Feuchtinger, A., Wu, Z., Schmidt, T., Rueckert, L., Becker, W., Huenges, J., Garz, A.K., Gohlke, B.O., Zolg, D.P., Kayser, G., Vooder, T., Preissner, R., Hahne, H., Tõnisson, N., Kramer, K., Götze, K., Bassermann, F., Schlegl, J., Ehrlich, H.C., Aiche, S., Walch, A., Greif, P.A., Schneider, S., Felder, E.R., Ruland, J., Médard, G., Jeremias, I., Spiekermann, K., Kuster, B., 2017. The target landscape of clinical kinase drugs. Science 358, eaan4368. https://doi.org/10.1126/science.aan4368

Kolde, R., Laur, S., Adler, P., Vilo, J., 2012. Robust rank aggregation for gene list integration and meta-analysis. Bioinformatics 28, 573–580. https://doi.org/10.1093/bioinformatics/btr709

Kotler, E., Shani, O., Goldfeld, G., Lotan-Pompan, M., Tarcic, O., Gershoni, A., Hopf, T.A., Marks, D.S., Oren, M., Segal, E., 2018. A systematic p53 mutation library links differential functional impact to cancer mutation pattern and evolutionary conservation. Molecular Cell 71, 178–190.e8. https://doi.org/10.1016/j.molcel.2018.06.012

Krantz, D.H., Luce, R.D., Suppes, P., Tversky, A., 1971. Foundations of Measurement, Vol. I: Additive and Polynomial Representations. Academic Press, New York.

Landrum, M.J., Lee, J.M., Benson, M., Brown, G.R., Chao, C., Chitipiralla, S., Gu, B., Hart, J., Hoffman, D., Jang, W., Karapetyan, K., Katz, K., Liu, C., Maddipatla, Z., Malheiro, A., McDaniel, K., Ovetsky, M., Riley, G., Zhou, G., Holmes, J.B., Kattman, B.L., Maglott, D.R., 2018. ClinVar: improving access to variant interpretations and supporting evidence. Nucleic Acids Research 46, D1062–D1067. https://doi.org/10.1093/nar/gkx1153

Liberzon, A., Birger, C., Thorvaldsdóttir, H., Ghandi, M., Mesirov, J.P., Tamayo, P., 2015. The Molecular Signatures Database hallmark gene set collection. Cell Systems 1, 417–425. https://doi.org/10.1016/j.cels.2015.12.004

Luce, R.D., Krantz, D.H., Suppes, P., Tversky, A., 1990. Foundations of Measurement, Vol. III: Representation, Axiomatization, and Invariance. Academic Press, San Diego.

Mavrommatis, C., Belsky, D.W., Ying, K., Moqri, M., Campbell, A., Richmond, A., Gladyshev, V.N., Chandra, T., McCartney, D.L., Marioni, R.E., 2025. An unbiased comparison of 14 epigenetic clocks in relation to 174 incident disease outcomes. Nature Communications 16, 11164. https://doi.org/10.1038/s41467-025-66106-y [UNVERIFIED: middle-author list (Ying, Moqri, Campbell, Richmond, Gladyshev, Chandra) drawn from open-access listing, not confirmed against the trusted source, which abbreviated it; year/venue/volume/article-number/DOI verified]

Metz, J.T., Johnson, E.F., Soni, N.B., Merta, P.J., Kifle, L., Hajduk, P.J., 2011. Navigating the kinome. Nature Chemical Biology 7, 200–202. https://doi.org/10.1038/nchembio.530

Notin, P., Kollasch, A., Ritter, D., van Niekerk, L., Paul, S., Spinner, H., Rollins, N., Shaw, A., Orenbuch, R., Weitzman, R., Frazer, J., Dias, M., Franceschi, D., Gal, Y., Marks, D.S., 2023. ProteinGym: large-scale benchmarks for protein fitness prediction and design. In: Advances in Neural Information Processing Systems 36 (NeurIPS 2023), Datasets and Benchmarks Track, pp. 64331–64379.

Patil, G.P., Taillie, C., 2004. Multiple indicators, partially ordered sets, and linear extensions: multi-criterion ranking and prioritization. Environmental and Ecological Statistics 11, 199–228. https://doi.org/10.1023/B:EEST.0000027209.93218.d9

Saul, D., Kosinsky, R.L., Atkinson, E.J., Doolittle, M.L., Zhang, X., LeBrasseur, N.K., Pignolo, R.J., Robbins, P.D., Niedernhofer, L.J., Ikeno, Y., Jurk, D., Passos, J.F., Hickson, L.T.J., Xue, A., Monroe, D.G., Tchkonia, T., Kirkland, J.L., Farr, J.N., Khosla, S., 2022. A new gene set identifies senescent cells and predicts senescence-associated pathways across tissues. Nature Communications 13, 4827. https://doi.org/10.1038/s41467-022-32552-1

Sayed, N., Huang, Y., Nguyen, K., Krejciova-Rajaniemi, Z., Grawe, A.P., Gao, T., Tibshirani, R., Hastie, T., Alpert, A., Cui, L., Kuznetsova, T., Rosenberg-Hasson, Y., Ostan, R., Monti, D., Lehallier, B., Shen-Orr, S.S., Maecker, H.T., Dekker, C.L., Wyss-Coray, T., Franceschi, C., Jojic, V., Haddad, F., Montoya, J.G., Wu, J.C., Davis, M.M., Furman, D., 2021. An inflammatory aging clock (iAge) based on deep learning tracks multimorbidity, immunosenescence, frailty and cardiovascular aging. Nature Aging 1, 598–615. https://doi.org/10.1038/s43587-021-00082-y

Scicluna, B.P., van Vught, L.A., Zwinderman, A.H., Wiewel, M.A., Davenport, E.E., Burnham, K.L., Nürnberg, P., Schultz, M.J., Horn, J., Cremer, O.L., Bonten, M.J., Hinds, C.J., Wong, H.R., Knight, J.C., van der Poll, T., MARS Consortium, 2017. Classification of patients with sepsis according to blood genomic endotype: a prospective cohort study. Lancet Respiratory Medicine 5, 816–826. https://doi.org/10.1016/S2213-2600(17)30294-1

Shavlakadze, T., Xiong, K., Donne, R., Glass, D.J., et al., 2026. Multitissue, multi-time point transcriptomic atlas of aging in mice and rats. Science Advances. https://doi.org/10.1126/sciadv.ady8401  [UNVERIFIED: complete author list and volume/eLocator to confirm at proof]

Singleton, J., Booth, R., 2022. Towards an axiomatic approach to truth discovery. Autonomous Agents and Multi-Agent Systems 36, 42. https://doi.org/10.1007/s10458-022-09569-3

Suppes, P., Krantz, D.H., Luce, R.D., Tversky, A., 1989. Foundations of Measurement, Vol. II: Geometrical, Threshold, and Probabilistic Representations. Academic Press, San Diego.

Szpilrajn, E., 1930. Sur l'extension de l'ordre partiel. Fundamenta Mathematicae 16, 386–389. https://doi.org/10.4064/fm-16-1-386-389

Uitdehaag, J.C.M., Zaman, G.J.R., 2011. A theoretical entropy score as a single value to express inhibitor selectivity. BMC Bioinformatics 12, 94. https://doi.org/10.1186/1471-2105-12-94

Yusipov, I., Kondakova, E., Kalyakulina, A., Krivonosov, M., Lobanova, N., Bacalini, M.G., Franceschi, C., Vedunova, M., Ivanchenko, M., 2022. Accelerated epigenetic aging and inflammatory/immunological profile (ipAGE) in patients with chronic kidney disease. GeroScience 44, 817–834. https://doi.org/10.1007/s11357-022-00540-4
