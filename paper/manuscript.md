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
partial order the indicators jointly induce (Patil and Taillie 2004, Brüggemann and Patil 2011). We do not
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
assay technologies (Davis 68 compounds, Klaeger 222, Anastassiadis 178, Metz 704). The measures split into
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
samples, GSE40279 and GSE87571). Position clocks are affine-related (R² 0.991 to 0.995), and rate clocks
relate to position through a signed logarithm (R² 0.93 to 0.95), a Weber-Fechner-like relation that is a
descriptive observation, not a scale-type proof (Section 2). All clocks fail an ordinal-stability check
against cell-type composition (absolute correlation up to 0.372), and GrimAge, a two-stage composite, is
structurally near-orthogonal to the others (cosine of coefficient vectors between −0.17 and 0.11). The
external anchor here, chronological age, is the training target of the position clocks, so for them the
orientation check is confirmatory rather than independent. This is stated as a limitation.

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
transactivation readouts, seven computational pathogenicity predictors (AlphaMissense, EVE, ESM1b, REVEL,
CADD, PrimateAI, BayesDel), and mammalian proliferation screens (Giacomelli, Kotler), with ClinVar as the
external clinical anchor (core panel 2,314 variants). Step 1 earns its place here. The eight
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
(strong disagreement, weak reproducibility):** five published blood signatures (Hallmark inflammatory,
Hallmark interferon-γ, SRS7, SRSq19, MARS8) on GSE65682 (802 samples), anchor 28-day mortality. The
Hallmark and SRS families are near-independent (cross-family mean +0.02), so two published sepsis-severity
signatures barely agree, but the detailed structure replicates only weakly across the discovery and
validation cohorts (Spearman 0.32), which we report as a limitation rather than a success. **Inflammaging
clocks (thin):** only two composite clocks are publicly reconstructable on shared samples (SImAge, and
ipAGE recomputed), so this is a qualitative check, not a full instance. Sepsis *endotype labels* (SRS
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
Taillie 2004, Brüggemann and Patil 2011), Szpilrajn's extension theorem, and the average-rank canonical
extension (De Loof et al. 2008). The clearest prior application of measurement theory in biology (Houle et
al. 2011) treats single measurements and scale types rather than competing measures of one concept as
extensions of a shared partial order. Our contribution is the protocol, the constitutive/artefactual
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
linear extensions of ≽𝓜 (Patil and Taillie 2004, De Loof et al. 2008) is the centroid of the admissible
class. It agrees with every admissible measure on the comparable pairs, and on each incomparable pair it
returns the fraction of extensions ranking one object above the other. It is weight-free and
parameter-free. Exact enumeration is #P-complete (Brightwell and Winkler, 1991), so it is estimated by
near-uniform sampling (Bubley and Dyer, 1999). It is canonical relative to the chosen base set 𝓜. This
justifies step 7.

## References

*Author-year list below. The methodological references are verified in `related-work.md`. The
domain/data-source references are the standard citations for each dataset and need a final page/year
verification and conversion to JTB's reference style before submission. Two items flagged in
`related-work.md` to confirm against the primary source: Foundations of Measurement Vol. II year (1989
vs 1990), and Michell (1997) page range.*

**Foundations (measurement theory, order theory, aggregation).**
- Arrow, K.J. (1963). *Social Choice and Individual Values*, 2nd ed. Wiley.
- Bland, J.M. and Altman, D.G. (1995). Comparing methods of measurement. *Lancet* 346, 1085–1087.
- Brightwell, G. and Winkler, P. (1991). Counting linear extensions. *Order* 8, 225–242.
- Brüggemann, R. and Patil, G.P. (2011). *Ranking and Prioritization for Multi-indicator Systems*. Springer.
- Bubley, R. and Dyer, M. (1999). Faster random generation of linear extensions. *Discrete Mathematics*
  201, 81–88.
- Campbell, D.T. and Fiske, D.W. (1959). Convergent and discriminant validation by the
  multitrait-multimethod matrix. *Psychological Bulletin* 56(2), 81–105.
- De Loof, K., De Baets, B., De Meyer, H. and Brüggemann, R. (2008). A Hitchhiker's Guide to Poset Ranking.
  *Comb. Chem. High Throughput Screen.* 11(9), 734–744.
- De Loof, K., De Baets, B. and De Meyer, H. (2011). Approximation of Average Ranks in Posets.
  *MATCH Commun. Math. Comput. Chem.* 66(1), 219–229.
- Dwork, C., Kumar, R., Naor, M. and Sivakumar, D. (2001). Rank aggregation methods for the web. *WWW '01*.
- Fattore, M. (2016). Partially Ordered Sets and the Measurement of Multidimensional Ordinal Deprivation.
  *Social Indicators Research* 128(2), 835–858.
- Halfon, E. and Reggiani, M.G. (1986). On ranking chemicals for environmental hazard. *Environ. Sci.
  Technol.* 20(11), 1173–1179.
- Houle, D., Pélabon, C., Wagner, G.P. and Hansen, T.F. (2011). Measurement and Meaning in Biology.
  *Quarterly Review of Biology* 86(1), 3–34.
- Krantz, D.H., Luce, R.D., Suppes, P. and Tversky, A. (1971). *Foundations of Measurement, Vol. I*.
  Academic Press.
- Kolde, R., Laur, S., Adler, P. and Vilo, J. (2012). Robust rank aggregation for gene list integration.
  *Bioinformatics* 28(4), 573–580.
- Lerche, D., Sørensen, P.B. and Brüggemann, R. (2003). Improved estimation of the ranking probabilities
  in partial orders. *J. Chem. Inf. Comput. Sci.* 43(5), 1471–1480.
- Luce, R.D., Krantz, D.H., Suppes, P. and Tversky, A. (1990). *Foundations of Measurement, Vol. III*.
  Academic Press.
- Michell, J. (1999). *Measurement in Psychology: A Critical History of a Methodological Concept*.
  Cambridge University Press.
- Patil, G.P. and Taillie, C. (2004). Multiple indicators, partially ordered sets, and linear extensions.
  *Environmental and Ecological Statistics* 11(2), 199–228.
- Singleton, J. and Booth, R. (2022). Towards an axiomatic approach to truth discovery. *Autonomous Agents
  and Multi-Agent Systems* 36(2), 42.
- Suppes, P., Krantz, D.H., Luce, R.D. and Tversky, A. (1989). *Foundations of Measurement, Vol. II*.
  Academic Press.
- Szpilrajn, E. (1930). Sur l'extension de l'ordre partiel. *Fundamenta Mathematicae* 16, 386–389.

**Domain precedents and data sources.**
- Belsky, D.W. et al. (2022). DunedinPACE, a DNA methylation biomarker of the pace of aging. *eLife* 11,
  e73420.
- Bosc, N. et al. (2017). Large-scale comparison of kinase selectivity metrics. *BMC Bioinformatics* 18, 17.
- Cheng, J. et al. (2023). Accurate proteome-wide missense variant effect prediction with AlphaMissense.
  *Science* 381, eadg7492.
- Davenport, E.E. et al. (2016). Genomic landscape of the individual host response and outcomes in sepsis
  (SRS). *Lancet Respiratory Medicine* 4, 259–271.
- Fan, C. et al. (2006). Concordance among gene-expression-based predictors for breast cancer.
  *New England Journal of Medicine* 355(6), 560–569.
- Fortuno, C. et al. (2021). Specifications of the ACMG/AMP variant interpretation guidelines for germline
  TP53. *Human Mutation* 42, 223–236.
- Furman, D. et al. (2021). An inflammatory aging clock (iAge). *Nature Aging* 1, 598–615.
- Giacomelli, A.O. et al. (2018). Mutational processes shape the landscape of TP53 mutations in human
  cancer. *Nature Genetics* 50, 1381–1387.
- Kato, S. et al. (2003). Understanding the function-structure and function-mutation relationships of p53.
  *PNAS* 100, 8424–8429.
- Kotler, E. et al. (2018). A systematic p53 mutation library links differential functional impact to
  cancer mutation pattern. *Molecular Cell* 71, 178–190.
- Mavrommatis, T. et al. (2025). An unbiased comparison of 14 epigenetic clocks. *Nature Communications*
  16, 11164.
- Notin, P. et al. (2023). ProteinGym: large-scale benchmarks for protein fitness prediction. *NeurIPS
  Datasets and Benchmarks*.
- Scicluna, B.P. et al. (2017). Classification of patients with sepsis according to blood genomic endotype
  (MARS). *Lancet Respiratory Medicine* 5, 816–826.
- Uitdehaag, J.C.M. and Zaman, G.J.R. (2011). A theoretical entropy score as a single value to express
  inhibitor selectivity. *BMC Bioinformatics* 12, 94.
- *Additional data sources to cite: NHANES (CDC, 2015–2016 and 2017–2018 cycles) and its Linked Mortality
  File, ClinVar (Landrum et al. 2018), dbNSFP (Liu et al.), EVE (Frazer et al. 2021), ESM1b (Brandes et
  al. 2023), SImAge and ipAGE (Kalyakulina et al. 2023 and 2022), MSigDB Hallmark gene sets (Liberzon et
  al. 2015). The four kinase datasets: Davis et al. 2011, Klaeger et al. 2017, Anastassiadis et al. 2011,
  Metz et al. 2011.*
