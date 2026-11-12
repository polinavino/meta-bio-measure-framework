# A measurement-theoretic protocol for competing quantifications of latent biological concepts

*Polina Vinogradova*
Independent Researcher.
Correspondence: polina.vino@gmail.com. ORCID: 0000-0003-3271-3841.

> **Target:** *Journal of Theoretical Biology* (hybrid, traditional/subscription route, no APC).
> Elsevier "Your Paper Your Way" allows flexible formatting at initial submission. A formatting pass to
> JTB style (reference format, figure files) happens after the scientific draft is complete. This file
> is the canonical submission draft. `main.md` and `formal-spine.md` are the working notes it draws on.
> Quantitative claims draw on stored outputs in `analysis/outputs/`, `tp53/analysis/outputs/`,
> `inflammation/analysis/outputs/` and `cvdrisk/analysis/outputs/`. `bash analysis/run_framework.sh`
> reproduces the local workflow, figures, registered numerical checks and post-remediation audit.
> Figure values are parsed from stored outputs. `paper/check_claims.py` checks registered quantities;
> it does not independently establish every biological interpretation or reference.

---

## Highlights
- Under admissibility, competing measures are order-preserving extensions of a shared partial order.
- So disagreement is expected, and the output should say which pairs the concept orders, not pick a winner.
- A protocol tests reproducible disagreement while controlling implementation, coverage and sampling artifacts.
- Demonstrated on seven domains: kinase and serotonin selectivity, ageing, exposure, TP53, inflammation,
  cardiovascular risk.
- Replication establishes stability; concept-level attribution requires additional assumptions and controls.

## Abstract
Many biological concepts are quantified by several formulas that can produce distinct rankings. We give a
measurement-theoretic protocol for comparing those formulas without assuming that one must be the correct
measure. Conditional on admissibility, order-preserving measures can reverse only pairs the concept leaves
incomparable; admissibility itself is a substantive assumption, not a consequence of agreement or replication.
The protocol declares the biological object, formula implementation, measurement panel, missing-data convention
and eligibility gate, computes weak consensus orders with ties preserved, and tests their stability across
cohorts and analysis choices. We distinguish reproducible disagreement from its attribution: shared biases,
threshold effects and coverage artifacts can also replicate. Public-data demonstrations across seven domains
show how implementation and population choices change apparent agreement clusters, how separation patterns
depend on the chosen axis, and how strong conditional rank recovery can coexist with poor eligibility
retention. When a scalar summary is useful, uniform average ranks of consensus-order equivalence classes are
one explicit aggregation convention, estimated with dispersed-chain convergence and Monte Carlo precision
checks. The output is the observed comparative structure, its sensitivity and uncertainty, and the evidence
still needed for concept-level interpretation.

**Keywords:** measurement theory, partially ordered sets, convergent validity, rank aggregation,
biomarkers, variant effect prediction

---

## 1. Introduction

For many biological concepts, choosing how to measure them is part of defining them. To say that compound A
is more selective than B, that an intervention slowed someone's ageing, or that a genetic variant is more
damaging than another, one must first choose a number to stand for the concept. For most interesting
concepts several defensible numbers are in routine use. These are different formulas, not implementations
of one agreed formula, and on the same data they rank the same objects differently. Faced with this, the
literature typically benchmarks the candidates against some reference and declares a winner.

We distinguish outcome-specific prediction benchmarking from comparison of biological definitions. An
established view of measurement provides a framework for the latter. In representational
measurement theory (Krantz et al., 1971, Suppes et al., 1989, Luce et al., 1990), a measure is an
order-preserving map from a concept's comparative structure into the real numbers. When several indicators
target one concept, the orders consistent with all of them are exactly the *extensions* of the
partial order the indicators jointly induce (Patil and Taillie, 2004, Brüggemann and Patil, 2011). We
bring this view to biological measurement as a usable protocol.

Two consequences of the view are immediate and classical. First, any two admissible measures must agree
wherever the concept's order is determined, and can disagree only on pairs the concept leaves
**incomparable** (Section 2, Lemma 1). Second, whenever incomparabilities exist, order preservation permits a whole *class* of admissible ordinal families, one per weak extension
(Section 2), and does not single out a unique one. This is a non-uniqueness
result, not an impossibility of satisfying all criteria; additional substantive criteria may distinguish
particular measures. The output should therefore report
which pairs the concept orders, which it leaves incomparable, and, if a single number is unavoidable, which
explicit aggregation convention is being used.

This reframing matters because the alternative is misleading. Benchmarking competing measures to a winner
treats concept-level ambiguity as estimation error to be resolved, when it may instead be a correct
description of a concept that does not totally order its objects. A TP53 variant with partial,
promoter-selective loss of function, a former smoker between clearly-exposed and clearly-unexposed, a
compound of intermediate selectivity: the concept itself may not order these objects. Reporting a single
ranking hides these cases.

**What is new here** is not the formalism, which is textbook and which we cite as foundation, but:

1. **A protocol** (Section 3) that operationalizes this view for a working analyst who has several measures
   and some cohorts, including the controls needed to avoid over-interpretation.
2. **A separation of reproducibility from causal interpretation.** The situation is a *convergent-validity*
   problem (Campbell and Fiske, 1959). Such divergence can be treated as method variance to be removed, but that is not the only possible
   interpretation. We distinguish reproducible disagreement from a *constitutive* interpretation, which additionally
   assumes admissibility and requires controls for systematic artifacts (Section 3, steps 3 and 4).
3. **A cross-domain demonstration** on public data (Section 4) spanning kinase pharmacology, epigenetic
   ageing, serotonin pharmacology, methylation biomarkers, TP53 variant effect, inflammatory burden, and
   cardiovascular risk.
   Several of the seven domains yield null or negative results.

Persistent disagreement establishes a stable difference between the chosen operationalisations. It supports
a concept-level partial-order interpretation only under the stated admissibility assumptions and after
plausible shared artifacts have been examined.

Section 2 states the measurement-theoretic foundations and what follows from each. Section 3 gives the
protocol. Section 4 works through the seven domains. Section 5 reports what the demonstrations establish,
with scope. Section 6 discusses reporting practice and limits. Appendix A collects the formal statements.

---

## 2. Foundations

We fix the object and record which established result justifies which step of the protocol. Everything in
this section is classical or definitional. Fuller statements and proofs are in Appendix A.

**The object.** Fix a latent concept κ (selectivity, biological age, variant impact, inflammatory burden) and
a finite set of objects Ω. Domain knowledge licenses *comparative judgments* a ≽ b ("a is at least as κ as
b"). The empirical structure is (Ω, ≽) with ≽ a preorder that is **not** assumed total. Write a ∥ b when
neither direction is supported, meaning a and b are *incomparable*. A *measure* f : Ω → ℝ is **admissible** if
a ≽ b implies f(a) ≥ f(b) and a ≻ b implies f(a) > f(b), an order-preserving homomorphism, which is the
representation condition of measurement theory. Two measures belong to the same exact **ordinal family** if they induce the
same total preorder. Empirical high-correlation groups are called **agreement clusters** below; they need
not be members of one exact ordinal family. A family is an equivalence class under strictly monotone reparameterization, that is, an
*ordinal* class (Lemma 2). We adopt the multiple-indicators-to-partial-order-to-linear-extensions setup of
Patil and Taillie (2004).

**Disagreement lies in the incomparable set (Lemma 1).** If two admissible measures rank a pair in
opposite strict directions, that pair is incomparable under κ. Equivalently, admissible measures agree on
the strict order of every comparable pair, and all inter-measure reversals lie among the incomparabilities.
This justifies protocol step 2.

**The order is estimated by consensus, and the estimate over-declares comparability.** The true ≽ is not
observed. Given the set 𝓜 of measures experts actually use, define the consensus order a ≽𝓜 b iff
f(a) ≥ f(b) for every f in 𝓜. Under the assumption that every reasonable measure is admissible, ≽ ⊆ ≽𝓜,
so the consensus can mark a pair as comparable when the measures agree though κ leaves
it incomparable. Replication tests stability conditional on the formulas and samples, not the admissibility
assumption or the absence of shared bias. An external anchor can check direction, but its independence from
inputs, training targets and fitted sign choices must be stated. Active-target count, for example, is the
S-score numerator and is a shared-input proxy, not an independent external anchor.

**Order preservation alone does not select a unique measure.** On the finite quotient of (Ω, ≽) by
its equivalence relation, Szpilrajn's extension theorem gives a linear extension and hence an admissible
measure. The admissible ordinal families correspond to the weak extensions, with a unique family iff ≽ is
total. If a pair is incomparable, both strict directions occur in admissible extensions. This establishes
non-uniqueness under order preservation alone, without excluding additional domain criteria that select a
particular construction. Benchmarking can answer a stated prediction or decision question; it cannot, from
rank agreement alone, identify the correct definition of a concept.

**An explicit scalar-summary convention.** When a single score is useful, we report the average rank of each
equivalence class over the linear extensions of the consensus quotient poset, under the uniform distribution
on those extensions. Objects tied on every admitted measure form one class and receive the same score. The
aggregate preserves comparable class orderings but may order incomparable classes by their mean positions.
The uniform law is an additional aggregation convention, not a consequence of admissibility, and does not
imply that each incomparable pair has probability one-half in each direction. The mutual-rank-probability
matrix retains more information about those pairs. No measure-specific weights are required, but the base
set, class treatment and distribution over extensions remain explicit choices. Counting linear extensions
is #P-complete (Brightwell and Winkler, 1991). We estimate average class ranks with an insertion Gibbs
sampler, checking four dispersed chains, split R-hat and batch-means Monte Carlo error; same-start chain
correlations are insufficient. Exact enumeration checks small cases. Failed diagnostics suppress aggregate
benchmarks (Section 3, step 7).

**A scale-type caution.** Families are ordinal classes only (Lemma 2). Two admissible measures being
linearly related on one dataset does not establish a shared interval scale. A linear relation is expected,
for instance, when both are regressions on a shared external target. A non-affine monotone relation between
two measures, such as rate versus position clocks, is a descriptive observation rather than a scale-type
proof. The protocol reports relations between families descriptively and does not assign interval or ratio
structure it cannot establish.

---

## 3. The protocol

Input: a latent concept κ, a set Ω of objects, candidate measures f₁ … f_k each claiming to quantify κ,
and ideally at least two independent cohorts. Each step names the established result that justifies it
(Section 2, Appendix A).

1. **Declare the object, type, scale and implementation; check orientation.** State the biological unit,
   target/feature panel, outcome definition, eligibility gate, exact formula, weights and parameter choices.
   Separate unmeasured, censored and measured-low entries, and state how undefined scores are handled.
   Compare published and modified implementations explicitly. Fix signs by declared conventions or on a
   discovery sample, then freeze them for validation. Record whether an anchor is external metadata, a shared
   input, a training target, or the outcome used to orient the scores. Correlation with a dependent anchor is
   a consistency check, not independent validation.
2. **Compute the weak consensus order, preserving ties.** A pair is incomparable iff at least one admitted
   measure ranks it strictly in each direction. A tie versus a strict comparison remains weakly comparable
   and is reported separately. Collapse identical indicator vectors into equivalence classes before
   computing linear extensions. Verify invariance to object-row permutations and strictly increasing
   reparameterisations.
3. **Test reproducibility without identifying its cause by replication alone.** Freeze implementations,
   eligibility and any fitted signs. On matched objects, compare pair relations directly; on disjoint
   populations, compare correlation structures and separation curves without claiming replication of the
   same incomparable pairs. Bootstrap independent objects or participants, retaining dependence among their
   scores. Replication establishes stability conditional on these choices. Shared inputs, assay limits,
   selection and systematic biases can also replicate; a constitutive interpretation needs additional
   evidence that the measures are admissible and the disagreement survives plausible artifact controls.
4. **Locate disagreement with explicit controls.** Examine separation, axis position, local density,
   distance to decision thresholds and rank boundaries. Declare which scores share inputs with the axis.
   Withhold a score that defines the proxy axis from the discordance panel as a sensitivity analysis; where
   possible also use an external proxy. Consensus-axis declines are descriptive associations, not established
   upper bounds on a constitutive component. Even an external-axis residual does not by itself identify a
   biological cause. Pairwise regressions are descriptive unless uncertainty accounts for shared objects.
5. **Describe agreement clusters and evaluate criteria.** Distinguish approximate correlation clusters from
   exact ordinal families. Separate properties proved by construction from empirical robustness and recovery.
   Declare apparatus context and analyst choices; test plausible alternative contexts without demanding
   invariance where the scientific question changes. Report the domain, scale and boundary conditions under
   which each guarantee holds. Boundary locality and stability across boundary crossings are different
   properties.
6. **Report structure, sensitivity and uncertainty.** Give the base measures, comparable and incomparable
   fractions, tie-only disagreements, eligibility counts, and leave-one-measure-out sensitivity. Report
   published-form and missingness scenarios. For reduced panels report eligibility retention, rank recovery
   conditional on eligibility, and pair recovery that gives lost eligibility zero credit. Use a benchmark
   winner only for a declared prediction or decision objective, with an appropriate endpoint.
7. **If a scalar summary is useful, state its convention and validate its computation.** Average ranks of
   equivalence classes under uniform linear extensions are one option. Use exact enumeration for small
   posets or a sampler with dispersed starts, checks of rescaling/row invariance, split R-hat, Monte Carlo
   errors and increasing run lengths. Here the numerical acceptance rules are maximum split R-hat <1.01
   and maximum batch-means Monte Carlo SE <0.0025 times the number of classes. These are empirical checks,
   not a proof of finite-time mixing. Suppress comparisons when diagnostics fail and avoid distinguishing
   close proxies beyond numerical precision. Distance from this aggregate measures agreement with the
   chosen convention, not overcommitment or biological error.

Steps 1 to 4 distinguish implementation, eligibility, sampling and axis effects. Steps 5 to 7 make the
remaining choices and numerical uncertainty explicit (Section 5).

---

## 4. Worked demonstrations

We apply the protocol in seven domains, all on public data. The domains differ in maturity and in how well
the protocol's assumptions hold. Kinase selectivity is the running example and the most
developed. TP53 variant effect and inflammatory burden are the newest and most data-rich. Every number
below is copied from a stored analysis output in the project repository. Table 1 summarises the instances and
the coverage of the protocol's steps across them.

**Table 1.** The ten illustrative instances considered by the protocol, grouped into the seven domains of this section.
*Incomparable* is the fraction of object pairs the consensus poset leaves unordered. *Best proxy* is the
existing measure whose ranking correlates most strongly (Spearman) with the uniform average class rank over
linear extensions. A dash marks a quantity not computed for that instance: the poset and average-rank steps
were run where the object set was large enough to support them, and the remaining instances were used for the
agreement-cluster and orientation steps only. The inflammatory-index poset is computed on a fixed random
subsample because pair counts scale quadratically. Every value is copied from a stored analysis output. Reproducibility summaries have different estimands:
matched-pair overlap is reported for kinase, whereas matrix-vector correlations compare dependent
indicator-pair correlations across disjoint cohorts; these should not be read as interchangeable tests.

| Instance | Measures | Objects (poset) | Incomparable | Best proxy | Cross-cohort reproducibility |
|---|---|---|---|---|---|
| Kinase selectivity (4.1) | 4 | 200 gated compounds | 37.9% | Gini, 0.949 | matched-panel incomparable-set Jaccard 0.294 |
| Blood clocks (4.2) | 5 | 1,385 samples | — | — | not tested |
| Retinal clocks (4.2) | 5 | 85 samples | — | — | single cohort |
| Serotonin selectivity (4.3) | 4 | 285 eligible of 297 ligands | 75.4% | — | not tested |
| Smoking signatures (4.4) | 4 | 464 and 689 samples | — | — | divergence replicates in both |
| TP53 variant effect (4.5) | 9 | 2,314 variants | 56.6% | median-of-8 promoters, 0.955 | platforms agree at the extremes |
| Inflammatory indices (4.6) | 7 | 600 of 5,327 mortality-linked participants | 72.8% | SII, 0.844 | matrix-vector Spearman 0.995 |
| Sepsis gene-set proxies (4.6) | 5 | 479 samples | 87.1% | inflammatory, 0.613 | frozen-sign matrix-vector Spearman 0.976 |
| Inflammaging clocks (4.6) | 5 | 343 samples | 54.1% | — | qualitative only |
| Cardiovascular risk (4.7) | 3 | 600 of 2,994 adults | 15.5% | — | mean Δr 0.004 |


### 4.1 Kinase-inhibitor selectivity
Concept: concentration of an inhibitor's activity across a declared protein-kinase panel. The four sources
span three assay technologies: Davis et al. (2011, 68 compounds), Klaeger et al. (2017, 222), Anastassiadis
et al. (2011, 178) and Metz et al. (2011, 704). Following the revised companion analysis, Davis uses 363
reference targets rather than 433 constructs, Klaeger 217 source-annotated protein kinases rather than 343
captured targets, and Anastassiadis 291 non-duplicate entries rather than 300. Targets, constructs and
co-captured proteins are distinct biological units; curation is part of the measurement definition.

The primary Klaeger consensus uses the S-score, baseline-subtracted log-affinity entropy and Gini, and the
top-off-target log-affinity gap. It retains the declared parameter-median summaries of the companion study,
with average ranks for ties and eligibility max(pKd)>6 applied before ranking the 200 compounds.
The sweeps use thresholds 5.5–8.0 and baselines 5.0–6.5 in 0.25 steps, and reference off-target indices
1–5. This parameter-median comparison retains the companion convention of assigning zero to undefined
zero-total-weight entropy/Gini configurations; it is therefore a convention-sensitive descriptive analysis.
The fixed-default, floor-baseline controls below have defined weights on the eligible population. Under
those implementations the consensus leaves **37.9% of pairs incomparable**. Gini (0.949) and entropy (0.944)
are close proxies for the checked uniform average-class-rank aggregate; the ratio correlates at 0.817.
These are agreement summaries conditional on this base set. Omitting the ratio reduces incomparability to
12.9%, demonstrating material base-set dependence.

Implementation changes the result. At fixed default parameters, replacing log-affinity entropy by published
association-constant entropy changes Klaeger incomparability from 39.7% to 42.2% when floor-coded entries have
zero weight, and to 65.3% when they retain floor weight. These fixed-default estimates answer a different
question from the primary parameter-median estimate. Published and modified entropies therefore cannot be
treated as interchangeable representatives of one empirical agreement cluster. The four-dataset controls
report both forms wherever the assay is an affinity scale; exponentiating single-concentration percentage
inhibition is not presented as association-constant entropy.

Active-target count is a shared-input proxy: it exactly determines the S-score at the same threshold.
With all four measures, strict discordance is 77.2% for equal-count pairs and 16.2% for widely separated
pairs; withholding the S-score changes those figures to 32.7% and 16.2% (Figure 1). Thus the magnitude of the
separation pattern depends strongly on the discordance panel. The fitted addition of axis position changes
R² from 0.1526 to 0.1534, a descriptive result rather than an estimate of a concept-level component.
Zero-active compounds remain unstable under the declared all-configuration conventions (rank standard
deviation 76.6 versus 32.4), but those conventions include undefined-score assignments; this is evidence for
explicit eligibility and undefinedness handling, not a universal property of low-signal biology.

The current companion candidate is gated exponential entropy: weights are zero at or below the assay floor
x₀ and exp((x−x₀)/T) above it, with an activity-eligibility gate max(x)>τ*. Its four current criteria are
activity eligibility (D1), boundary locality (D2), monotonicity under weak-target addition (D3), and
orientation on equally strong targets (D4). The D3 construction guarantee requires T≥τ*−x₀; T=τ*−x₀ is the
default. Boundary locality means that moving a gate without crossing any value leaves the score unchanged;
it does not promise stable rankings when values are crossed. The historical bounded-gap and distributional-
consistency criteria, and the earlier smooth-hinge candidate, are not current validation claims.

Panel recovery is evaluated separately from these construction properties. On Davis 20-target random panels,
log-affinity entropy has mean conditional Spearman 0.917 but mean eligibility retention only 0.758 and
unconditional pair recovery 0.509. The candidate on the same panels has conditional Spearman 0.854 and pair
recovery 0.484. Pair recovery scores lost eligibility as zero and a new tie as one-half, on pairs strictly
ordered and eligible on the full panel. High conditional rank recovery alone therefore overstates how much
of the original comparison remains usable.

A matched-platform control reuses 46 CID-resolved Davis/Klaeger compounds, 201 shared curated genes, and the
40 compounds eligible on both platforms. Across 780 pairs, the incomparable-set Jaccard is 0.294 (matched-
compound bootstrap 95% CI [0.205, 0.383]). Direction agreement is 0.954 among the 348 pairs strictly comparable
on both platforms. This directly measures pair-set overlap, without inferring that either assay reveals the
true concept order. Full results are in `kinase_controls.txt`, `kinase_platform.txt` and `protocol_controls.txt`.

### 4.2 Epigenetic clocks (biological age)
Concept: biological age. Measures: DNA-methylation clocks of three declared *types*, which the protocol's
step 1 keeps separate rather than pooling. Position clocks report where a sample is on the ageing
trajectory, deviation clocks report departure from it, and rate clocks report speed along it (1,385 blood
samples, GSE40279 and GSE87571). Rate clocks such as DunedinPACE (Belsky et al., 2022) are one declared
type, and independent comparisons confirm that clocks disagree and differ in predictive value
(Mavrommatis et al., 2025). The single-linear clocks are affine-related to one another (R² 0.992 to
0.995, over the position clocks Horvath and Hannum and the deviation clock PhenoAge), and the rate clock
relates to them through a signed logarithm (R² 0.93 to 0.95), a Weber-Fechner-like relation that is a
descriptive observation, not a scale-type proof (Section 2). Every clock tested against cell-type
composition fails an ordinal-stability check (absolute correlation up to 0.372, over Horvath, Hannum,
PhenoAge and GrimAge, with DunedinPACE not tested), and the four clocks with a single published
coefficient vector are near-orthogonal to each other (pairwise cosine between −0.17 and 0.11), so they
agree on age while sharing almost no CpGs. GrimAge is a two-stage composite with no single coefficient
vector, so it cannot enter that comparison. The
external anchor here, chronological age, is the training target of the position clocks, so for them the
orientation check is confirmatory rather than independent.

The clock idea is not specific to blood or to methylation. On a rat retinal bulk RNA-seq aging series
(GSE314970, 85 samples over 6 to 27 months, Shavlakadze et al., 2026), five competing transcriptomic aging
measures disagree sharply: the SenMayo and Fridman senescence signatures (Saul et al., 2022, Fridman and
Tainsky, 2008), a Hallmark inflammatory inflammaging proxy, and a fitted elastic-net clock. Against true age
they range from Spearman 0.98 for the fitted clock down to 0.13 for the Fridman senescence-up signature. The
curated senescence signatures are weak age predictors on their own. The five measures form a senescence
and inflammaging agreement cluster plus a distinct fitted-clock direction, with the Fridman senescence-down signature
anti-correlated with the rest of the senescence block. The mean-rank consensus tracks age (Spearman 0.755) and
mean age rises monotonically across consensus tertiles (9.7, 18.3, 22.7 months), and discordance falls with
separation on the consensus axis as elsewhere, though even well-separated pairs still disagree about half the
time. This is a single-cohort, rat, single-sex extension with the clock trained in-distribution, so it is
qualitative, on par with the inflammaging clocks. It shows the clock agreement-cluster structure and the near-tie pattern
recur in a new tissue and a new measurement modality. Full analysis is in the companion clock repository
`epigenetic-clock-desiderata`, subdirectory `eye_aging/`.

### 4.3 Serotonin-receptor selectivity
Concept: concentration of observed binding across a declared serotonin-receptor panel. This is distinct from
selectivity for a designated 5-HT2A target over 5-HT2B. The ChEMBL-derived companion subset contains 297
ligands with both receptors measured and at least four receptors tested, but their observed panels range
from 4 to 12 of 13 available receptor columns. After requiring an observed target above pKi 6, 285 remain.
We compare the observed-target S-score, log-affinity entropy and Gini, and a top-two log-affinity gap. Untested
receptors remain unknown; floor imputation is a sensitivity scenario. The companion's quotient of its top two
pKi values is a separate implementation, not the kinase log-affinity gap.

On the 285-ligand subset the observed-only consensus leaves 75.4% of pairs incomparable, versus 56.1% under
floor imputation. Entropy's relation to measured-panel size is strong (absolute Spearman 0.943). The broader
949-ligand, at-least-four-tested-receptors set retains 898 eligible compounds; on its largest identical
observed receptor panel, 332 eligible compounds remain and entropy/Gini correlate at 0.996. The apparent
entropy-outlier pattern in variable panels therefore cannot be attributed confidently to a biological
measure family. The observed-only mean-rank consensus correlates with observed active-target count at 0.740
on the broader eligible set; that count shares assay inputs and coverage with the scores. A separate
36-compound curated psychedelic panel in the companion repository is descriptive context only. Coverage,
implementation and selection preclude a general cross-domain transfer claim from this instance.

### 4.4 Methylation exposure biomarkers
Concept: environmental exposure or disease signature. Measures: competing smoking-exposure signatures (AHRR,
Breitling et al., 2011, EpiSmoke, Shenker et al., 2013, Joehanes, Joehanes et al., 2016, and EpiTob,
Chamberlain et al., 2022) across two cohorts (GSE50660, Tsaprouni et al., 2014, n=464, and GSE42861,
Liu et al., 2013, n=689). The signatures agree almost
perfectly on never-smokers (96%) and well on current smokers (82%), and diverge on former smokers (71% pooled,
78% and 66% per cohort). The divergence reproduces in both cohorts separately, indicating stable structure conditional on the analysis,
but it is artefactual rather than constitutive: it is a threshold effect at the classification cutoff, and
which group agrees worst changes when the cutoff changes, so it is not evidence of a concept-level
"consequential middle" (Section 5). Freezing never-smoker thresholds from GSE50660 and applying them to
GSE42861 changes never/former/current agreement to 0.819/0.649/0.885; cohort-specific threshold results
are therefore descriptive, not a frozen-rule validation. A cross-tissue arm (fetal alcohol spectrum disorder) shows why the
reliability gate is needed. Buccal-derived signatures (Portales-Casamar et al., 2016, Lussier et al., 2018)
scored on blood transfer poorly (AUC 0.68 to 0.79), and the blood-native episignature
(van der Laan et al., 2025, AUC 0.93 to 0.96) is non-independent of its evaluation cohort (shared laboratory
and tissue, confirmed by the depositing author), so its apparent performance may be optimistic and does not supply independent validation.

### 4.5 TP53 variant effect Concept: the functional impact of a TP53 missense variant. Measures: eight
Kato/Ishioka yeast transactivation readouts (Kato et al., 2003), seven computational pathogenicity predictors,
namely AlphaMissense (Cheng et al., 2023), EVE (Frazer et al., 2021), ESM1b (Brandes et al., 2023), REVEL
(Ioannidis et al., 2016), CADD (Kircher et al., 2014, with the v1.6 scores of Rentzsch et al., 2021),
PrimateAI (Sundaram et al., 2018) and BayesDel (Feng, 2017), the last four taken from dbNSFP v4.1a (Liu et
al., 2020), and mammalian proliferation screens (Giacomelli et al., 2018, Kotler et al., 2018, via ProteinGym,
Notin et al., 2023), with ClinVar (Landrum et al., 2018) as the external clinical anchor and the ClinGen TP53
specifications (Fortuno et al., 2021) for the field's functional rules (core panel 2,314 variants). Step 1 is
necessary here. The eight transactivation readouts and AlphaMissense are correctly oriented against ClinVar
(AUC 0.96 to 0.99), while all four proliferation screens are mis-oriented (AUC 0.016 to 0.064) under the naive
"higher fitness is worse" convention, so they must be flipped before use. The fifteen core-panel measures form
two agreement clusters by measurement basis, the eight transactivation readouts (within-cluster mean 0.78) and the seven
predictors (mean 0.76), agreeing less across the two blocks (mean 0.50) than within either. The proliferation
screens are held out of the agreement-cluster analysis, and the three Giacomelli conditions serve instead as the
independent mammalian platform in the reproducibility test below. The consensus poset over the nine base
measures (the eight readouts plus AlphaMissense) is 56.6% incomparable. It contains 2,313 equivalence classes
among 2,314 variants; the tied pair receives one shared class score. With four dispersed insertion-Gibbs
chains passing the numerical criteria, the field's median-of-eight-promoters rule correlates with the
uniform average-class-rank aggregate at Spearman 0.955, and computational predictors range from 0.572 to
0.729. The aggregate discriminates the retained ClinVar labels at AUC 0.996; this is a descriptive benchmark
on the selected labelled variants, not independent proof of biological validity. The earlier same-start
adjacent-swap estimates were insufficiently mixed and are superseded.

Strict discordance among the eight promoters falls from 0.924 to 0.045 across consensus-separation quintiles,
and from 0.606 to 0.256 across AlphaMissense-separation quintiles (Figure 1). The fitted position coefficients
are +0.362 and +0.119 respectively. AlphaMissense is outside this promoter discordance panel, but its residual
association does not identify a constitutive component or prove a partial-loss-of-function mechanism.
Two experimental *platforms*, yeast transactivation and
mammalian proliferation, agree on the pooled extremes (+0.63) and not in the middle (−0.05), but that
contrast is carried by one end: tertile by tertile the platforms agree on deleterious variants (+0.56) and no
better than chance on functional ones (+0.01) or in the middle (−0.05). The shared pattern is asymmetric and consistent with resolution concentrated at one end; it does not
establish a symmetric skeleton with an undecided middle. The gap
from the middle does survive the fair comparison (mean of the two end tertiles minus the middle +0.34,
bootstrap 95% CI [+0.25, +0.42]), so there is reproducible cross-platform structure, but it sits at the
deleterious end alone. Finally, all
seven predictors score 0.90 to 0.99 against curated ClinVar but only 0.79 to 0.89 against the independent
functional truth on the variants of uncertain significance where prediction is actually used, so benchmark performance depends on the selected population and endpoint (Section 5). These comparisons
change both variant spectrum and reference labels; they do not isolate training leakage as the cause. A hotspot-codon holdout showed no inflation
for any predictor, an honest null.

### 4.6 Inflammatory burden Concept: systemic inflammatory burden, quantified three ways of differing
strength and implementation fidelity. **Clinical indices (stable correlation structure):** seven routine indices (NLR, PLR, MLR, SII, SIRI, CRP, CRP-to-albumin)
on two independent NHANES cycles, anchor all-cause mortality. Two agreement clusters emerge by measurement basis,
leukocyte-count ratios and CRP-protein, weakly correlated with each other (between-cluster mean 0.19), and CRP
and its albumin ratio are highly correlated (Spearman 0.998), but their rank vectors are not identical.
Omitting either is therefore a sensitivity analysis, not removal of an exact duplicate. The consensus discriminates mortality at
the extremes (AUC 0.69) but not the middle (0.51), and its inter-index correlation structure replicates across
the two cycles at Spearman **0.995** (participant-bootstrap 95% CI [0.969, 0.999]).
These cohort summaries use 6,197 and 5,882 complete-index participants; the mortality-linked subsets contain
5,327 and 5,104. The 600-person poset subsample is drawn from the former mortality-linked population.
Its closest single index to the checked aggregate is SII (0.844), while the equal-weight mean-rank score
agrees at 0.983. The correlation structure can reproduce because the indices share inputs and algebra;
replication does not distinguish that mechanism from a concept-level partial order. Disagreement falls steeply
with separation on the consensus axis (0.984 for near-ties to 0.274 for well-separated pairs) but is flat on
the external-to-formula axis examined here, chronological age (0.735 to 0.694, R² 0.0013, Figure 1). Age is too
weak a proxy for inflammatory burden to settle which reading is right, so we report both. **Sepsis gene-set proxies
(strong disagreement, orientation-sensitive reproducibility):** five blood gene sets drawn from published signatures (the Hallmark inflammatory and
interferon-γ sets, Liberzon et al., 2015, the SRS7 Sepsis Response Signature, Davenport et al., 2016, the
SRSq19 SepstratifieR signature, Cano-Gamez et al., 2022, and the MARS8 endotype signature, Scicluna et al.,
2017) on GSE65682, scored on the 479 of 802 samples carrying 28-day mortality, which is also the anchor. These scores are mean-z gene-set proxies, not reconstructions of the published SRSq or MARS classifiers.
Gene-wise centers/scales and score directions are fitted on discovery and frozen on validation. With that rule the detailed correlation
structure replicates at Spearman 0.976, versus the earlier 0.32 obtained by fitting directions independently
in each cohort. The participant-bootstrap 95% interval is [0.684, 0.988], reflecting uncertainty in only ten
indicator-pair correlations. The full frozen-sign consensus is 87.1% incomparable, with the inflammatory
proxy closest to its checked aggregate (0.613). Hallmark and SRS proxy blocks remain weakly correlated
(between-cluster mean +0.05). Four of five marginal mortality directions match; MARS8 flips and its pooled
raw AUC is 0.505. The sign instability concerns outcome association, while the fixed-definition structure
is considerably more stable. Neither observation establishes concept-level incomparability. **Inflammaging clocks (thin):** only two composite clocks are publicly reconstructable on
shared samples, SImAge (Kalyakulina et al., 2023) and a recomputed ipAGE (Yusipov et al., 2022). Two further
inflammaging clocks, iAge (Sayed et al., 2021) and IMM-AGE (Alpert et al., 2019), are not reconstructable
here, so this is a qualitative check, not a full instance. Sepsis *endotype labels* (SRS versus MARS) are
partitions rather than scalar orders, so the right comparison for them is cluster agreement rather than
average rank. Only the signature scores are analysed with the protocol.

### 4.7 Cardiovascular risk
Concept: ten-year risk of a cardiovascular event. Measures: three published scores in routine clinical use,
the ACC/AHA Pooled Cohort Equations (Goff et al., 2014), and the Framingham general-CVD laboratory and office
models (D'Agostino et al., 2008), the last substituting body-mass index for the lipid panel. Cohorts: NHANES
2015-16 (n=2,994) and 2017-18 (n=2,913), adults aged 40 to 79 with complete risk factors. The anchor is
cause-specific cardiovascular death from the NHANES Linked Mortality File, with all-cause death reported
alongside because the cardiovascular event count is small (30 and 10). All three are correctly oriented
against both anchors (cardiovascular-death AUC 0.755 to 0.815).

This domain behaves unlike the others, and the contrast is the reason for including it. The three scores do
not split into distinct agreement clusters: they agree on ranking almost perfectly (Spearman 0.932 to 0.946 in both cycles) and
that structure is the most stable in this work, reproducing across cycles to a mean absolute difference of
0.004. The consensus poset still leaves 15.5% and 15.7% of pairs incomparable, the lowest of any instance in
Table 1. The scores differ in their numerical risks: median ten-year values are
7.5% under the Pooled Cohort Equations and 15.3% under the Framingham office model. Their endpoints also
differ: PCE predicts hard ASCVD whereas the Framingham models predict broader general CVD (Goff et al., 2014;
D'Agostino et al., 2008). Applying a common numerical cutoff therefore mixes endpoint and formula changes
and cannot establish a calibration error or interchangeable treatment-eligibility decisions.

At a shared cutoff of 0.075, 26.2% and 26.6% of adults fall on different sides under at least one score;
at 0.10 the figures are 25.9% and 26.9%. Among the two Framingham models, which share an endpoint, the 0.075
classification discrepancies are 10.0% and 10.3%, and the 0.10 discrepancies are 11.1% and 11.6%. These are
classification sensitivity results, not estimated changes in clinical management. Mortality during observed
follow-up supplies an orientation check; it does not validate ten-year incident-event calibration. Strict
indicator-pair discordance declines from 0.444 to 0.000 along the consensus axis, a descriptive result linked
to that same score panel. This example illustrates why the protocol first requires matching outcomes and
measurement types. Analysis is in `cvdrisk/`.

## 5. What the demonstrations establish

With the controls applied, the demonstrations support the following conditional findings.

1. **Comparative structure depends on implementation and eligibility.** Under the primary parameter-median
   kinase implementation, 37.9% of pairs are incomparable; fixed-default implementation and floor scenarios
   give different estimates. Serotonin panel size strongly predicts entropy. These effects must be resolved
   before assigning biological meaning to agreement clusters.
2. **Separation associations depend on the axis and the discordance panel.** In TP53, strict discordance
   falls from 0.924 to 0.045 on the consensus axis and from 0.606 to 0.256 on the external-to-panel
   AlphaMissense axis. For inflammatory indices the consensus decline is 0.984 to 0.274, whereas age gives
   0.735 to 0.694. In kinase, removing the S-score from an active-count-axis comparison changes equal-count
   discordance from 77.2% to 32.7%, with wide-separation discordance remaining 16.2%. These are descriptive
   sensitivity results, not estimates or bounds on a constitutive component. The TP53 position residual
   warrants further biological investigation without identifying its cause.

![Strict reversals against separation, with external-proxy and shared-input controls.](figures/fig1_near_tie.png)

**Figure 1.** Strict inter-measure reversals versus pair separation. TP53 and inflammatory-index points are
separation-quantile bins; consensus axes share the scores whose disagreement is plotted, whereas AlphaMissense
and age are outside those discordance panels. Age is a weak inflammatory-burden proxy. Kinase bins are actual
active-count differences (0, 1–2, 3–5, 6–15, ≥16), not separation quantiles. Active count defines the default
S-score, so both the four-measure and S-score-withheld results are shown. Ties versus strict comparisons are
reported separately in the stored outputs and are excluded from the reversal curves. Generated by
`paper/make_figures.py` from those outputs.

3. **Replication establishes stability conditional on the comparison.** NHANES correlation-structure
   reproducibility is 0.995 with bootstrap interval [0.969, 0.999]. For sepsis proxies, freezing discovery gene-wise preprocessing and
   directions gives 0.976 with interval [0.684, 0.988]; fitting signs independently had produced an artificial
   impression of weak replication. In contrast, matched-compound kinase incomparable-set overlap is only
   Jaccard 0.294 despite high direction agreement on pairs comparable on both assays. Matrix-vector
   correlations, scalar-consensus agreement and direct pair-set overlap measure different things.
4. **Aggregate benchmarks require convergence and an explicit convention.** The checked uniform class-rank
   aggregate correlates with the TP53 median-of-eight rule at 0.955, with kinase Gini at 0.949, and with the
   inflammatory-index mean-rank score at 0.983. These support useful descriptive approximations, not a claim
   that an aggregate identifies the correct biology. The same-start sampler had substantially overstated
   some agreements: SII's aggregate correlation changes from 0.983 to 0.844 after adequate mixing.
5. **Recovery must include the objects that disappear.** Davis 20-target log-entropy recovery has conditional
   Spearman 0.917 but eligibility retention 0.758 and pair recovery 0.509. Coverage, conditional accuracy,
   unconditional pair recovery and their panel uncertainty answer different questions.
6. **Orientation checks have differing independence.** Concept-positive scores have the expected signs
   against the declared proxies, but active count shares inputs with selectivity, age trains position clocks,
   and discovery mortality fits the sepsis signs. Smoking metadata and clinical labels supply different
   kinds of evidence, with their own selection and curation limits. Pooling opposite extremes also creates a
   range-restriction advantage over a middle tertile; the per-tertile checks do not establish a universal
   symmetric extremes-versus-middle pattern.

The framework does not identify constitutive disagreement from replication alone. Its contribution is an
explicit comparison and control procedure, with evidence for stability kept separate from attribution to a
latent concept. Prediction benchmarks remain appropriate for a declared endpoint; comparing curated ClinVar
labels with functional classifications on VUS changes both the reference and the variant population, so the
observed AUC gap does not isolate training overlap or prove benchmark inflation. The formalism is borrowed
from measurement and partial-order theory; the empirical contribution is the controlled application.

## 6. Discussion

**Reporting practice.** The protocol implies concrete changes to how competing measures are reported. Declare
each measure's type and scale, and check its orientation against a declared anchor, with dependence recorded, before any comparison
(step 1). Report the reproducible incomparable set and the per-measure or per-cluster criterion verdicts rather than a single
ranking (step 6). Validate a measure in its native domain, tissue, or panel rather than assuming transfer (the
reliability gate, Sections 4.4 and 4.6). When a single number is unavoidable, report the
chosen average-rank aggregate, its numerical uncertainty and each measure's agreement with it (step 7).

**The conceptual claim and its limits.** Partial orders offer a useful representation of competing
biological quantifications, conditional on substantive admissibility assumptions. The observable findings
here are disagreements among selected implementations and their sensitivity to populations, panels,
missingness and computational choices. Replication strengthens confidence that a difference is stable; it
does not establish that the latent concept itself leaves the pair unordered. The strongest defensible
output is therefore a conditional comparative structure with explicit uncertainty and unresolved causes.

**Limitations.** Four are substantive. First, the base set is a choice. The consensus order, the incomparable
set, and the canonical aggregate are all defined relative to the set of measures admitted, and there is no
closed enumeration of "all measures" of a biological concept. We mitigate this with a stated selection rule
(measures in routine use, computable on the data, distinct in the order they induce), with robustness checks,
and by reporting the base set explicitly, but the choice itself remains.
Second, replication estimands and evidence differ. Population correlation structures can be stable
while individual pair relations transfer less well, and the sepsis proxy interval is wide despite a high
fixed-sign point estimate. Official SRSq/MARS classifiers and an independent prospective validation of the
proxy scores are not reconstructed here. Third, anchor independence varies: active-target count shares inputs with the selectivity scores,
chronological age is a training target of position clocks, and sepsis mortality fits discovery signs.
These checks are confirmatory or shared-input checks rather than independent validation. Fourth, the framework assumes a
concept-level order exists and that reasonable measures respect it. This is one substantive assumption, argued
rather than proved. It also restricts scope to scalar orders, so categorical *endotype* systems fall
outside the protocol (Section 4.6).

**Scope and relation to prior work.** The formalism is entirely borrowed. It is representational
measurement theory (Krantz et al., 1971), the partial-order theory of multi-indicator systems (Patil and
Taillie 2004, Brüggemann and Patil, 2011), Szpilrajn's extension theorem, and the average-rank canonical
extension (De Loof et al., 2008). The clearest prior application of measurement theory in biology (Houle et
al., 2011) treats single measurements and scale types rather than competing measures of one concept as
extensions of a shared partial order. The partial-order approach has precedent in environmental chemistry
(Halfon and Reggiani, 1986) and in social indicators (Fattore, 2016), and it contrasts with
rank-aggregation methods that force a consensus (Dwork et al., 2001, Kolde et al., 2012) rather than
preserve incomparability. An early biological example is the concordance among breast-cancer
gene-expression predictors (Fan et al., 2006), which share few components yet agree where the concept is
determined. Our contribution is the protocol, its implementation, eligibility and replication controls,
and the cross-domain demonstration on public data. For a working analyst, it gives a
reusable procedure for the common situation where several defensible measures of one latent concept do not
agree.

---

## Declarations

**Declaration of competing interest.** The author declares no competing financial or personal interests.

**Funding.** This research received no external funding.

**Declaration of generative AI in scientific writing.** During the preparation of this work the author
used Claude (Anthropic) to assist with data acquisition, analysis scripting, and drafting. The author
reviewed and edited all output and takes full responsibility for the accuracy and integrity of the work.

**Data and code availability.** All analyses use publicly available data and are fully reproducible from
the project repository. For each domain the repository contains a README documenting data provenance and
the measure-selection rule, the analysis scripts, and their stored outputs. *Repository URL and an archived
DOI (Zenodo) to be added on submission.*

**CRediT author contributions.** Polina Vinogradova: conceptualization, methodology, software, formal
analysis, investigation, writing – original draft, writing – review and editing.

---

## Appendix A — formal justification

This appendix states which established result justifies which step of the protocol, and fixes the object
so the steps are well-defined. Every result is classical or definitional.

**Definition 1 (the object).** Fix a latent concept κ and a finite object set Ω. The empirical structure
is (Ω, ≽), where ≽ is a preorder (reflexive and transitive) that is not assumed total. Write a ~ b when
both a ≽ b and b ≽ a, write a ≻ b for strict order, and write a ∥ b when neither a ≽ b nor b ≽ a
(a and b are incomparable).

**Definition 2 (admissible measure).** A measure is a map f : Ω → ℝ. It is admissible if a ≽ b implies
f(a) ≥ f(b) and a ≻ b implies f(a) > f(b). This is the order-preserving representation condition of
representational measurement theory. Write 𝓕(≽) for the admissible measures.

**Definition 3 (families, weak extensions, linear extensions).** A measure f induces a total preorder on Ω:
write a ⪰_f b iff f(a) ≥ f(b), and write ≻_f for its strict part. Two measures are in the same *family* when
they induce the same total preorder. A total preorder on Ω is a *weak extension* of ≽ when it preserves the
concept's weak and strict comparisons alike: a ≽ b implies a ⪰_f b, and a ≻ b implies a ≻_f b. A *linear extension* in the standard order-theoretic sense is defined on the quotient Ω/~,
where each equivalence class becomes one object. It preserves class comparisons and introduces no ties
between distinct classes. On finite Ω, f ↦ ⪰_f maps 𝓕(≽) onto the weak extensions of ≽, and the
linear extensions order the equivalence classes without further ties. The distinction matters because an admissible measure may
assign two incomparable objects the same value, whereas the aggregate below averages over linear extensions
only. Throughout, ≽ and ≻ are the concept's order and ⪰_f and ≻_f the order a measure induces. This is the
multiple-indicators-to-partial-order setup of Patil and Taillie (2004).

**Consensus (how ≽ is estimated).** Given the set 𝓜 of measures in use, define a ≽𝓜 b iff f(a) ≥ f(b) for
every f in 𝓜. If every reasonable measure is admissible for the true ≽, then ≽ ⊆ ≽𝓜, so the consensus
over-declares comparability. The comparable pairs of ≽𝓜 are an upper bound on the truly comparable pairs,
and the incomparable pairs of ≽𝓜 are a lower bound on the truly incomparable pairs. Reproducibility across independent cohorts (protocol step 3) tests stability of this estimate under
resampling populations, conditional on fixed implementations. It does not establish admissibility or
exclude systematic artifacts. The independence of an anchor (step 1) must be examined separately.

**Lemma 1 (disagreement lies in the incomparable set).** If f, g in 𝓕(≽) strictly disagree on a pair,
meaning f(a) > f(b) while g(a) < g(b), then a ∥ b. *Proof.* If a ≽ b then admissibility of g gives
g(a) ≥ g(b), a contradiction. If b ≽ a then admissibility of f gives f(b) ≥ f(a), a contradiction. So
a ∥ b. This is the contrapositive of admissibility. For the estimated order ≽𝓜 it holds by construction,
and its role is only to justify step 2.

**Lemma 2 (families are ordinal classes).** If f, g in 𝓕(≽) are in the same family, there is a strictly
increasing φ with g = φ ∘ f, and conversely any strictly increasing φ gives a same-family measure. *Proof.*
Define φ(f(a)) := g(a). It is well-defined and strictly increasing because f and g induce the same order. So a
family is an equivalence class under strictly monotone reparameterisation, an ordinal class.
Linear-relatedness on a dataset therefore does not establish a shared interval scale, and a non-affine but
strictly increasing relation between two measures keeps them in the same family. Family membership is
diagnosed by rank agreement, not by the regression form relating two measures.

**Proposition (non-uniqueness).** On finite (Ω, ≽): 𝓕(≽) is non-empty, the families are in bijection with
the weak extensions of ≽, and there is a unique family iff ≽ is total. *Proof.* Apply Szpilrajn's extension theorem to the finite quotient Ω/~. Numbering its classes in the
resulting order and assigning each object's class number gives an admissible representation on Ω. Conversely a weak extension is realised by numbering its levels. A finite
poset has a unique weak extension iff it is total. If a ∥ b, both a ≻ b and b ≻ a extensions
exist. So when incomparabilities exist there are at least two order-distinct admissible measures, and no
requirement that only entails order-preservation can single one out. This justifies step 6.

**The aggregate convention.** Form the quotient of ≽𝓜 by mutual comparability (identical indicator
vectors on the observed objects). On this finite poset, assign every linear extension equal probability,
and average each class's rank over those extensions. Copy the resulting score to every object in the class.
This preserves strict comparisons and ties of the consensus preorder. It can order incomparable classes;
the uniform distribution over extensions is an additional explicit convention. Mutual-rank probabilities
retain pairwise uncertainty. Exact enumeration is used for small quotients. For larger ones an insertion
Gibbs step selects a class uniformly, removes it, and reinserts it uniformly into the slots allowed by the
other classes. Conditional on their order those slots have equal probability under the uniform extension
law, so this update satisfies detailed balance. Adjacent swaps of incomparable classes are reachable, and
self transitions give aperiodicity. Those stationary-distribution properties do not guarantee adequate
finite-run mixing; the dispersed-start, split-R-hat and precision controls in Section 3 are required.
Numerical acceptance is empirical and does not remove uncertainty about the base set or the aggregation
convention.

## References

*Elsevier name-date (Harvard) style; DOIs verified against primary/publisher records. `[UNVERIFIED]` tags mark
residual gaps to confirm at proof. Journal house style (including semicolon-separated in-text citation lists)
is applied by the journal at proof.*

Alpert, A., Pickman, Y., Leipold, M., Rosenberg-Hasson, Y., Ji, X., Gaujoux, R., Rabani, H., Starosvetsky, E.,
Kveler, K., Schaffert, S., Furman, D., Caspi, O., Rosenschein, U., Khatri, P., Dekker, C.L., Maecker, H.T.,
Davis, M.M., Shen-Orr, S.S., 2019. A clinically meaningful metric of immune age derived from high-dimensional
longitudinal monitoring. Nature Medicine 25, 487–495. https://doi.org/10.1038/s41591-019-0381-y

Anastassiadis, T., Deacon, S.W., Devarajan, K., Ma, H., Peterson, J.R., 2011. Comprehensive assay of kinase
catalytic activity reveals features of kinase inhibitor selectivity. Nature Biotechnology 29, 1039–1045.
https://doi.org/10.1038/nbt.2017

Arrow, K.J., 1963. Social Choice and Individual Values, 2nd ed. Wiley, New York.

Belsky, D.W., Caspi, A., Corcoran, D.L., Sugden, K., Poulton, R., Arseneault, L., Baccarelli, A., Chamarti,
K., Gao, X., Hannon, E., Harrington, H.L., Houts, R., Kothari, M., Kwon, D., Mill, J., Schwartz, J., Vokonas,
P., Wang, C., Williams, B.S., Moffitt, T.E., 2022. DunedinPACE, a DNA methylation biomarker of the pace of
aging. eLife 11, e73420. https://doi.org/10.7554/eLife.73420

Bland, J.M., Altman, D.G., 1995. Comparing methods of measurement: why plotting difference against standard
method is misleading. Lancet 346, 1085–1087. https://doi.org/10.1016/S0140-6736(95)91748-9

Bosc, N., Meyer, C., Bonnet, P., 2017. The use of novel selectivity metrics in kinase research. BMC
Bioinformatics 18, 17. https://doi.org/10.1186/s12859-016-1413-y

Brandes, N., Goldman, G., Wang, C.H., Ye, C.J., Ntranos, V., 2023. Genome-wide prediction of disease variant
effects with a deep protein language model. Nature Genetics 55, 1512–1522.
https://doi.org/10.1038/s41588-023-01465-0

Breitling, L.P., Yang, R., Korn, B., Burwinkel, B., Brenner, H., 2011. Tobacco-smoking-related differential
DNA methylation: 27K discovery and replication. American Journal of Human Genetics 88, 450–457.
https://doi.org/10.1016/j.ajhg.2011.03.003

Brightwell, G., Winkler, P., 1991. Counting linear extensions. Order 8, 225–242.
https://doi.org/10.1007/BF00383444

Brüggemann, R., Patil, G.P., 2011. Ranking and Prioritization for Multi-indicator Systems: Introduction to
Partial Order Applications. Springer, New York.

Bubley, R., Dyer, M., 1999. Faster random generation of linear extensions. Discrete Mathematics 201, 81–88.
https://doi.org/10.1016/S0012-365X(98)00333-1

Campbell, D.T., Fiske, D.W., 1959. Convergent and discriminant validation by the multitrait-multimethod
matrix. Psychological Bulletin 56, 81–105. https://doi.org/10.1037/h0046016

Cano-Gamez, E., Burnham, K.L., Goh, C., Allcock, A., Malick, Z.H., Overend, L., Kwok, A., Smith, D.A.,
Peters-Sengers, H., Antcliffe, D., McKechnie, S., Scicluna, B.P., van der Poll, T., Gordon, A.C., Hinds, C.J.,
Davenport, E.E., Knight, J.C., 2022. An immune dysfunction score for stratification of patients with acute
infection based on whole-blood gene expression. Science Translational Medicine 14, eabq4433.
https://doi.org/10.1126/scitranslmed.abq4433

Chamberlain, J.D., Nusslé, S., Chapatte, L., Kinnaer, C., Petrovic, D., Pradervand, S., Bochud, M., Harris,
S.E., Corley, J., Cox, S.R., Gonseth Nusslé, S., 2022. Blood DNA methylation signatures of lifestyle
exposures: tobacco and alcohol consumption. Clinical Epigenetics 14, 155.
https://doi.org/10.1186/s13148-022-01376-7

Cheng, J., Novati, G., Pan, J., Bycroft, C., Žemgulytė, A., Applebaum, T., Pritzel, A., Wong, L.H., Zielinski,
M., Sargeant, T., Schneider, R.G., Senior, A.W., Jumper, J., Hassabis, D., Kohli, P., Avsec, Ž., 2023.
Accurate proteome-wide missense variant effect prediction with AlphaMissense. Science 381, eadg7492.
https://doi.org/10.1126/science.adg7492

D'Agostino, R.B. Sr., Vasan, R.S., Pencina, M.J., Wolf, P.A., Cobain, M., Massaro, J.M., Kannel, W.B., 2008.
General cardiovascular risk profile for use in primary care: the Framingham Heart Study. Circulation 117,
743–753. https://doi.org/10.1161/CIRCULATIONAHA.107.699579

Davenport, E.E., Burnham, K.L., Radhakrishnan, J., Humburg, P., Hutton, P., Mills, T.C., Rautanen, A., Gordon,
A.C., Garrard, C., Hill, A.V.S., Hinds, C.J., Knight, J.C., 2016. Genomic landscape of the individual host
response and outcomes in sepsis: a prospective cohort study. Lancet Respiratory Medicine 4, 259–271.
https://doi.org/10.1016/S2213-2600(16)00046-1

Davis, M.I., Hunt, J.P., Herrgard, S., Ciceri, P., Wodicka, L.M., Pallares, G., Hocker, M., Treiber, D.K.,
Zarrinkar, P.P., 2011. Comprehensive analysis of kinase inhibitor selectivity. Nature Biotechnology 29,
1046–1051. https://doi.org/10.1038/nbt.1990

De Loof, K., De Baets, B., De Meyer, H., Brüggemann, R., 2008. A hitchhiker's guide to poset ranking.
Combinatorial Chemistry & High Throughput Screening 11, 734–744. https://doi.org/10.2174/138620708786306032

De Loof, K., De Baets, B., De Meyer, H., 2011. Approximation of average ranks in posets. MATCH Communications
in Mathematical and in Computer Chemistry 66, 219–229.

Dwork, C., Kumar, R., Naor, M., Sivakumar, D., 2001. Rank aggregation methods for the web. In: Proceedings of
the 10th International Conference on World Wide Web (WWW '01). ACM, New York, pp. 613–622.
https://doi.org/10.1145/371920.372165

Fan, C., Oh, D.S., Wessels, L., Weigelt, B., Nuyten, D.S.A., Nobel, A.B., van't Veer, L.J., Perou, C.M., 2006.
Concordance among gene-expression-based predictors for breast cancer. New England Journal of Medicine 355,
560–569. https://doi.org/10.1056/NEJMoa052933

Fattore, M., 2016. Partially ordered sets and the measurement of multidimensional ordinal deprivation. Social
Indicators Research 128, 835–858. https://doi.org/10.1007/s11205-015-1059-6

Feng, B.J., 2017. PERCH: a unified framework for disease gene prioritization. Human Mutation 38, 243–251.
https://doi.org/10.1002/humu.23158

Fortuno, C., Lee, K., Olivier, M., Pesaran, T., Mai, P.L., de Andrade, K.C., Attardi, L.D., Crowley, S.,
Evans, D.G., Feng, B.J., Foreman, A.K.M., Frone, M.N., Huether, R., James, P.A., McGoldrick, K., Mester, J.,
Seifert, B.A., Slavin, T.P., Witkowski, L., Zhang, L., Plon, S.E., Spurdle, A.B., Savage, S.A., 2021.
Specifications of the ACMG/AMP variant interpretation guidelines for germline TP53 variants. Human Mutation
42, 223–236. https://doi.org/10.1002/humu.24152

Frazer, J., Notin, P., Dias, M., Gomez, A., Min, J.K., Brock, K., Gal, Y., Marks, D.S., 2021. Disease variant
prediction with deep generative models of evolutionary data. Nature 599, 91–95.
https://doi.org/10.1038/s41586-021-04043-8

Fridman, A.L., Tainsky, M.A., 2008. Critical pathways in cellular senescence and immortalization revealed by
gene expression profiling. Oncogene 27, 5975–5987. https://doi.org/10.1038/onc.2008.213

Giacomelli, A.O., Yang, X., Lintner, R.E., McFarland, J.M., Duby, M., Kim, J., Howard, T.P., Takeda, D.Y., Ly,
S.H., Kim, E., Gannon, H.S., Hurhula, B., Sharpe, T., Goodale, A., Fritchman, B., Steelman, S., Vazquez, F.,
Tsherniak, A., Aguirre, A.J., Doench, J.G., Piccioni, F., Roberts, C.W.M., Meyerson, M., Getz, G.,
Johannessen, C.M., Root, D.E., Hahn, W.C., 2018. Mutational processes shape the landscape of TP53 mutations in
human cancer. Nature Genetics 50, 1381–1387. https://doi.org/10.1038/s41588-018-0204-y

Goff, D.C. Jr., Lloyd-Jones, D.M., Bennett, G., Coady, S., D'Agostino, R.B., Gibbons, R., Greenland, P.,
Lackland, D.T., Levy, D., O'Donnell, C.J., Robinson, J.G., Schwartz, J.S., Shero, S.T., Smith, S.C. Jr.,
Sorlie, P., Stone, N.J., Wilson, P.W.F., 2014. 2013 ACC/AHA guideline on the assessment of cardiovascular
risk. Circulation 129 (25 Suppl 2). https://doi.org/10.1161/01.cir.0000437741.48606.98

Halfon, E., Reggiani, M.G., 1986. On ranking chemicals for environmental hazard. Environmental Science &
Technology 20, 1173–1179. https://doi.org/10.1021/es00154a010

Houle, D., Pélabon, C., Wagner, G.P., Hansen, T.F., 2011. Measurement and meaning in biology. Quarterly Review
of Biology 86, 3–34. https://doi.org/10.1086/658408

Ioannidis, N.M., Rothstein, J.H., Pejaver, V., Middha, S., McDonnell, S.K., Baheti, S., Musolf, A., Li, Q.,
Holzinger, E., Karyadi, D., Cannon-Albright, L.A., Teerlink, C.C., Stanford, J.L., Isaacs, W.B., Xu, J.,
Cooney, K.A., Lange, E.M., Schleutker, J., Carpten, J.D., Powell, I.J., Cussenot, O., Cancel-Tassin, G.,
Giles, G.G., MacInnis, R.J., Maier, C., Hsieh, C.L., Wiklund, F., Catalona, W.J., Foulkes, W.D., Mandal, D.,
Eeles, R.A., Kote-Jarai, Z., Bustamante, C.D., Schaid, D.J., Hastie, T., Ostrander, E.A., Bailey-Wilson, J.E.,
Radivojac, P., Thibodeau, S.N., Whittemore, A.S., Sieh, W., 2016. REVEL: an ensemble method for predicting the
pathogenicity of rare missense variants. American Journal of Human Genetics 99, 877–885.
https://doi.org/10.1016/j.ajhg.2016.08.016

Joehanes, R., Just, A.C., Marioni, R.E., Pilling, L.C., Reynolds, L.M., Mandaviya, P.R., et al., 2016.
Epigenetic signatures of cigarette smoking. Circulation: Cardiovascular Genetics 9, 436–447.
https://doi.org/10.1161/CIRCGENETICS.116.001506

Kalyakulina, A., Yusipov, I., Kondakova, E., Bacalini, M.G., Franceschi, C., Vedunova, M., Ivanchenko, M.,
2023. Small immunological clocks identified by deep learning and gradient boosting. Frontiers in Immunology
14, 1177611. https://doi.org/10.3389/fimmu.2023.1177611

Kato, S., Han, S.Y., Liu, W., Otsuka, K., Shibata, H., Kanamaru, R., Ishioka, C., 2003. Understanding the
function-structure and function-mutation relationships of p53 tumor suppressor protein by high-resolution
missense mutation analysis. Proceedings of the National Academy of Sciences of the United States of America
100, 8424–8429. https://doi.org/10.1073/pnas.1431692100

Kircher, M., Witten, D.M., Jain, P., O'Roak, B.J., Cooper, G.M., Shendure, J., 2014. A general framework for
estimating the relative pathogenicity of human genetic variants. Nature Genetics 46, 310–315.
https://doi.org/10.1038/ng.2892

Klaeger, S., Heinzlmeir, S., Wilhelm, M., Polzer, H., Vick, B., Koenig, P.A., Reinecke, M., Ruprecht, B.,
Petzoldt, S., Meng, C., Zecha, J., Reiter, K., Qiao, H., Helm, D., Koch, H., Schoof, M., Canevari, G., Casale,
E., Re Depaolini, S., Feuchtinger, A., Wu, Z., Schmidt, T., Rueckert, L., Becker, W., Huenges, J., Garz, A.K.,
Gohlke, B.O., Zolg, D.P., Kayser, G., Vooder, T., Preissner, R., Hahne, H., Tõnisson, N., Kramer, K., Götze,
K., Bassermann, F., Schlegl, J., Ehrlich, H.C., Aiche, S., Walch, A., Greif, P.A., Schneider, S., Felder,
E.R., Ruland, J., Médard, G., Jeremias, I., Spiekermann, K., Kuster, B., 2017. The target landscape of
clinical kinase drugs. Science 358, eaan4368. https://doi.org/10.1126/science.aan4368

Kolde, R., Laur, S., Adler, P., Vilo, J., 2012. Robust rank aggregation for gene list integration and
meta-analysis. Bioinformatics 28, 573–580. https://doi.org/10.1093/bioinformatics/btr709

Kotler, E., Shani, O., Goldfeld, G., Lotan-Pompan, M., Tarcic, O., Gershoni, A., Hopf, T.A., Marks, D.S.,
Oren, M., Segal, E., 2018. A systematic p53 mutation library links differential functional impact to cancer
mutation pattern and evolutionary conservation. Molecular Cell 71, 178–190.e8.
https://doi.org/10.1016/j.molcel.2018.06.012

Krantz, D.H., Luce, R.D., Suppes, P., Tversky, A., 1971. Foundations of Measurement, Vol. I: Additive and
Polynomial Representations. Academic Press, New York.

Landrum, M.J., Lee, J.M., Benson, M., Brown, G.R., Chao, C., Chitipiralla, S., Gu, B., Hart, J., Hoffman, D.,
Jang, W., Karapetyan, K., Katz, K., Liu, C., Maddipatla, Z., Malheiro, A., McDaniel, K., Ovetsky, M., Riley,
G., Zhou, G., Holmes, J.B., Kattman, B.L., Maglott, D.R., 2018. ClinVar: improving access to variant
interpretations and supporting evidence. Nucleic Acids Research 46, D1062–D1067.
https://doi.org/10.1093/nar/gkx1153

Liberzon, A., Birger, C., Thorvaldsdóttir, H., Ghandi, M., Mesirov, J.P., Tamayo, P., 2015. The Molecular
Signatures Database hallmark gene set collection. Cell Systems 1, 417–425.
https://doi.org/10.1016/j.cels.2015.12.004

Liu, X., Li, C., Mou, C., Dong, Y., Tu, Y., 2020. dbNSFP v4: a comprehensive database of transcript-specific
functional predictions and annotations for human nonsynonymous and splice-site SNVs. Genome Medicine 12, 103.
https://doi.org/10.1186/s13073-020-00803-9

Luce, R.D., Krantz, D.H., Suppes, P., Tversky, A., 1990. Foundations of Measurement, Vol. III: Representation,
Axiomatization, and Invariance. Academic Press, San Diego.

Lussier, A.A., Morin, A.M., MacIsaac, J.L., Salmon, J., Weinberg, J., Reynolds, J.N., Pavlidis, P., Chudley,
A.E., Kobor, M.S., 2018. DNA methylation as a predictor of fetal alcohol spectrum disorder. Clinical
Epigenetics 10, 5. https://doi.org/10.1186/s13148-018-0439-6

Mavrommatis, C., Belsky, D.W., Ying, K., Moqri, M., Campbell, A., Richmond, A., Gladyshev, V.N., Chandra, T.,
McCartney, D.L., Marioni, R.E., 2025. An unbiased comparison of 14 epigenetic clocks in relation to 174
incident disease outcomes. Nature Communications 16, 11164. https://doi.org/10.1038/s41467-025-66106-y
[UNVERIFIED: middle-author list (Ying, Moqri, Campbell, Richmond, Gladyshev, Chandra) drawn from open-access
listing, not confirmed against the trusted source, which abbreviated it; year/venue/volume/article-number/DOI
verified]

Metz, J.T., Johnson, E.F., Soni, N.B., Merta, P.J., Kifle, L., Hajduk, P.J., 2011. Navigating the kinome.
Nature Chemical Biology 7, 200–202. https://doi.org/10.1038/nchembio.530

Notin, P., Kollasch, A., Ritter, D., van Niekerk, L., Paul, S., Spinner, H., Rollins, N., Shaw, A., Orenbuch,
R., Weitzman, R., Frazer, J., Dias, M., Franceschi, D., Gal, Y., Marks, D.S., 2023. ProteinGym: large-scale
benchmarks for protein fitness prediction and design. In: Advances in Neural Information Processing Systems 36
(NeurIPS 2023), Datasets and Benchmarks Track, pp. 64331–64379.

Patil, G.P., Taillie, C., 2004. Multiple indicators, partially ordered sets, and linear extensions:
multi-criterion ranking and prioritization. Environmental and Ecological Statistics 11, 199–228.
https://doi.org/10.1023/B:EEST.0000027209.93218.d9

PDSP Ki Database, 2026. NIMH Psychoactive Drug Screening Program Kᵢ Database. University of North Carolina at
Chapel Hill. https://pdspdb.unc.edu/databases/kiDownload/search.php (accessed 4 May 2026).

Portales-Casamar, E., Lussier, A.A., Jones, M.J., MacIsaac, J.L., Edgar, R.D., Mah, S.M., Barhdadi, A.,
Provost, S., Lemieux-Perreault, L.P., Cynader, M.S., Chudley, A.E., Dubé, M.P., Reynolds, J.N., Pavlidis, P.,
Kobor, M.S., 2016. DNA methylation signature of human fetal alcohol spectrum disorder. Epigenetics & Chromatin
9, 25. https://doi.org/10.1186/s13072-016-0074-4

Rentzsch, P., Schubach, M., Shendure, J., Kircher, M., 2021. CADD-Splice—improving genome-wide variant effect
prediction using deep learning-derived splice scores. Genome Medicine 13, 31.
https://doi.org/10.1186/s13073-021-00835-9

Saul, D., Kosinsky, R.L., Atkinson, E.J., Doolittle, M.L., Zhang, X., LeBrasseur, N.K., Pignolo, R.J.,
Robbins, P.D., Niedernhofer, L.J., Ikeno, Y., Jurk, D., Passos, J.F., Hickson, L.J., Xue, A., Monroe, D.G.,
Tchkonia, T., Kirkland, J.L., Farr, J.N., Khosla, S., 2022. A new gene set identifies senescent cells and
predicts senescence-associated pathways across tissues. Nature Communications 13, 4827.
https://doi.org/10.1038/s41467-022-32552-1

Sayed, N., Huang, Y., Nguyen, K., Krejciova-Rajaniemi, Z., Grawe, A.P., Gao, T., Tibshirani, R., Hastie, T.,
Alpert, A., Cui, L., Kuznetsova, T., Rosenberg-Hasson, Y., Ostan, R., Monti, D., Lehallier, B., Shen-Orr,
S.S., Maecker, H.T., Dekker, C.L., Wyss-Coray, T., Franceschi, C., Jojic, V., Haddad, F., Montoya, J.G., Wu,
J.C., Davis, M.M., Furman, D., 2021. An inflammatory aging clock (iAge) based on deep learning tracks
multimorbidity, immunosenescence, frailty and cardiovascular aging. Nature Aging 1, 598–615.
https://doi.org/10.1038/s43587-021-00082-y

Scicluna, B.P., van Vught, L.A., Zwinderman, A.H., Wiewel, M.A., Davenport, E.E., Burnham, K.L., Nürnberg, P.,
Schultz, M.J., Horn, J., Cremer, O.L., Bonten, M.J., Hinds, C.J., Wong, H.R., Knight, J.C., van der Poll, T.,
MARS Consortium, 2017. Classification of patients with sepsis according to blood genomic endotype: a
prospective cohort study. Lancet Respiratory Medicine 5, 816–826.
https://doi.org/10.1016/S2213-2600(17)30294-1

Shavlakadze, T., Xiong, K., Donne, R., Glass, D.J., et al., 2026. Multitissue, multi-time point transcriptomic
atlas of aging in mice and rats. Science Advances. https://doi.org/10.1126/sciadv.ady8401  [UNVERIFIED:
complete author list and volume/eLocator to confirm at proof]

Shenker, N.S., Polidoro, S., van Veldhoven, K., Sacerdote, C., Ricceri, F., Birrell, M.A., Belvisi, M.G.,
Brown, R., Vineis, P., Flanagan, J.M., 2013. Epigenome-wide association study in the European Prospective
Investigation into Cancer and Nutrition (EPIC-Turin) identifies novel genetic loci associated with smoking.
Human Molecular Genetics 22, 843–851. https://doi.org/10.1093/hmg/dds488

Singleton, J., Booth, R., 2022. Towards an axiomatic approach to truth discovery. Autonomous Agents and
Multi-Agent Systems 36, 42. https://doi.org/10.1007/s10458-022-09569-3

Sundaram, L., Gao, H., Padigepati, S.R., McRae, J.F., Li, Y., Kosmicki, J.A., Fritzilas, N., Hakenberg, J.,
Dutta, A., Shon, J., Xu, J., Batzoglou, S., Li, X., Farh, K.K.H., 2018. Predicting the clinical impact of
human mutation with deep neural networks. Nature Genetics 50, 1161–1170.
https://doi.org/10.1038/s41588-018-0167-z

Suppes, P., Krantz, D.H., Luce, R.D., Tversky, A., 1989. Foundations of Measurement, Vol. II: Geometrical,
Threshold, and Probabilistic Representations. Academic Press, San Diego.

Szpilrajn, E., 1930. Sur l'extension de l'ordre partiel. Fundamenta Mathematicae 16, 386–389.
https://doi.org/10.4064/fm-16-1-386-389

Tsaprouni, L.G., Yang, T.P., Bell, J., Dick, K.J., Kanoni, S., Nisbet, J., Viñuela, A., Grundberg, E., Nelson,
C.P., Meduri, E., Buil, A., Cambien, F., Hengstenberg, C., Erdmann, J., Schunkert, H., Goodall, A.H.,
Ouwehand, W.H., Dermitzakis, E., Spector, T.D., Samani, N.J., Deloukas, P., 2014. Cigarette smoking reduces
DNA methylation levels at multiple genomic loci but the effect is partially reversible upon cessation.
Epigenetics 9, 1382–1396. https://doi.org/10.4161/15592294.2014.969637

Uitdehaag, J.C.M., Zaman, G.J.R., 2011. A theoretical entropy score as a single value to express inhibitor
selectivity. BMC Bioinformatics 12, 94. https://doi.org/10.1186/1471-2105-12-94

van der Laan, L., Rooney, K., Trajkova, S., Haghshenas, S., Relator, R., Lauffer, P., Vos, N., Levy, M.A.,
Brunetti-Pierri, N., Bedeschi, M.F., et al., 2025. Discovery of a DNA methylation episignature as a molecular
biomarker for fetal alcohol syndrome. Genetics in Medicine 27, 101586.
https://doi.org/10.1016/j.gim.2025.101586

Yusipov, I., Kondakova, E., Kalyakulina, A., Krivonosov, M., Lobanova, N., Bacalini, M.G., Franceschi, C.,
Vedunova, M., Ivanchenko, M., 2022. Accelerated epigenetic aging and inflammatory/immunological profile
(ipAGE) in patients with chronic kidney disease. GeroScience 44, 817–834.
https://doi.org/10.1007/s11357-022-00540-4

Zdrazil, B., Felix, E., Hunter, F., Manners, E.J., Blackshaw, J., Corbett, S., de Veij, M., Ioannidis, H.,
Mendez Lopez, D., Mosquera, J.F., Magarinos, M.P., Bosc, N., Arcila, R., Kizilören, T., Gaulton, A., Bento,
A.P., Adasme, M.F., Monecke, P., Landrum, G.A., Leach, A.R., 2024. The ChEMBL Database in 2023: a drug
discovery platform spanning multiple bioactivity data types and time periods. Nucleic Acids Research 52,
D1180–D1192. https://doi.org/10.1093/nar/gkad1004
