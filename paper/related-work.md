# Related work — verified references and positioning

All citations below were verified (authors/year/venue/claim). **Accuracy fixes** applied vs earlier
drafts are flagged ⚑. "Could not fully verify" items are listed at the end — confirm before submission.
The honest bottom line (Part 3) is unchanged: the formalism is prior art; the contribution is the
protocol + the constitutive/artifactual reproducibility test + the cross-domain biological demonstration.

## Foundations the protocol builds on (cite as such — not as novelty)

- **Krantz, Luce, Suppes & Tversky, *Foundations of Measurement*** — Vol. I *Additive and Polynomial
  Representations* (Krantz, Luce, Suppes, Tversky, 1971); Vol. II *Geometrical, Threshold, and
  Probabilistic Representations* (Suppes, Krantz, Luce, Tversky, **1989**); Vol. III *Representation,
  Axiomatization, and Invariance* (Luce, Krantz, Suppes, Tversky, 1990). Academic Press. ⚑ Author order
  rotates by volume — cite each volume with its own order. Measurement = homomorphism from an empirical
  relational structure to a numerical one; scale types via admissible transformations.
- **Szpilrajn, E. (1930).** "Sur l'extension de l'ordre partiel." *Fundamenta Mathematicae* 16:386–389.
  Every partial order extends to a linear order — the "each measure is a linear extension" backbone.
- **Patil, G.P. & Taillie, C. (2004).** "Multiple indicators, partially ordered sets, and linear
  extensions: multi-criterion ranking and prioritization." *Environmental and Ecological Statistics*
  11(2):199–228. ⚑ (subtitle was missing before.) The foundational precedent: disagreeing indicators →
  partial order; composites force arbitrary weights; linear extensions rank without weights.
- **Brüggemann, R. & Patil, G.P. (2011).** *Ranking and Prioritization for Multi-indicator Systems:
  Introduction to Partial Order Applications.* Springer. ⚑ Book-length Hasse-diagram methodology;
  weight-free weak orders.
- **Average rank over linear extensions (the canonical weight-free summary):** ⚑ correct authorship —
  **De Loof, K., De Baets, B., De Meyer, H. & Brüggemann, R. (2008),** "A Hitchhiker's Guide to Poset
  Ranking," *Comb. Chem. High Throughput Screen.* 11(9):734–744; **De Loof, De Baets & De Meyer (2011),**
  "Approximation of Average Ranks in Posets," *MATCH Commun. Math. Comput. Chem.* 66(1):219–229;
  **Lerche, D., Sørensen, P.B. & Brüggemann, R. (2003),** "Improved Estimation of the Ranking
  Probabilities in Partial Orders …," *J. Chem. Inf. Comput. Sci.* 43(5):1471–1480, doi:10.1021/ci0300036.
  ⚑ (Earlier cite "De Loof, Brüggemann, De Baets" was wrong — De Meyer is the consistent co-author; and
  Lerche 2003 has a third author, **Sørensen**, previously omitted.)
- **Bubley, R. & Dyer, M. (1999).** "Faster random generation of linear extensions." *Discrete
  Mathematics* 201:81–88. The near-uniform linear-extension sampler **we use** to estimate average ranks
  (`analysis/avg_rank.py`). **Must cite** (currently missing).
- **Brightwell, G. & Winkler, P. (1991).** "Counting linear extensions." *Order* 8:225–242.
  #P-completeness — why average ranks must be sampled/approximated, not enumerated.

## Framing / contrast

- **Campbell, D.T. & Fiske, D.W. (1959).** "Convergent and discriminant validation by the
  multitrait-multimethod matrix." *Psychological Bulletin* 56(2):81–105. ⚑ exact title. Convergent
  validity is the frame; psychometrics treats mid-range divergence as *method variance to remove* — our
  departure is to call part of it *constitutive*.
- **Arrow, K.J. (1951; 2nd ed. 1963).** *Social Choice and Individual Values.* Impossibility of a
  canonical aggregator. **Singleton, J. & Booth, R. (2022),** "Towards an axiomatic approach to truth
  discovery," *Auton. Agents Multi-Agent Syst.* 36(2):42 — impossibility for truth discovery (already
  cited by the kinase repo).
- **Bland, J.M. & Altman, D.G. (1986),** *Lancet* 327:307–310 (limits of agreement); ⚑ **proportional
  bias** specifically is **Bland & Altman (1995),** *Lancet* 346:1085–1087 — cite the 1995 paper if the
  "agree-at-extremes/diverge-mid-range" pattern is invoked.

## Measurement theory in biology (bucket A)
- **Houle, D., Pélabon, C., Wagner, G.P. & Hansen, T.F. (2011).** "Measurement and Meaning in Biology."
  *Quart. Rev. Biol.* 86(1):3–34. The key precedent for RMT-in-biology; legitimizes our move but treats
  single measurements/scale types, not competing measures of one concept as extensions of a shared poset.

## Posets in bioinformatics / cheminformatics / indicators (bucket B)
- **Halfon, E. & Reggiani, M.G. (1986),** *Environ. Sci. Technol.* 20(11):1173–1179 — Hasse-diagram
  ranking of chemicals without weights (direct ancestor).
- **Fattore, M. (2016),** "Partially Ordered Sets and the Measurement of Multidimensional Ordinal
  Deprivation," *Soc. Indic. Res.* 128(2):835–858 (PARSEC) — strongest social-science analogue.
- Brüggemann/Voigt PyHasse & poset-in-chemistry literature (e.g. PMID 21534922).

## Canonical/weight-free ranking machinery (bucket C)
- De Loof et al. (2008/2011), Lerche et al. (2003), Bubley & Dyer (1999), Brightwell & Winkler (1991) —
  above.

## "Which measure to trust" in biology (bucket D)
- **Fan, C. et al. (2006),** "Concordance among Gene-Expression–Based Predictors for Breast Cancer,"
  *NEJM* 355(6):560–569 — different signatures, few shared genes, concordant patient predictions:
  empirical instance of "agree where the order is determined." Ideal motivating example.
- **Uitdehaag, J.C.M. & Zaman, G.J.R. (2011),** *BMC Bioinformatics* 12:94 (selectivity entropy) and
  **Bosc, N. et al. (2017),** *BMC Bioinformatics* 18:17 (comparison of kinase selectivity metrics) —
  the prior art for the kinase running example; documents metric disagreement directly.
- **Belsky, D.W. et al. (2022),** *eLife* 11:e73420 (DunedinPACE) + multi-clock comparisons — the
  biological-age example and its reliability contrasts (motivates the reproducibility test).

## No-canonical-aggregator / rank aggregation (bucket E)
- Arrow (1951); Singleton & Booth (2022) — above.
- **Dwork, Kumar, Naor & Sivakumar (2001),** WWW'01 — Kemeny consensus NP-hard; best aggregator
  intractable/non-unique.
- **Kolde, R. et al. (2012),** "Robust rank aggregation," *Bioinformatics* 28(4):573–580 — the dominant
  bioinformatics rank-aggregation method; contrast paradigm (force consensus vs preserve incomparability).
- **Marbach, D. et al. (2012),** *Nat. Methods* 9(8):796–804 (wisdom-of-crowds networks); **Li, Y. et
  al. (2016),** truth-discovery survey — the "combine everything" ethos our protocol qualifies.

## What is genuinely left to contribute (honest)
Not the formalism (all prior art). The novelty is: **(1)** a reusable *protocol* packaging these pieces
for the recurring biology situation of several competing measures of one latent concept; **(2)** the
**constitutive-vs-artifactual test via cross-cohort reproducibility** (no prior work partitions poset
incomparabilities this way — apparently distinctive); **(3)** the **cross-domain biological
demonstration** (one protocol + one universal near-tie law across selectivity, age, exposure — see
`analysis/transfer.py`).

## Additional verified references (resolved this session)
- **Sawamura, J., Morishita, S. & Ishigooka, J. (2014).** "Interpretation for scales of measurement
  linking with abstract algebra." *J. Clinical Bioinformatics* 4:9, doi:10.1186/2043-9113-4-9. (Bucket A —
  measurement-scale rigor in bioinformatics.)
- **Nardo, M., Saisana, M., Saltelli, A., Tarantola, S., Hoffman, A. & Giovannini, E.** *Handbook on
  Constructing Composite Indicators: Methodology and User Guide* — OECD/JRC book **(2008)**, ISBN
  978-92-64-04345-9, doi:10.1787/9789264043466-en; earlier **OECD Statistics Working Paper 2005/3**,
  doi:10.1787/533411815016. (Bucket B — the "weight-and-aggregate" practice the poset approach reacts
  against.) ⚑ correct initials: Saltelli **A.**, Tarantola **S.**
- **Michell, J. (1997),** "Quantitative science and the definition of measurement in psychology,"
  *Br. J. Psychol.* 88(3):355–383; **Michell, J. (1999),** *Measurement in Psychology: A Critical History
  of a Methodological Concept*, Cambridge UP. (Bucket A — the "quantitative imperative": don't assume
  quantitative structure without testing for it; supports our ordinal-only caution.)
- **Mavrommatis, C., Belsky, D.W., … McCartney, D.L. & Marioni, R.E. (2025).** "An unbiased comparison
  of 14 epigenetic clocks in relation to 174 incident disease outcomes." *Nature Communications*
  16:11164, doi:10.1038/s41467-025-66106-y (Generation Scotland, n≈18,859). (Bucket D — clocks disagree /
  differ in predictive value.)

## Still genuinely uncertain (confirm against your copy before submission)
- *Foundations of Measurement* Vol. II **year 1989 vs 1990** — reliable catalogs disagree (Academic Press
  1989 vs some listings/Psychometrika review 1990). ⚑ Vol. II author order is **Suppes, Krantz, Luce,
  Tversky** (not the Vol. I order) — corrected above.
- Michell (1997) page range **355–383** — matches every secondary source; DOI/volume/issue firm, but not
  read off the paywalled primary page.
(All other previously-uncertain items are now resolved above / inline.)
