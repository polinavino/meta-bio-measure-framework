# meta-bio-measure-framework

A methods/Perspective paper: **a protocol for comparing competing measures of a latent biological
concept.** When several formally distinct measures each claim to quantify one concept (kinase
*selectivity*, *biological age*, *exposure*, a *disease signature*) and disagree, the field usually
benchmarks them to pick a winner. This project separates the choice of definition from prediction benchmarks and gives a comparison
protocol, grounded in established measurement theory and demonstrated across seven domains.

## In plain language

Competing biological formulas can rank the same objects differently. This project gives a protocol for
finding out how much of that disagreement survives changes in implementation, eligibility, panels,
missing-data conventions and cohorts. Stable disagreement is an observable result; attributing it to the
biological concept itself additionally requires assumptions about which formulas are admissible and
controls for systematic artifacts.

Where all admitted measures weakly agree, the consensus retains that comparison. A strict reversal makes
a pair incomparable; a tie in one measure and a strict comparison in another does not. If a scalar summary
is needed, we average class ranks over uniform linear extensions of the consensus quotient poset. This
is an explicit aggregation convention that can resolve incomparable pairs. It has no measure-specific
weights, but does not eliminate choices about the base set or the distribution over rankings.

## Current analysis status (November 2026)

The post-submission kinase revisions prompted a framework-wide correction. Current analyses preserve ties,
report strict reversals separately, freeze data-derived sepsis preprocessing and directions on validation, and use four
initially dispersed insertion-Gibbs chains for average ranks. Split R-hat and batch-means Monte Carlo
precision must pass before aggregate benchmarks are reported. Identical indicator vectors are collapsed
to classes before sampling; tied objects receive the same class score.

- Kinase: 37.9% incomparable under the declared parameter-median implementations; fixed-default formula
  and floor scenarios are separate sensitivity analyses. Active count is a shared-input proxy, with a
  S-score-withheld control. Matched-panel incomparable-set Jaccard is 0.294 across 40 shared eligible compounds.
- TP53: median-of-eight agreement with the checked aggregate is 0.955, replacing the insufficiently mixed
  estimate of 0.987. Population and endpoint shifts also limit interpretation of the ClinVar/VUS comparison.
- Inflammatory indices: the checked aggregate agrees with SII at 0.844; mean rank agrees at 0.983. CRP/CAR
  correlate at 0.998 but are not exactly rank-identical.
- Sepsis: these are gene-set proxies, not official classifier reconstructions. Frozen discovery signs give
  matrix-structure reproducibility 0.976, with participant-bootstrap interval [0.684, 0.988].
- Serotonin: observed receptor panels vary; missing receptors remain unknown in the primary analysis.
  Entropy strongly tracks panel size, and a fixed-panel control changes the apparent clustering.
- CVD: common numerical cutoffs are classification sensitivity checks. The PCE and Framingham endpoints
  differ, so discrepancies are not interpreted as treatment changes or calibration error.

`paper/manuscript.md` is authoritative. Outputs in `analysis/outputs/reaudit.txt` and the domain output
directories document validation; `paper/check_claims.py` checks registered quantitative claims against
stored results. Older working notes are historical wherever they conflict.

## Repo map
- **`paper/manuscript.md`** — **the canonical submission draft** (Journal of Theoretical Biology).
  Every quantitative claim in it is copied from a stored analysis output. Section numbers cited
  elsewhere in this README refer to this file.
- **`paper/main.md`**, **`paper/formal-spine.md`** — working notes the manuscript draws on. `formal-spine.md`
  is the source for Appendix A (classical results cited, not claimed). Where either disagrees with
  `manuscript.md`, the manuscript wins.
- **`paper/references-formatted.md`** — reference list in Elsevier name-date style, with a DOI-verification
  log at the end. A superset of the manuscript's list.
- **`paper/related-work.md`** — prior-art notes.
- **`synthesis/repo-findings.md`** — verified ground truth from the four source repos; **§6** has the
  controlled re-analyses, prior-art map, and decision log. **Start here** for provenance.
- **`analysis/`** — the controlled re-analyses behind §4.1 to §4.4 and §5. `run_all.sh` regenerates the current core scientific
  results in `analysis/outputs/` (the initial audit report is historical), which is the source of truth for those numbers. The scripts read the four
  source repos by absolute path, so those must be cloned as siblings.
- **`tp53/`** — TP53 variant effect (§4.5), computed end-to-end on public data (ProteinGym DMS +
  AlphaMissense + NCI/Kato transactivation + ClinVar). `tp53/README.md` writeup, `tp53/analysis/run_all.sh`
  regenerates `tp53/analysis/outputs/`. See `synthesis/repo-findings.md` §6.7.
- **`inflammation/`** — inflammatory burden (§4.6), three sub-topics on public data: clinical
  indices/NHANES (strong), sepsis signatures/GSE65682 (strong disagreement, weak reproducibility),
  inflammaging clocks (thin). `inflammation/README.md` documents where the candidate measures come from
  and the base-set/completeness caveat. `inflammation/analysis/run_all.sh` regenerates its outputs.
- **`paper/framework-core.md`**, **`meta-paper-plan.md`** — SUPERSEDED / historical (old "new-math" framing).

### The seven domains and where each one lives
Four of the original instances were computed in separate sibling repositories, which have to be archived
alongside this one for the paper's data-availability statement to hold.

| § | Domain | Instance repository |
|---|---|---|
| 4.1 | Kinase selectivity | `../selectivity` (preprint on ChemRxiv) |
| 4.2 | Epigenetic ageing, incl. the rat retina | `../epigenetic-clock-desiderata`, rat retina in its `eye_aging/` |
| 4.3 | Serotonin-receptor selectivity | `../psychedelic-selectivity` |
| 4.4 | Methylation exposure biomarkers | `../methylation-biomarker-agreement` |
| 4.5 | TP53 variant effect | `tp53/` (this repo) |
| 4.6 | Inflammatory burden | `inflammation/` (this repo) |
| 4.7 | Cardiovascular risk | `cvdrisk/` (this repo) |

## Historical status (2026-09; see November corrections above)
Re-scoped from a new-mathematics claim to a methods/Perspective paper after three adversarial referees
and controlled re-analyses. Core empirical claims are honest and controlled (near-ties, not a
"consequential middle"; detection-floor instability; cross-cohort reproducibility). Formalism is
borrowed prior art (Patil–Taillie; RMT; Szpilrajn; average-rank canonical extension; Campbell–Fiske).
**Target venue decided: Journal of Theoretical Biology**, with `paper/manuscript.md` written to it and
`paper/references-formatted.md` in its reference style. The venue survey below is retained as the
reasoning behind that choice and as a fallback list; it is not a live recommendation.

## Analysis controls imported from the kinase revisions

The protocol now declares the biological target unit, exact formula, affinity weighting, missingness,
undefined scores, eligibility and parameter contexts before interpreting disagreements. Published and
modified entropy forms are compared explicitly. Boundary locality, robustness to boundary crossings,
orientation and panel recovery are distinct properties; construction guarantees state their domain and
scale conditions. The current kinase candidate is gated exponential entropy, with T >= tau* - x0 for the
weak-addition guarantee. Historical softplus/hinge experiments retain their old criterion numbering only
as records of exploration.

Recovery reports eligibility retention, conditional rank agreement and pair recovery that counts lost
eligibility as zero credit. Base-set omissions, participant/compound bootstraps, fixed-panel serotonin
controls and matched-compound kinase comparisons expose sensitivity. Replication establishes stability
conditional on these choices, not the cause of disagreement. Uniform extension averaging is a stated
convention, not a uniquely valid biological measure.

Run the full local workflow with `bash analysis/run_framework.sh`. Python with numpy/pandas/scipy/sklearn,
R with dplyr, and `clang++` for the small native sampler are required; source data and sibling repositories
must already be present. No source-domain repository is modified by this workflow. Input/code SHA256 hashes
and sibling HEAD revisions are recorded by `analysis/provenance.py`.

---

## Venue survey (historical: the reasoning behind choosing JTB, kept as a fallback list)

Preprint to **bioRxiv** (and cross-list arXiv q-bio / math.ST for the formal-methods audience) **first,
regardless of target** — it establishes priority and is expected by every venue below. Prerequisite
note: the kinase preprint is live (ChemRxiv); the clock manuscript does not yet exist (repo/README only)
— the paper currently cites the repo, so no venue is blocked, but a posted clock preprint would
strengthen submission.

| Venue | Type to submit | Fit | Prestige | Accept odds | Main risk |
|---|---|---|---|---|---|
| **PLOS Computational Biology** | Perspective (or Methods) | High | Good | Medium–High | Must read as broadly useful to comp-bio, not niche |
| **Patterns** (Cell Press) | Perspective / Methodology | High if generality is demonstrated | High | Medium | "Does the frame earn its cross-domain generality?" |
| **Briefings in Bioinformatics** | Review/Perspective | High | Medium | High | Reads as a review unless the protocol is foregrounded |
| **GigaScience / GigaByte** | Technical/Methods + code | Medium | Medium | Medium | Wants a packaged tool + reproducible pipeline |
| **Bioinformatics** (OUP) | Application Note / Methods | Medium | Good | Medium | Expects usable software, not a protocol on paper |
| **J. Mathematical Psychology** / measurement-theory venue | Research article | Formal core only | Niche | Low–Medium | Biology examples are out of scope; audience mismatch |

### Per-venue notes

**PLOS Computational Biology — Perspective (recommended, see below).**
- *Pros:* reputable, open access, broad comp-bio readership; explicitly welcomes conceptual/opinion
  pieces that reframe how a field thinks; tolerant of formal content in an appendix; strong for building
  a distinctive public identity ("applied measurement theory for computational biology"). The
  honest-empirical-core + protocol shape fits a Perspective well.
- *Cons:* editors gate on "is this broadly interesting to computational biologists"; the formal-methods
  framing must be pitched as *useful*, not as self-contained theory; empirical core is modest.

**Patterns (Cell Press) — Perspective/Methodology (highest upside).**
- *Pros:* its remit is exactly cross-disciplinary data-science methodology and general principles across
  domains; higher profile; a "here is a domain-independent protocol" paper is on-brand. Best venue for
  impact + identity *if* the generality is convincingly demonstrated.
- *Cons:* the referee's sharpest objection lands here — does the unification *earn* its generality, or is
  it four notes plus a shared metaphor? Would likely require the cross-domain transfer demonstration we
  have not done (derive where measures disagree in one domain from another's structure). More work,
  higher rejection risk.

**Briefings in Bioinformatics — Perspective/Review (reliable fallback).**
- *Pros:* practitioner audience that would actually use the protocol; welcomes methodological
  perspectives and surveys; high acceptance odds; fast-ish. Good if speed and a receptive audience
  matter more than prestige.
- *Cons:* lower prestige; gravitational pull toward "review" — must keep the protocol (not the survey)
  as the centerpiece or it reads as yet-another-review.

**GigaScience / Bioinformatics (tool-first venues).**
- *Pros:* if we package the protocol as a small reusable package (compute the consensus poset,
  reproducibility split, average-rank extension, desiderata table) with the `analysis/` scripts as a
  worked example, either becomes a strong fit and adds lasting utility.
- *Cons:* both effectively require the software artifact; without it the paper is a poor match. This is a
  strategic choice (build the tool) more than a venue choice.

**Measurement-theory / math-psych venues.**
- *Pros:* the formal core would be appreciated; signals the formal-methods identity strongly.
- *Cons:* the biological instances — the whole point — are out of scope; wrong audience for impact in
  comp-bio. Only sensible as a *second*, formal companion, not the flagship.

### Decision
**Journal of Theoretical Biology**, and `paper/manuscript.md` is written to it. JTB takes the formal
core in the body rather than only in an appendix, which the venues above would not, and it is free
hybrid, so no APC. The survey above ranked PLOS Computational Biology first on comp-bio reach; that
remains the strongest alternative if JTB declines.

Fallbacks in order: **PLOS Computational Biology** (Perspective), **Briefings in Bioinformatics**
(fast, receptive), **Patterns** (only if a genuine cross-domain transfer result is added), and
**GigaScience/Bioinformatics** (only if the protocol is packaged as software). Reserve a math-psych
venue for a possible formal companion paper.
