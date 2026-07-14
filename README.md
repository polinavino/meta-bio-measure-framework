# meta-bio-measure-framework

A methods/Perspective paper: **a protocol for comparing competing measures of a latent biological
concept.** When several formally distinct measures each claim to quantify one concept (kinase
*selectivity*, *biological age*, *exposure*, a *disease signature*) and disagree, the field usually
benchmarks them to pick a winner. This paper argues that is mis-specified and gives a better procedure,
justified with established measurement theory and demonstrated across four domains.

## In plain language

Lots of things in biology are measured by more than one formula at the same time. How "selective" a
drug is has at least four competing scores; someone's "biological age" has a dozen different "clocks";
smoking exposure or a disease has several rival molecular signatures. These aren't different versions of
one agreed formula — they are genuinely different formulas, and run on the same data they often **rank
the same things differently.**

The usual reaction is to treat this as a contest: benchmark the formulas and declare one the winner. We
think that's the wrong move, because the disagreement usually isn't a sign that one formula is broken.
It's a sign that **the concept itself doesn't fully decide the answer.** At the extremes everyone agrees
— a drug that hits one target and nothing else *is* selective; a newborn *is* biologically younger than
a 90-year-old. But in the messy middle — a drug with two strong targets and five weak ones, a
middle-aged person — the concept is genuinely silent about who ranks higher, and each formula quietly
makes its own arbitrary call. So "which formula is correct?" is the wrong question. The useful questions
are: *where does the concept actually decide, where is it silent, and is that pattern real (does it show
up again in a fresh dataset) or just noise?*

**Combining multiple measures (the part worth explaining).** Once you accept that no single formula is
"the right one," the natural thing to do is not to pick one but to **combine all the reasonable ones into
a consensus.** Think of the competing formulas as a panel of judges ranking contestants:

- Wherever **all** the judges agree that A beats B, that ordering is solid — keep it.
- Wherever the judges **split** on A vs B, don't crown one judge as correct. Instead, consider *every*
  possible full ranking that is consistent with everything the judges *did* agree on, and take each
  contestant's **average position** across all of those rankings.

That average is a single consensus ranking with two nice properties: it agrees with the formulas exactly
where they agree, and on the cases they disagree about it reports the *balance of opinion* (how the
formulas split) instead of forcing a coin-flip. It needs no weighting (you never have to declare one
formula more important) and has no tuning knobs. When you genuinely need one number, this consensus is
the honest one to report — it commits to nothing beyond what the measures actually agree on. (The formal
name is the *average rank over the linear extensions of the consensus order*; the panel-of-judges picture
is exactly what that computes.)

A single formula is still fine and often more convenient — you can compute it on one item in isolation,
and it can have an intuitive physical meaning. But if you're able to look at the whole space of
reasonable formulas at once, the consensus is the less arbitrary summary.

## Repo map
- **`paper/main.md`** — the manuscript draft (methods/Perspective).
- **`paper/formal-spine.md`** — Appendix A: formal justification of the protocol (corrected; classical
  results cited, not claimed).
- **`synthesis/repo-findings.md`** — verified ground truth from the five source repos; **§6** has the
  controlled re-analyses, prior-art map, and decision log. **Start here.**
- **`analysis/`** — scripts reproducing the controlled findings (§6.2) + their README.
- **`paper/framework-core.md`**, **`meta-paper-plan.md`** — SUPERSEDED / historical (old "new-math" framing).

## Status (2026-07)
Re-scoped from a new-mathematics claim to a methods/Perspective paper after three adversarial referees
+ controlled re-analyses. Core empirical claims are honest and controlled (near-ties, not a
"consequential middle"; detection-floor instability; cross-cohort reproducibility). Formalism is
borrowed prior art (Patil–Taillie; RMT; Szpilrajn; average-rank canonical extension; Campbell–Fiske).
Not yet tailored to a target venue (see below).

## Refinements from the kinase cross-dataset validation — FOLDED IN ✓

The kinase candidate was validated across all four datasets (orientation correct everywhere, D4 100%,
T-robust, gate-necessity shown by the un-gated inversion). That work surfaced protocol/framework
refinements. **All four are written into the paper:** #1 → `paper/main.md` §2 step 1 + the
"Sharpening the desiderata" note (G2); #2, #3 → same note (G3, G1) + `synthesis/repo-findings.md` §2.5
table; #4 → `paper/formal-spine.md` §1.1. **#1 and #4 are now DEMONSTRATED across all four domains**
(not just kinase-asserted) via `analysis/external_anchor.py` — every standard measure correctly oriented
against its external anchor, the un-gated candidate control caught, and consensus–anchor agreement +0.65
to +0.94 (stronger at the extremes); see main.md §5 point 6. The notes below are the rationale/record.

1. **Orientation check against an external anchor → sharpen G2 (monotonicity in the concept).**
   The kinase candidate passed stability, monotonicity, and panel-convergence while being *oriented
   backwards* (it measured promiscuity, not selectivity); the internal desiderata missed it, a crude
   external proxy (n_active) caught it. Lesson: **verify G2 against a theory-light external proxy of the
   concept's direction, not just internal consistency** — n_active (selectivity), chronological age
   (clocks), self-reported pack-years (smoking). This also breaks a circularity: `≽` is currently
   estimated purely as the consensus of the measures, so a collectively mis-oriented set would go
   uncaught; an external anchor is the independent check. (Motivating example: the un-gated inversion,
   −0.943 vs the distribution family.)

2. **Distinguish apparatus-fixed vs analyst-chosen nuisance parameters → sharpen G3 (ordinal stability
   under nuisance).** The D3 "is the floor a free baseline or a fixed physical constant?" confusion is
   general. **G3 should require robustness only to analyst-chosen parameters; apparatus-fixed ones
   (detection floor, assay type) are declared context, not things to vary.** Conflating them cuts both
   ways — dinging a measure for a parameter that isn't free, or "passing" it by testing the wrong knob.
   Empirical wrinkle worth stating: boundary-sensitivity is itself data-dependent (Anastassiadis +0.91
   vs Klaeger −0.93), because it depends on how many objects sit near the boundary.

3. **Reliability gate (G1) can be load-bearing for correctness, not just a "flag unreliable" nicety.**
   Without the gate the kinase measure doesn't merely get noisy at low signal — it *inverts*. So G1 is
   more central than "report NA below threshold": gating the low-signal region can be what keeps a
   measure correctly oriented at all.

4. **(Bonus) External anchor also validates the consensus order.** The same external proxy (idea 1) can
   check that the consensus `≽` agrees with a crude independent proxy at the extremes (where both should
   be unambiguous) — an external-validity check on the skeleton, complementing the internal cross-cohort
   reproducibility test.

Also: the cross-dataset validation turns "a valid closed-form candidate exists" (option i) into
"validated across four assay technologies," which strengthens the (i)-vs-(ii) framing in §5.

Placement when written up: #1 → §2 protocol (an orientation step) + the §2 "Sharpening the desiderata" note (G2); #2, #3 → that same §2 note (G3, G1);
#4 → the `≽`-estimation / reproducibility discussion (formal-spine §1.1). Running illustration: the
kinase inversion + cross-dataset validation (`selectivity/candidate_validation.py`).

---

## Candidate venues

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

### Recommendation
**Target PLOS Computational Biology as a Perspective**, with bioRxiv+arXiv preprint first. Rationale:
it is the best balance of (a) reach into the intended comp-bio audience, (b) tolerance for the
formal-methods angle that is the author's differentiator, (c) realistic acceptance odds for an honest
methodological Perspective with a modest empirical core, and (d) identity-building value. It does not
demand the cross-domain-transfer demonstration that **Patterns** would likely require.

Keep **Patterns** as the stretch option: if we later add a genuine cross-domain transfer result (and/or
package the tool), reposition there for higher impact. Keep **Briefings in Bioinformatics** as the
fast, receptive fallback. Treat **GigaScience/Bioinformatics** as contingent on deciding to build a
software package. Reserve a math-psych venue only for a possible formal companion paper.

*(Decision not yet made; no paper-tailoring done. Confirm the target before we adapt length, format,
and framing to it.)*
