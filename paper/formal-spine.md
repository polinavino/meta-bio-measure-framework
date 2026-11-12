# Formal justification — synchronized working copy

The authoritative version is Appendix A of `manuscript.md`. This copy follows the November 2026
corrections: quotient preorders before sampling, distinguish weak and linear extensions, and separate
replication from admissibility and causal attribution.

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
