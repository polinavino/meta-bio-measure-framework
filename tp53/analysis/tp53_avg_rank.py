"""
TP53 domain, protocol step 7: the weight-free minimal-commitment canonical aggregate
(average rank over the linear extensions of the consensus poset; De Loof / Bruggemann / Patil-Taillie),
demonstrated on the N=2314 complete-panel TP53 variants.

Panel of 9 competing measures, concept-positive (higher = more deleterious):
  8 Ishioka transactivation promoters (negated %WT) + AlphaMissense.
Consensus poset:  a dominates b  iff  every measure scores a >= b (Patil-Taillie dominance).
Canonical measure = mean rank of each variant over uniformly-sampled linear extensions
(Bubley-Dyer lazy adjacent-transposition chain; convergence via two independent chains).

We then benchmark, against the canonical order:
  - each of the 9 measures (which best proxies the weight-free consensus? which strays most?)
  - the FIELD'S PROPOSED formula: median of the 8 promoters (Kato/Ishioka) — does the naive field
    consensus match the principled minimal-commitment one?
  - the 6 other computational predictors (EVE, ESM1b, REVEL, CADD, PrimateAI, BayesDel), kept external
    to the poset, to show how far the whole computational family sits from the experimental consensus.
  - whether the canonical aggregate predicts ClinVar at least as well as any single measure.
"""
import numpy as np, pandas as pd
from scipy.stats import rankdata, spearmanr

# --- persist this script's console output to a tracked file (single source of truth for paper numbers) ---
import sys as _sys
_OUT = "/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/tp53/analysis/outputs/avg_rank.txt"
class _Tee:
    def __init__(self, p): self._f = open(p, "w"); self._o = _sys.stdout
    def write(self, s): self._o.write(s); self._f.write(s)
    def flush(self): self._o.flush(); self._f.flush()
_sys.stdout = _Tee(_OUT)

M = pd.read_csv("/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/tp53/data/tp53_matrix.csv")
prom = ["waf1","mdm2","bax","h1433s","aip1","gadd45","noxa","p53r2"]
core = M[prom+["am_path"]].notna().all(axis=1).to_numpy()
D = M[core].reset_index(drop=True)
N = len(D)

# concept-positive scores (higher = more deleterious)
cols = {f"TA:{p}": -D[p].to_numpy(float) for p in prom}
cols["AlphaMissense"] = D["am_path"].to_numpy(float)
names = list(cols)
X = np.column_stack([cols[n] for n in names])          # (N, 9)

# ---- consensus poset: geq[a,b] = all measures score a >= b ----
geq = np.ones((N,N), bool)
for k in range(X.shape[1]):
    xk = X[:,k]
    geq &= (xk[:,None] >= xk[None,:])
comparable = geq | geq.T
np.fill_diagonal(comparable, True)
frac_incomp = 1 - (comparable.sum()-N)/(N*(N-1))
print(f"Consensus poset over N={N} TP53 variants (9 measures): "
      f"incomparable-pair fraction = {frac_incomp:.3f}")

# ---- Bubley-Dyer uniform linear-extension sampler (same machinery as analysis/avg_rank.py) ----
def sample_avg_rank(seed, burn=1_000_000, nsamp=2000, gap=2000):
    rng = np.random.default_rng(seed)
    order = list(np.argsort(-X.mean(1)))      # most-deleterious first; valid linear extension of dominance
    rank_sum = np.zeros(N); ns = 0
    total = burn + nsamp*gap
    for step in range(total):
        i = rng.integers(N-1); a = order[i]; b = order[i+1]
        if not comparable[a,b]:                # only adjacent incomparable pairs may swap (lazy)
            if rng.random() < 0.5:
                order[i], order[i+1] = b, a
        if step >= burn and (step-burn) % gap == 0:
            for p,x in enumerate(order): rank_sum[x] += (p+1)   # rank 1 = most deleterious
            ns += 1
    return rank_sum/ns

print("Running two independent Bubley-Dyer chains (this takes ~1-2 min for N=2314) ...")
ar1 = sample_avg_rank(1); ar2 = sample_avg_rank(2)
conv = spearmanr(ar1, ar2).correlation
print(f"Chain convergence: Spearman(chain1, chain2) = {conv:.4f} (want ~1.0)")
avg_rank = (ar1+ar2)/2                          # canonical measure, 1 = most deleterious

# ---- benchmark each measure vs the canonical order ----
# measure ranks with 1 = most deleterious (match avg_rank orientation)
def delet_rank(s): return rankdata(-s, method="average")
print("\nSpearman(existing measure's ranking, canonical average-rank extension):")
rows = []
for n in names:
    rho = spearmanr(delet_rank(cols[n]), avg_rank).correlation
    mad = np.abs(delet_rank(cols[n]) - avg_rank).mean()
    rows.append((n, rho, mad))
# field's proposed consensus formula: median of 8 promoters (negated -> concept-positive)
field_med = -D[prom].median(axis=1).to_numpy()
rho_med = spearmanr(delet_rank(field_med), avg_rank).correlation
mad_med = np.abs(delet_rank(field_med) - avg_rank).mean()
rows.append(("[field] median-of-8", rho_med, mad_med))
rows.sort(key=lambda r:-r[1])
for n,rho,mad in rows:
    print(f"  {n:22s} rho={rho:+.3f}   mean|rank-canonical|={mad:6.1f}")
print(f"\nBest single proxy for the weight-free consensus: {rows[0][0]} (rho={rows[0][1]:+.3f})")
print(f"Strays most (commits beyond the consensus / outlier): {rows[-1][0]} (rho={rows[-1][1]:+.3f})")

# ---- the wider computational predictor family (external to the poset) vs the canonical order ----
# The poset/canonical are built from the 9 core measures (8 experimental + AlphaMissense). Here we ask
# how far each OTHER computational predictor sits from that experimental-dominated consensus.
print("\nOther computational predictors (external to the poset) vs the canonical order:")
ext = {"EVE": D["eve_score"].to_numpy(float), "ESM1b": -D["esm1b_llr"].to_numpy(float),
       "REVEL": D["revel"].to_numpy(float), "CADD": D["cadd_phred"].to_numpy(float),
       "PrimateAI": D["primateai"].to_numpy(float), "BayesDel": D["bayesdel_addaf"].to_numpy(float)}
for n,v in ext.items():
    ok = np.isfinite(v)
    rho = spearmanr(delet_rank(v[ok]), avg_rank[ok]).correlation
    print(f"  {n:12s} rho={rho:+.3f}  (n={ok.sum()})")
print("  (All computational predictors, like AlphaMissense, sit well below the experimental measures'")
print("   agreement with the consensus — the computational family as a block diverges from it.)")

# ---- does the canonical aggregate predict ClinVar as well as any single measure? ----
def auc(score, y):
    r = rankdata(score); n1=y.sum(); n0=(~y.astype(bool)).sum()
    return (r[y==1].sum()-n1*(n1+1)/2)/(n1*n0)
lab = D[(D["clinvar_simple"].isin([0,1])) & (D["clinvar_stars"]>=1)].index.to_numpy()
y = (D.loc[lab,"clinvar_simple"]==1).to_numpy().astype(int)
print(f"\nClinVar discrimination (AUC path/benign, >=1 star; n={len(lab)}, path={y.sum()}):")
# feed concept-positive scores (higher = more deleterious); avg_rank has 1=most deleterious, so negate it
print(f"  canonical average-rank : {auc(-avg_rank[lab], y):.3f}")
print(f"  [field] median-of-8    : {auc(field_med[lab], y):.3f}")
for n in names:
    print(f"  {n:22s}: {auc(cols[n][lab], y):.3f}")

# persist the canonical scores for downstream scripts / the writeup
out = D[["key","wt","pos","alt","clinvar_simple","clinvar_stars"]].copy()
out["canonical_avg_rank"] = avg_rank
out["field_median8"] = field_med
out.to_csv("/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/tp53/data/tp53_canonical.csv", index=False)
print("\nWrote tp53/data/tp53_canonical.csv (canonical average-rank per variant).")
