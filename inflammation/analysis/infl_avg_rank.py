"""
Clinical inflammatory indices, protocol step 7 (consensus poset + average-rank canonical aggregate).

7 concept-positive indices (higher = more inflammation). Consensus poset: a dominates b iff every index
scores a >= b. Canonical measure = mean rank over uniformly-sampled linear extensions (Bubley-Dyer;
same sampler as ../../analysis/avg_rank.py and ../../tp53/). Computed on a fixed random subsample of the
NHANES 2015-16 complete panel (the O(N^2) poset is impractical on all ~5000). Benchmarks each index and
the equal-weight mean-rank score against the canonical, and checks mortality discrimination.
Deterministic (seeded subsample + chains).
"""
import numpy as np, pandas as pd
from scipy.stats import rankdata, spearmanr

# --- persist console output to a tracked file (single source of truth for README numbers) ---
import sys as _sys
_OUT = "/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/inflammation/analysis/outputs/avg_rank.txt"
class _Tee:
    def __init__(self, p): self._f = open(p, "w"); self._o = _sys.stdout
    def write(self, s): self._o.write(s); self._f.write(s)
    def flush(self): self._o.flush(); self._f.flush()
_sys.stdout = _Tee(_OUT)

IDX = ["NLR","PLR","MLR","SII","SIRI","CRP","CAR"]
D = pd.read_csv("/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/inflammation/data/indices_2015.csv")
D = D[D[IDX].notna().all(axis=1) & D["mortstat"].notna()].reset_index(drop=True)
sub = np.random.default_rng(0).choice(len(D), size=600, replace=False)
D = D.iloc[sub].reset_index(drop=True); N = len(D)
X = np.column_stack([D[c].to_numpy(float) for c in IDX])   # concept-positive (higher = more inflammation)

geq = np.ones((N,N), bool)
for k in range(X.shape[1]):
    xk = X[:,k]; geq &= (xk[:,None] >= xk[None,:])
comparable = geq | geq.T; np.fill_diagonal(comparable, True)
frac_incomp = 1 - (comparable.sum()-N)/(N*(N-1))
print(f"Consensus poset over N={N} NHANES participants (7 indices): incomparable fraction = {frac_incomp:.3f}")

def sample_avg_rank(seed, burn=400_000, nsamp=2000, gap=500):
    rng = np.random.default_rng(seed)
    order = list(np.argsort(-X.mean(1)))
    rank_sum = np.zeros(N); ns = 0
    for step in range(burn + nsamp*gap):
        i = rng.integers(N-1); a=order[i]; b=order[i+1]
        if not comparable[a,b] and rng.random() < 0.5:
            order[i],order[i+1] = b,a
        if step>=burn and (step-burn)%gap==0:
            for p,x in enumerate(order): rank_sum[x]+=(p+1)
            ns+=1
    return rank_sum/ns
print("Running two Bubley-Dyer chains ...")
ar1=sample_avg_rank(1); ar2=sample_avg_rank(2)
print(f"Chain convergence: Spearman = {spearmanr(ar1,ar2).correlation:.4f}")
avg_rank = (ar1+ar2)/2                                     # 1 = most inflamed

def drank(s): return rankdata(-s, method="average")        # 1 = most inflamed
print("\nSpearman(index ranking, canonical average-rank):")
rows=[(c, spearmanr(drank(X[:,k]),avg_rank).correlation) for k,c in enumerate(IDX)]
mean_rank = rankdata(-rankdata(X,axis=0).mean(1))          # equal-weight mean-rank score
rows.append(("[equal-weight mean rank]", spearmanr(mean_rank,avg_rank).correlation))
for c,r in sorted(rows,key=lambda t:-t[1]): print(f"  {c:24s} rho={r:+.3f}")

# mortality discrimination of canonical vs each index
def auc(score,yb): r=rankdata(score);n1=yb.sum();n0=(yb==0).sum();return (r[yb==1].sum()-n1*(n1+1)/2)/(n1*n0)
y=D["mortstat"].to_numpy().astype(int)
print(f"\nMortality AUC (n={N}, deaths={y.sum()}):  canonical={auc(-avg_rank,y):.3f}")
for k,c in enumerate(IDX): print(f"  {c:5s} {auc(X[:,k],y):.3f}")
