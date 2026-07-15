"""
Clinical inflammatory indices, protocol step 3 (reproducibility across independent cohorts).

Two NHANES cycles (2015-16, 2017-18) are independent samples of the US population. We can't match people
across cycles, so (as in the smoking domain) we test whether the DISAGREEMENT STRUCTURE replicates:
  (a) the inter-index correlation matrix (family structure),
  (b) mean inter-index agreement,
  (c) the near-tie slope (discordance falls with separation),
  (d) the anchor orientation (each index vs mortality).
Structure that reproduces across independent cohorts is constitutive, not sampling noise. Deterministic.
"""
import numpy as np, pandas as pd
from scipy.stats import rankdata, spearmanr

# --- persist console output to a tracked file (single source of truth for README numbers) ---
import sys as _sys
_OUT = "/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/inflammation/analysis/outputs/reproducibility.txt"
class _Tee:
    def __init__(self, p): self._f = open(p, "w"); self._o = _sys.stdout
    def write(self, s): self._o.write(s); self._f.write(s)
    def flush(self): self._o.flush(); self._f.flush()
_sys.stdout = _Tee(_OUT)

IDX = ["NLR","PLR","MLR","SII","SIRI","CRP","CAR"]
DATA = "/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/inflammation/data"
def load(cyc):
    d = pd.read_csv(f"{DATA}/indices_{cyc}.csv")
    return d[d[IDX].notna().all(axis=1)].reset_index(drop=True)
A, B = load("2015"), load("2017")
def auc(score,y): r=rankdata(score);n1=y.sum();n0=(y==0).sum();return (r[y==1].sum()-n1*(n1+1)/2)/(n1*n0)

print(f"=== Cross-cohort reproducibility: NHANES 2015-16 (N={len(A)}) vs 2017-18 (N={len(B)}) ===\n")

SA, SB = A[IDX].corr("spearman"), B[IDX].corr("spearman")
iu = np.triu_indices(len(IDX),1)
va, vb = SA.to_numpy()[iu], SB.to_numpy()[iu]
print("(a) inter-index correlation structure:")
print(f"    Spearman(2015 corr-vector, 2017 corr-vector) = {spearmanr(va,vb).correlation:.3f}  "
      f"(21 index pairs); mean|Δr| = {np.abs(va-vb).mean():.3f}")
print(f"(b) mean inter-index agreement:  2015 {va.mean():.3f}  |  2017 {vb.mean():.3f}")

print("\n(c) near-tie slope (discordance by consensus-separation quintile), per cohort:")
def slope(D):
    D = D.iloc[np.random.default_rng(0).choice(len(D),size=min(2000,len(D)),replace=False)].reset_index(drop=True)
    N=len(D); Rk=np.column_stack([rankdata(D[c]) for c in IDX])
    i,j=np.triu_indices(N,1); s0=np.sign(Rk[j,0]-Rk[i,0]); disc=np.zeros(len(i),bool)
    for k in range(1,len(IDX)): disc|=(np.sign(Rk[j,k]-Rk[i,k])!=s0)
    ax=(rankdata(Rk.mean(1))-1)/(N-1); sep=np.abs(ax[i]-ax[j])
    qs=np.quantile(sep,[0,.2,.8,1.0])
    near=(sep>=qs[0])&(sep<qs[1]); far=(sep>=qs[2])
    return disc[near].mean(), disc[far].mean()
for cyc,D in [("2015",A),("2017",B)]:
    n,f = slope(D); print(f"    {cyc}: near-tie discordance {n:.3f}  ->  well-separated {f:.3f}")

print("\n(d) anchor orientation (index vs mortality AUC), per cohort:")
print("    index   2015    2017")
for c in IDX:
    a=A[A["mortstat"].notna()]; b=B[B["mortstat"].notna()]
    aa=auc(a[c].to_numpy(float),a["mortstat"].to_numpy().astype(int))
    bb=auc(b[c].to_numpy(float),b["mortstat"].to_numpy().astype(int))
    print(f"    {c:5s}  {aa:.3f}  {bb:.3f}")
print("\nReading: matching correlation structure, near-tie slopes, and anchor orientation across two")
print("independent cohorts => the disagreement structure is a property of the indices, not one sample.")
