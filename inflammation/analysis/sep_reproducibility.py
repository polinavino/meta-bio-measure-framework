"""
Sepsis signatures, protocol step 3 (reproducibility across independent cohorts). GSE65682 carries a
MARS discovery/validation split. We test whether the disagreement structure replicates: (a) inter-signature
correlation structure, (b) mortality-orientation signs, (c) near-tie slope. Gene standardization and orientation are discovery-frozen. Replication tests stability of the
chosen proxy scores, not their biological admissibility or a constitutive cause. Deterministic.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "analysis"))
from order_utils import pair_disagreement, average_rank
import numpy as np, pandas as pd
from scipy.stats import rankdata, spearmanr

import sys as _sys
_OUT = "/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/inflammation/analysis/outputs/sepsis_reproducibility.txt"
class _Tee:
    def __init__(self, p): self._f = open(p, "w"); self._o = _sys.stdout
    def write(self, s): self._o.write(s); self._f.write(s)
    def flush(self): self._o.flush(); self._f.flush()
_sys.stdout = _Tee(_OUT)

SIG = ["Inflammatory","IFNg","SRS7","SRSq19","MARS8"]
D = pd.read_csv("/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/inflammation/data/sepsis_scores.csv")
A = D[D["cohort"]=="discovery"].reset_index(drop=True)
B = D[D["cohort"]=="validation"].reset_index(drop=True)
print(f"=== Cross-cohort reproducibility: MARS discovery (N={len(A)}) vs validation (N={len(B)}) ===\n")

def signs(D):
    y=D["mortality"].to_numpy().astype(int)
    return {s:(1 if spearmanr(D[s],y).correlation>=0 else -1) for s in SIG}
sA,sB=signs(A),signs(B)
print("(a) mortality-orientation signs per cohort (should match):")
print("   ", {s:(sA[s],sB[s]) for s in SIG}, "->", "MATCH" if sA==sB else "DIFFER")

def corr(D,sg):
    O=pd.DataFrame({s:sg[s]*D[s] for s in SIG}); return O.corr("spearman")
SA,SB=corr(A,sA),corr(B,sA)  # freeze discovery signs on validation
iu=np.triu_indices(len(SIG),1)
va,vb=SA.to_numpy()[iu],SB.to_numpy()[iu]
print(f"\n(b) inter-signature correlation structure: Spearman(disc-vec, val-vec) = {spearmanr(va,vb).correlation:.3f}"
      f"  (10 pairs); mean|Δr| = {np.abs(va-vb).mean():.3f}")

print("\n(c) near-tie slope (discordance near-tie -> well-separated), per cohort:")
def slope(D,sg):
    Rk=np.column_stack([rankdata(sg[s]*D[s]) for s in SIG]); N=len(D)
    i,j=np.triu_indices(N,1)
    disc, tie_only = pair_disagreement(Rk, i, j)
    ax=(rankdata(Rk.mean(1))-1)/(N-1); sep=np.abs(ax[i]-ax[j]); qs=np.quantile(sep,[0,.2,.8,1.0])
    return disc[(sep>=qs[0])&(sep<qs[1])].mean(), disc[sep>=qs[2]].mean()
for lab,D,sg in [("discovery",A,sA),("validation",B,sA)]:
    n,f=slope(D,sg); print(f"    {lab:11s} near-tie {n:.3f} -> well-separated {f:.3f}")
print("\nOrientation is fitted on discovery and frozen on validation; sign flips are reported above.")
print("Replication here concerns the chosen scores and samples. Shared bias can replicate;")
print("neither a stable curve nor a stable correlation matrix proves concept-level incomparability.")
