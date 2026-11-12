"""
Sepsis signatures, protocol step 7 (consensus poset + average-rank canonical). 5 mortality-oriented
signatures on the mortality-labelled GSE65682 samples. checked insertion Gibbs sampler (same as ../../analysis/).
Benchmarks each signature against the canonical order and its mortality discrimination. Deterministic.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "analysis"))
from order_utils import pair_disagreement, average_rank
import numpy as np, pandas as pd
from scipy.stats import rankdata, spearmanr
from sepsis_orientation import discovery_signs

import sys as _sys
_OUT = "/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/inflammation/analysis/outputs/sepsis_avg_rank.txt"
class _Tee:
    def __init__(self, p): self._f = open(p, "w"); self._o = _sys.stdout
    def write(self, s): self._o.write(s); self._f.write(s)
    def flush(self): self._o.flush(); self._f.flush()
_sys.stdout = _Tee(_OUT)

SIG = ["Inflammatory","IFNg","SRS7","SRSq19","MARS8"]
D = pd.read_csv("/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/inflammation/data/sepsis_scores.csv")
D = D[D["mortality"].notna()].reset_index(drop=True); y=D["mortality"].to_numpy().astype(int)
signs = discovery_signs(D, SIG)
X = np.column_stack([signs[s]*D[s].to_numpy(float) for s in SIG]); N=len(D)

geq=np.ones((N,N),bool)
for k in range(X.shape[1]): xk=X[:,k]; geq &= (xk[:,None]>=xk[None,:])
comparable=geq|geq.T; np.fill_diagonal(comparable,True)
print(f"Consensus poset over N={N} samples (5 signatures): incomparable fraction = "
      f"{1-(comparable.sum()-N)/(N*(N-1)):.3f}")

avg = average_rank(X, "Sepsis consensus")
def dr(s): return rankdata(-s,method="average")
print("\nSpearman(signature ranking, canonical average-rank):")
for s,r in sorted([(SIG[k],spearmanr(dr(X[:,k]),avg).correlation) for k in range(len(SIG))],key=lambda t:-t[1]):
    print(f"  {s:12s} rho={r:+.3f}")
def auc(sc,yb): r=rankdata(sc);n1=yb.sum();n0=(yb==0).sum();return (r[yb==1].sum()-n1*(n1+1)/2)/(n1*n0)
print(f"\nMortality AUC:  canonical={auc(-avg,y):.3f}")
for k,s in enumerate(SIG): print(f"  {s:12s} {auc(X[:,k],y):.3f}")
