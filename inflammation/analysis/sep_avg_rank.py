"""
Sepsis signatures, protocol step 7 (consensus poset + average-rank canonical). 5 mortality-oriented
signatures on the mortality-labelled GSE65682 samples. Bubley-Dyer sampler (same as ../../analysis/).
Benchmarks each signature against the canonical order and its mortality discrimination. Deterministic.
"""
import numpy as np, pandas as pd
from scipy.stats import rankdata, spearmanr

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
signs={s:(1 if spearmanr(D[s],y).correlation>=0 else -1) for s in SIG}
X = np.column_stack([signs[s]*D[s].to_numpy(float) for s in SIG]); N=len(D)

geq=np.ones((N,N),bool)
for k in range(X.shape[1]): xk=X[:,k]; geq &= (xk[:,None]>=xk[None,:])
comparable=geq|geq.T; np.fill_diagonal(comparable,True)
print(f"Consensus poset over N={N} samples (5 signatures): incomparable fraction = "
      f"{1-(comparable.sum()-N)/(N*(N-1)):.3f}")

def sar(seed,burn=300_000,nsamp=2000,gap=400):
    rng=np.random.default_rng(seed); order=list(np.argsort(-X.mean(1))); rs=np.zeros(N); ns=0
    for step in range(burn+nsamp*gap):
        i=rng.integers(N-1); a=order[i]; b=order[i+1]
        if not comparable[a,b] and rng.random()<0.5: order[i],order[i+1]=b,a
        if step>=burn and (step-burn)%gap==0:
            for p,x in enumerate(order): rs[x]+=(p+1)
            ns+=1
    return rs/ns
print("Running two Bubley-Dyer chains ...")
ar1=sar(1); ar2=sar(2); print(f"Chain convergence: Spearman = {spearmanr(ar1,ar2).correlation:.4f}")
avg=(ar1+ar2)/2
def dr(s): return rankdata(-s,method="average")
print("\nSpearman(signature ranking, canonical average-rank):")
for s,r in sorted([(SIG[k],spearmanr(dr(X[:,k]),avg).correlation) for k in range(len(SIG))],key=lambda t:-t[1]):
    print(f"  {s:12s} rho={r:+.3f}")
def auc(sc,yb): r=rankdata(sc);n1=yb.sum();n0=(yb==0).sum();return (r[yb==1].sum()-n1*(n1+1)/2)/(n1*n0)
print(f"\nMortality AUC:  canonical={auc(-avg,y):.3f}")
for k,s in enumerate(SIG): print(f"  {s:12s} {auc(X[:,k],y):.3f}")
