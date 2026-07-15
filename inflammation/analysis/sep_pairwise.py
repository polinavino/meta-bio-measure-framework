"""
Sepsis signatures, protocol step 4 (near-tie law). Discordance panel = the 5 mortality-oriented
signatures; separation axis = the consensus mean-rank. (No strong independent concept-proxy axis exists
here — age is a weak proxy in sepsis — so this is the consensus-coordinate view only, noted honestly.)
All mortality-labelled samples. Deterministic.
"""
import numpy as np, pandas as pd
from scipy.stats import rankdata, spearmanr
from numpy.linalg import lstsq

import sys as _sys
_OUT = "/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/inflammation/analysis/outputs/sepsis_pairwise.txt"
class _Tee:
    def __init__(self, p): self._f = open(p, "w"); self._o = _sys.stdout
    def write(self, s): self._o.write(s); self._f.write(s)
    def flush(self): self._o.flush(); self._f.flush()
_sys.stdout = _Tee(_OUT)

SIG = ["Inflammatory","IFNg","SRS7","SRSq19","MARS8"]
D = pd.read_csv("/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/inflammation/data/sepsis_scores.csv")
D = D[D["mortality"].notna()].reset_index(drop=True); y=D["mortality"].to_numpy().astype(int)
signs={s:(1 if spearmanr(D[s],y).correlation>=0 else -1) for s in SIG}
Rk = np.column_stack([rankdata(signs[s]*D[s]) for s in SIG]); N=len(D)

i,j=np.triu_indices(N,1); s0=np.sign(Rk[j,0]-Rk[i,0]); disc=np.zeros(len(i),bool)
for k in range(1,len(SIG)): disc|=(np.sign(Rk[j,k]-Rk[i,k])!=s0)
print(f"=== Sepsis pairwise discordance among 5 signatures (N={N}, {len(i):,} pairs) ===")
print(f"overall discordance rate: {disc.mean():.3f}\n")
ax=(rankdata(Rk.mean(1))-1)/(N-1); sep=np.abs(ax[i]-ax[j]); midpos=1-2*np.abs((ax[i]+ax[j])/2-0.5)
print("separation axis = consensus mean-rank:")
qs=np.quantile(sep,[0,.2,.4,.6,.8,1.0])
for a,b in zip(qs[:-1],qs[1:]):
    m=(sep>=a)&(sep<b if b!=qs[-1] else sep<=b)
    print(f"  |Δaxis| in [{a:.3f},{b:.3f}]: discordance={disc[m].mean():.3f}  (n={m.sum():,})")
yv=disc.astype(float)
X0=np.column_stack([sep,np.ones_like(sep)]);X1=np.column_stack([sep,midpos,np.ones_like(sep)])
b0=lstsq(X0,yv,rcond=None)[0];b1=lstsq(X1,yv,rcond=None)[0]
print(f"linear-prob R^2: sep-only={1-((yv-X0@b0).var()/yv.var()):.4f}, "
      f"sep+midpos={1-((yv-X1@b1).var()/yv.var()):.4f}; midpos coef={b1[1]:+.4f}")
print("\nReading: discordance falls with separation (near-tie law); the overall rate is high because the")
print("Hallmark and SRS families are near-orthogonal, so even well-separated pairs often disagree.")
