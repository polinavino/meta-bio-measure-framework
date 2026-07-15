"""
Clinical inflammatory indices, protocol step 4 (locate disagreement, controlled) — the near-tie law.

Discordance panel = the 7 indices. Two separation axes: the consensus mean-rank (consensus coord) and
chronological age (external to the index formulas). Question: is inter-index disagreement driven by how
close two people are on the inflammation axis (near-ties), or by where they sit? Fixed random subsample
of the NHANES 2015-16 complete panel (pairwise blows up on all ~5000). Deterministic (seeded subsample).
"""
import numpy as np, pandas as pd
from scipy.stats import rankdata
from numpy.linalg import lstsq

# --- persist console output to a tracked file (single source of truth for README numbers) ---
import sys as _sys
_OUT = "/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/inflammation/analysis/outputs/pairwise.txt"
class _Tee:
    def __init__(self, p): self._f = open(p, "w"); self._o = _sys.stdout
    def write(self, s): self._o.write(s); self._f.write(s)
    def flush(self): self._o.flush(); self._f.flush()
_sys.stdout = _Tee(_OUT)

IDX = ["NLR","PLR","MLR","SII","SIRI","CRP","CAR"]
D = pd.read_csv("/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/inflammation/data/indices_2015.csv")
D = D[D[IDX].notna().all(axis=1)].reset_index(drop=True)
sub = np.random.default_rng(0).choice(len(D), size=2000, replace=False)
D = D.iloc[sub].reset_index(drop=True); N = len(D)
Rk = np.column_stack([rankdata(D[c]) for c in IDX])       # per-index ranks (concept-positive)

i,j = np.triu_indices(N,1)
s0 = np.sign(Rk[j,0]-Rk[i,0]); disc = np.zeros(len(i),bool)
for k in range(1,len(IDX)):
    disc |= (np.sign(Rk[j,k]-Rk[i,k]) != s0)
print(f"=== Pairwise discordance among 7 indices (N={N}, {len(i):,} pairs) ===")
print(f"overall discordance rate: {disc.mean():.3f}\n")

def analyze(axis, label):
    ax = (rankdata(axis)-1)/(N-1)
    sep = np.abs(ax[i]-ax[j]); midpos = 1-2*np.abs((ax[i]+ax[j])/2 - 0.5)
    print(f"--- separation axis = {label} ---")
    qs = np.quantile(sep,[0,.2,.4,.6,.8,1.0])
    for a,b in zip(qs[:-1],qs[1:]):
        m=(sep>=a)&(sep<b if b!=qs[-1] else sep<=b)
        print(f"    |Δaxis| in [{a:.3f},{b:.3f}]: discordance={disc[m].mean():.3f}  (n={m.sum():,})")
    y=disc.astype(float)
    X0=np.column_stack([sep,np.ones_like(sep)]); X1=np.column_stack([sep,midpos,np.ones_like(sep)])
    b0=lstsq(X0,y,rcond=None)[0]; b1=lstsq(X1,y,rcond=None)[0]
    r0=1-((y-X0@b0).var()/y.var()); r1=1-((y-X1@b1).var()/y.var())
    print(f"    linear-prob R^2: sep-only={r0:.4f}, sep+midpos={r1:.4f}; midpos coef={b1[1]:+.4f}\n")

analyze(Rk.mean(1), "consensus mean-rank (consensus coord)")
analyze(D["age"].to_numpy(float), "chronological age (external to the index formulas)")
print("Reading: if discordance falls steeply with separation and midpos adds little, disagreement is a")
print("near-tie phenomenon, not a 'consequential middle' — the same pattern as kinase and TP53.")
