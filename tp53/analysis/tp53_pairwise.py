"""
TP53 domain, protocol step 4 (locate disagreement, with controls) — the near-tie law.

Mirrors analysis/kinase_pairwise.py + kinase_concentration2.py: is inter-measure disagreement driven
by how CLOSE two variants are on the concept axis (near-ties), or by WHERE they sit (a 'consequential
middle')? Discordance panel = the 8 transactivation promoters. We test two separation axes:
  (a) the 8-promoter consensus (consensus coordinate; mild circularity, as in kinase_concentration)
  (b) AlphaMissense percentile (independent of the promoter panel, as in kinase_concentration2)
For each: bin pairs by |Δaxis| (separation) and, among near-ties, by middle-position; then a
linear-probability model asks whether middle-position adds anything beyond separation.
Deterministic.
"""
import numpy as np, pandas as pd
from scipy.stats import rankdata
from numpy.linalg import lstsq

# --- persist this script's console output to a tracked file (single source of truth for paper numbers) ---
import sys as _sys
_OUT = "/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/tp53/analysis/outputs/pairwise.txt"
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

# concept-positive promoter scores (higher = more deleterious = less transactivation)
P = np.column_stack([-D[p].to_numpy(float) for p in prom])       # (N,8)
Pr = np.column_stack([rankdata(P[:,k]) for k in range(8)])        # per-promoter ranks

i,j = np.triu_indices(N,1)
# discordant iff the 8 promoters do NOT unanimously agree which of the pair is more deleterious
s0 = np.sign(Pr[j,0]-Pr[i,0]); disc = np.zeros(len(i),bool)
for k in range(1,8):
    disc |= (np.sign(Pr[j,k]-Pr[i,k]) != s0)
print(f"=== TP53 pairwise discordance among 8 promoters (N={N}, {len(i):,} pairs) ===")
print(f"overall discordance rate: {disc.mean():.3f}\n")

def analyze(axis, label):
    ax = (rankdata(axis)-1)/(N-1)                 # 0..1 percentile on the concept axis
    sep = np.abs(ax[i]-ax[j])                      # separation (near-tie-ness)
    midpos = 1-2*np.abs((ax[i]+ax[j])/2 - 0.5)     # 1 = pair centered in the middle, 0 = at an end
    print(f"--- separation axis = {label} ---")
    print("  discordance by separation quintile (near-ties should dominate):")
    qs = np.quantile(sep,[0,.2,.4,.6,.8,1.0])
    for a,b in zip(qs[:-1],qs[1:]):
        m=(sep>=a)&(sep<=b if b==qs[-1] else sep<b)
        print(f"    |Δaxis| in [{a:.3f},{b:.3f}]: discordance={disc[m].mean():.3f}  (n_pairs={m.sum():,})")
    nt = sep <= np.quantile(sep,0.2)              # near-tie pairs (closest quintile)
    print("  among NEAR-TIE pairs, discordance by middle-position of the pair:")
    for lo,hi,lab in [(0,.33,"ends"),(.33,.66,"mid-ish"),(.66,1.01,"center")]:
        m=nt&(midpos>=lo)&(midpos<hi)
        print(f"    {lab:8s} (midpos {lo:.2f}-{hi:.2f}): discordance={disc[m].mean():.3f}  (n_pairs={m.sum():,})")
    y=disc.astype(float)
    X0=np.column_stack([sep,np.ones_like(sep)])
    X1=np.column_stack([sep,midpos,np.ones_like(sep)])
    b0=lstsq(X0,y,rcond=None)[0]; b1=lstsq(X1,y,rcond=None)[0]
    r0=1-((y-X0@b0).var()/y.var()); r1=1-((y-X1@b1).var()/y.var())
    print(f"  linear-prob model R^2: sep-only={r0:.4f}, sep+midpos={r1:.4f}; midpos coef={b1[1]:+.4f}")
    print(f"  => near-ties drive discordance; middle-position adds {'~0' if abs(b1[1])<0.03 else 'SOME'} beyond separation.\n")

analyze(P.mean(1), "8-promoter consensus (consensus coord)")
analyze(D["am_path"].to_numpy(float), "AlphaMissense (independent of the promoter panel)")
print("Conclusion (cf. kinase): separation (near-tie-ness) is the DOMINANT driver of TP53 disagreement")
print("(discordance ~0.93 at near-ties -> ~0.05 when well-separated). Unlike kinase, a MODEST residual")
print("middle effect survives: mid-deleteriousness near-ties disagree more than end near-ties. It is")
print("partly circular on the consensus axis (midpos coef +0.33) but shrinks to +0.10 on the independent")
print("AlphaMissense axis. Honest reading: the near-tie law holds and dominates, with a small genuine")
print("middle residual (partial-LOF / separation-of-function variants) that kinase lacked — a real")
print("cross-domain difference, not a 'consequential middle' of the retired P3 kind.")
