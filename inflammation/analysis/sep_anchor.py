"""
Sepsis signatures, protocol step 1/6 (orientation + anchor) and the endotype (clustering) caveat.

Anchor = 28-day mortality. We report each signature's RAW mortality AUC and its orientation sign (this
IS the orientation step; because we orient to the anchor, mortality is NOT re-used as an independent
validator — cross-cohort reproducibility plays that role, see sep_reproducibility.py). We then show the
mortality-oriented consensus across tertiles, and describe how the signatures relate to the given MARS
endotype labels — a partition, a different object from a scalar order (cluster comparison, not average
rank, is the right tool; noted honestly). Deterministic.
"""
import numpy as np, pandas as pd
from scipy.stats import rankdata, spearmanr

import sys as _sys
_OUT = "/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/inflammation/analysis/outputs/sepsis_anchor.txt"
class _Tee:
    def __init__(self, p): self._f = open(p, "w"); self._o = _sys.stdout
    def write(self, s): self._o.write(s); self._f.write(s)
    def flush(self): self._o.flush(); self._f.flush()
_sys.stdout = _Tee(_OUT)

SIG = ["Inflammatory","IFNg","SRS7","SRSq19","MARS8"]
def auc(score,y): r=rankdata(score);n1=y.sum();n0=(y==0).sum();return (r[y==1].sum()-n1*(n1+1)/2)/(n1*n0)
D = pd.read_csv("/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/inflammation/data/sepsis_scores.csv")
Dm = D[D["mortality"].notna()].reset_index(drop=True); y = Dm["mortality"].to_numpy().astype(int)
print(f"=== Sepsis anchor + endotype caveat (N={len(Dm)}, deaths={y.sum()}) ===\n")

print("(1) each signature vs mortality (raw AUC and the orientation sign it implies):")
signs={}
for s in SIG:
    a=auc(Dm[s].to_numpy(float),y); signs[s]=1 if a>=0.5 else -1
    print(f"  {s:12s} raw AUC={a:.3f}  -> orient {'+' if signs[s]>0 else '-'}  (|AUC-0.5|={abs(a-0.5):.3f})")
print("  (Signatures differ in how strongly they track outcome; some barely beat chance — they do not")
print("   agree on who is high-risk, the same disagreement seen in the family structure.)")

O = np.column_stack([signs[s]*Dm[s].to_numpy(float) for s in SIG])
cons = np.column_stack([rankdata(O[:,k]) for k in range(len(SIG))]).mean(1)
cpct=(rankdata(cons)-1)/(len(cons)-1); q1,q2=np.quantile(cpct,[1/3,2/3])
print(f"\n(2) mortality-oriented consensus: AUC={auc(cons,y):.3f}; death rate by consensus tertile:")
for lab,m in [("bottom",cpct<=q1),("middle",(cpct>q1)&(cpct<q2)),("top",cpct>=q2)]:
    print(f"    {lab:7s} n={m.sum():3d}  death rate={y[m].mean():.3f}")

print("\n(3) endotype caveat — signatures vs the given MARS labels (a PARTITION, not an order):")
De=D[D["endotype"].notna()].reset_index(drop=True)
Oe=np.column_stack([signs[s]*De[s].to_numpy(float) for s in SIG])
conse=np.column_stack([rankdata(Oe[:,k]) for k in range(len(SIG))]).mean(1)
for cls in ["Mars1","Mars2","Mars3","Mars4"]:
    m=(De["endotype"]==cls).to_numpy()
    print(f"    {cls}: n={m.sum():3d}  mean consensus-rank pct={((rankdata(conse)-1)/(len(conse)-1))[m].mean():.2f}")
print("    The consensus separates the MARS endotypes only partially; comparing two label systems (e.g.")
print("    SRS vs MARS) needs cluster-agreement (ARI), not the average-rank machinery — an adjacent problem.")
