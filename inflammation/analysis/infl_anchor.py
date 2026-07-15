"""
Clinical inflammatory indices, protocol step 1 (orientation vs external anchor) + anchor-validates-consensus.

External anchors (independent of the blood-count formulas): all-cause mortality (linked mortality file)
and chronological age. Higher inflammation is expected to track higher mortality and higher age.
(1) orientation: each concept-positive index should have mortality AUC > 0.5 and positive age correlation.
(2) consensus (mean rank of the 7 indices) validated: death rate rises across consensus tertiles; AUC vs
    mortality; extremes vs middle. NHANES 2015-16 complete panel with mortality. Deterministic.
"""
import numpy as np, pandas as pd
from scipy.stats import rankdata, spearmanr

# --- persist console output to a tracked file (single source of truth for README numbers) ---
import sys as _sys
_OUT = "/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/inflammation/analysis/outputs/anchor.txt"
class _Tee:
    def __init__(self, p): self._f = open(p, "w"); self._o = _sys.stdout
    def write(self, s): self._o.write(s); self._f.write(s)
    def flush(self): self._o.flush(); self._f.flush()
_sys.stdout = _Tee(_OUT)

IDX = ["NLR","PLR","MLR","SII","SIRI","CRP","CAR"]
def auc(score, y):
    r = rankdata(score); n1 = y.sum(); n0 = (y==0).sum()
    return (r[y==1].sum()-n1*(n1+1)/2)/(n1*n0)

D = pd.read_csv("/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/inflammation/data/indices_2015.csv")
D = D[D[IDX].notna().all(axis=1) & D["mortstat"].notna()].reset_index(drop=True)
y = D["mortstat"].to_numpy().astype(int)
age = D["age"].to_numpy(float)
print(f"=== Orientation + anchor validation (NHANES 2015-16, N={len(D)}, deaths={y.sum()}) ===\n")

print("(1) ORIENTATION — each index vs external anchors (expect mortality AUC>0.5, age r>0):")
for c in IDX:
    s = D[c].to_numpy(float)
    print(f"  {c:5s}  mortality AUC={auc(s,y):.3f}   age r={spearmanr(s,age).correlation:+.3f}   "
          f"{'OK' if auc(s,y)>0.5 else 'MIS-ORIENTED'}")

# consensus = mean rank of the 7 indices
R = np.column_stack([rankdata(D[c]) for c in IDX]); cons = R.mean(1)
cpct = (rankdata(cons)-1)/(len(cons)-1)
print("\n(2) ANCHOR VALIDATES THE CONSENSUS (consensus = mean rank of 7 indices):")
print(f"  consensus mortality AUC = {auc(cons,y):.3f}   |  consensus vs age r = {spearmanr(cons,age).correlation:+.3f}")
q1,q2 = np.quantile(cpct,[1/3,2/3])
print("\n  death rate by consensus tertile (expect monotone rise):")
for lab,m in [("bottom",cpct<=q1),("middle",(cpct>q1)&(cpct<q2)),("top",cpct>=q2)]:
    print(f"    {lab:7s} n={m.sum():4d}  death rate={y[m].mean():.3f}  mean age={age[m].mean():.1f}")
ext=(cpct<=q1)|(cpct>=q2); mid=~ext
print(f"\n  consensus mortality AUC:  extremes {auc(cons[ext],y[ext]):.3f} (n={ext.sum()})"
      f"  |  middle {auc(cons[mid],y[mid]):.3f} (n={mid.sum()})")
print("\n  (Age confounds inflammation and mortality; mean age per tertile is reported so the reader can")
print("   see the anchor is crude but directionally valid — the point is orientation, not a causal claim.)")
