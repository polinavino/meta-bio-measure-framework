"""
Clinical inflammatory indices, protocol step 5 (families by induced order).

7 competing indices of systemic inflammation, all concept-positive (higher = more inflammation):
  cell-count ratios: NLR, PLR, MLR, SII, SIRI ; CRP-protein: CRP, CAR (CRP/albumin).
Spearman (rank) correlation, so scale/skew-invariant. Question: do they split into a count-ratio family
and a CRP-protein family? Reports within/between-block agreement. NHANES 2015-16 complete panel.
Deterministic.
"""
import numpy as np, pandas as pd

# --- persist console output to a tracked file (single source of truth for README numbers) ---
import sys as _sys
_OUT = "/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/inflammation/analysis/outputs/families.txt"
class _Tee:
    def __init__(self, p): self._f = open(p, "w"); self._o = _sys.stdout
    def write(self, s): self._o.write(s); self._f.write(s)
    def flush(self): self._o.flush(); self._f.flush()
_sys.stdout = _Tee(_OUT)

IDX = ["NLR","PLR","MLR","SII","SIRI","CRP","CAR"]
COUNT = ["NLR","PLR","MLR","SII","SIRI"]; PROT = ["CRP","CAR"]
D = pd.read_csv("/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/inflammation/data/indices_2015.csv")
D = D[D[IDX].notna().all(axis=1)]
S = D[IDX].corr(method="spearman")

print(f"=== Inflammatory-index families (NHANES 2015-16, N={len(D)}) ===\n")
print("Spearman correlation (all concept-positive; higher = more inflammation):")
print("        " + " ".join(f"{c:>5s}" for c in IDX))
for r in IDX:
    print(f"  {r:5s} " + " ".join(f"{S.loc[r,c]:+.2f}" for c in IDX))

def blk(rows, cols, diag):
    v = S.loc[rows, cols].to_numpy()
    v = v[np.triu_indices_from(v,1)] if diag else v.ravel()
    return v.min(), v.max(), v.mean()
cc = blk(COUNT, COUNT, True); pp = blk(PROT, PROT, True); cp = blk(COUNT, PROT, False)
print(f"\nWithin count-ratio family (NLR/PLR/MLR/SII/SIRI): r = {cc[0]:.2f}..{cc[1]:.2f}, mean {cc[2]:.2f}")
print(f"Within CRP-protein family (CRP/CAR):              r = {pp[0]:.2f}..{pp[1]:.2f}, mean {pp[2]:.2f}")
print(f"Between the two families:                         r = {cp[0]:.2f}..{cp[1]:.2f}, mean {cp[2]:.2f}")
mean_off = (S.sum(1)-1)/(len(IDX)-1)
print("\nMean agreement of each index with the rest (ascending):")
for i in mean_off.sort_values().index:
    print(f"  {i:5s} {mean_off[i]:+.2f}")
print("\nReading: the count-ratio indices cohere; CRP/CAR sit apart (they measure an acute-phase protein,")
print("not a leukocyte ratio) — a family split by measurement basis, as in kinase and TP53.")
