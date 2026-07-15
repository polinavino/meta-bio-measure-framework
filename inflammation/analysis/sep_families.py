"""
Sepsis signatures, protocol step 5 (families) + the headline disagreement.

5 competing signature scores on GSE65682. Each mean-z sign is not a-priori meaningful for the
mixed-direction classifier sets, so we orient every score to positively predict 28-day mortality
(the anchor) before comparing — see NOTE in build_sepsis_scores.py. Question: do the signatures cohere
or split into families? Spearman on the mortality-labelled samples. Deterministic.
"""
import numpy as np, pandas as pd
from scipy.stats import spearmanr

import sys as _sys
_OUT = "/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/inflammation/analysis/outputs/sepsis_families.txt"
class _Tee:
    def __init__(self, p): self._f = open(p, "w"); self._o = _sys.stdout
    def write(self, s): self._o.write(s); self._f.write(s)
    def flush(self): self._o.flush(); self._f.flush()
_sys.stdout = _Tee(_OUT)

SIG = ["Inflammatory","IFNg","SRS7","SRSq19","MARS8"]
D = pd.read_csv("/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/inflammation/data/sepsis_scores.csv")
D = D[D["mortality"].notna()].reset_index(drop=True)
y = D["mortality"].to_numpy().astype(int)

# orient each signature to positively predict mortality (record the sign it needed)
signs = {s: (1 if spearmanr(D[s], y).correlation >= 0 else -1) for s in SIG}
O = pd.DataFrame({s: signs[s]*D[s] for s in SIG})
print(f"=== Sepsis signature families (GSE65682, N={len(D)} mortality-labelled, deaths={y.sum()}) ===\n")
print("orientation sign applied so each positively predicts mortality:", signs)
print("(a negative sign means the raw mean-z ran opposite to outcome — mean-z direction is not a-priori")
print(" meaningful for the mixed-direction SRS/MARS classifier genes.)\n")

S = O[SIG].corr("spearman")
print("Spearman among mortality-oriented signature scores:")
print("            " + " ".join(f"{s[:7]:>7s}" for s in SIG))
for r in SIG:
    print(f"  {r:11s} " + " ".join(f"{S.loc[r,c]:+.2f}" for c in SIG))

hall=["Inflammatory","IFNg"]; srs=["SRS7","SRSq19"]
print(f"\nHallmark pair (Inflammatory,IFNg) r = {S.loc['Inflammatory','IFNg']:+.2f}")
print(f"SRS pair (SRS7,SRSq19)           r = {S.loc['SRS7','SRSq19']:+.2f}")
cross = S.loc[hall, srs].to_numpy().ravel()
print(f"Hallmark vs SRS families         r = {cross.min():+.2f}..{cross.max():+.2f} (mean {cross.mean():+.2f})")
print(f"MARS8 vs the rest                r = {S.loc['MARS8',[c for c in SIG if c!='MARS8']].min():+.2f}"
      f"..{S.loc['MARS8',[c for c in SIG if c!='MARS8']].max():+.2f}")
print("\nReading: even after orienting all five to the same outcome, the Hallmark-inflammatory and SRS")
print("families remain nearly independent — 'inflammatory dysregulation' is not one order. Strong")
print("instance of the framework's premise: competing measures of one concept genuinely disagree.")
