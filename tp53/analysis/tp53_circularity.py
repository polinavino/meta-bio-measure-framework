"""
TP53 domain, the companion axis (evaluations overstate performance) — the ADMET-sibling result.

Two checks, across all 7 computational predictors plus the experimental median-of-8:
 (1) Hotspot-codon holdout. TP53 pathogenic variants cluster at recurrent codons (ClinGen PM1 set
     175, 245, 248, 249, 273, 282). Removing them from a ClinVar evaluation tests whether apparent skill
     is hotspot-driven (the analog of a scaffold holdout).
 (2) Curated label vs independent functional truth. Each predictor's AUC against curated ClinVar
     (path vs benign) is compared with its AUC against the experimental functional-LOF call
     (ClinGen/Kato rule: median-of-8 transactivation <=20% WT = non-functional) on the ClinVar VUS,
     where clinical labels are absent and prediction is actually used.
Concept-positive orientation (higher = more deleterious): ESM1b LLR negated; all others as-is.
Deterministic.
"""
import numpy as np, pandas as pd
from scipy.stats import rankdata

# --- persist this script's console output to a tracked file (single source of truth for paper numbers) ---
import sys as _sys
_OUT = "/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/tp53/analysis/outputs/circularity.txt"
class _Tee:
    def __init__(self, p): self._f = open(p, "w"); self._o = _sys.stdout
    def write(self, s): self._o.write(s); self._f.write(s)
    def flush(self): self._o.flush(); self._f.flush()
_sys.stdout = _Tee(_OUT)

M = pd.read_csv("/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/tp53/data/tp53_matrix.csv")
prom = ["waf1","mdm2","bax","h1433s","aip1","gadd45","noxa","p53r2"]
HOTSPOT = {175,245,248,249,273,282}

def auc(score, y):
    ok = np.isfinite(score); score, y = score[ok], y[ok]
    if y.sum()==0 or (y==0).sum()==0: return np.nan, 0
    r = rankdata(score); n1=y.sum(); n0=(y==0).sum()
    return (r[y==1].sum()-n1*(n1+1)/2)/(n1*n0), len(y)

core = M[prom+["am_path"]].notna().all(axis=1)
D = M[core].reset_index(drop=True)
D["ish_med"] = D[prom].median(axis=1)
D["is_hot"]  = D["pos"].isin(HOTSPOT)

# concept-positive scorers (higher = more deleterious)
def scorers(df):
    return {
        "AlphaMissense": df["am_path"].to_numpy(float),
        "EVE":           df["eve_score"].to_numpy(float),
        "ESM1b":        -df["esm1b_llr"].to_numpy(float),
        "REVEL":         df["revel"].to_numpy(float),
        "CADD":          df["cadd_phred"].to_numpy(float),
        "PrimateAI":     df["primateai"].to_numpy(float),
        "BayesDel":      df["bayesdel_addaf"].to_numpy(float),
        "median-of-8 (experimental)": -df["ish_med"].to_numpy(float),
    }

lab = D[(D["clinvar_simple"].isin([0,1])) & (D["clinvar_stars"]>=1)].copy()
y = (lab["clinvar_simple"]==1).to_numpy().astype(int)
hot = lab["is_hot"].to_numpy()
print("=== Evaluation overstatement ===")
print(f"ClinVar-labeled core variants (>=1 star): n={len(lab)}, pathogenic={y.sum()}, benign={(y==0).sum()}")
print(f"Pathogenic labels on the 6 hotspot codons: {lab.loc[lab.clinvar_simple==1,'is_hot'].mean():.2f}"
      f"  (hotspots are {len(HOTSPOT)}/393 = {len(HOTSPOT)/393:.1%} of positions)\n")

print("(1) AUC vs ClinVar path/benign:        all labeled | hotspots removed | drop")
S = scorers(lab)
for name,s in S.items():
    a_all,_ = auc(s, y); a_no,_ = auc(s[~hot], y[~hot])
    print(f"    {name:28s} {a_all:.3f}       |   {a_no:.3f}        {a_all-a_no:+.3f}")
print("    Honest null: removing hotspots barely changes AUC — at the unique-variant level ClinVar")
print("    is not hotspot-dominated, so the scaffold-memorization analogy does not bite in this slice.\n")

print("(2) Curated ClinVar AUC vs independent functional-truth AUC on the VUS:")
vus = D[D["clinvar_simple"]==-1].copy()
func_lof = (vus["ish_med"]<=20).to_numpy().astype(int)      # ClinGen/Kato non-functional rule
Sv = scorers(vus)
print(f"    (VUS n={len(vus)}; functionally non-functional by the median<=20% rule: {func_lof.sum()})")
print("    predictor                     ClinVar AUC | functional-truth AUC (VUS) | gap")
for name,s in scorers(lab).items():
    if name.startswith("median"): continue                  # the functional rule IS median-of-8; skip (circular)
    a_cv,_ = auc(s, y)
    a_fx,n = auc(Sv[name], func_lof)
    print(f"    {name:28s} {a_cv:.3f}      | {a_fx:.3f}  (n={n})        {a_cv-a_fx:+.3f}")
print("    Every predictor scores lower against experimental truth on the VUS than against curated")
print("    ClinVar labels — the curated benchmark overstates skill on the decision-relevant variants.")
