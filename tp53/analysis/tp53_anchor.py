"""
TP53 domain, protocol step 1 (orientation vs an external anchor) + the anchor-validates-consensus check.

External anchor = ClinVar germline clinical significance (Pathogenic=1 / Benign=0), >=1 review star.
This is INDEPENDENT of the assay/predictor measures (clinical/pedigree evidence), so it is a
legitimate theory-light anchor for the concept "deleteriousness".

(1) ORIENTATION CHECK (G2). Each concept-positive measure (higher = more deleterious) should score
    pathogenic > benign, i.e. AUC > 0.5. We also test the DMS assays under their ProteinGym
    "higher = higher proliferative fitness" convention with the a-priori expectation that in these
    loss/dominance selection screens higher fitness = more deleterious (+), and let the anchor confirm
    or flag each — the operational content of step 1.

(2) ANCHOR VALIDATES THE CONSENSUS. Consensus = mean rank of the 9 core measures (8 promoters + AM).
    We check pathogenic-fraction rises monotonically across consensus tertiles, VUS concentrate in the
    middle tertile, and agreement with ClinVar is stronger at the consensus extremes than the middle
    (the concept decides the extremes).
Deterministic.
"""
import numpy as np, pandas as pd
from scipy.stats import rankdata

# --- persist this script's console output to a tracked file (single source of truth for paper numbers) ---
import sys as _sys
_OUT = "/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/tp53/analysis/outputs/anchor.txt"
class _Tee:
    def __init__(self, p): self._f = open(p, "w"); self._o = _sys.stdout
    def write(self, s): self._o.write(s); self._f.write(s)
    def flush(self): self._o.flush(); self._f.flush()
_sys.stdout = _Tee(_OUT)

M = pd.read_csv("/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/tp53/data/tp53_matrix.csv")
prom = ["waf1","mdm2","bax","h1433s","aip1","gadd45","noxa","p53r2"]

def auc(score, y):
    """P(score_pathogenic > score_benign); rank-based, ties averaged."""
    r = rankdata(score); n1 = y.sum(); n0 = (~y.astype(bool)).sum()
    return (r[y==1].sum() - n1*(n1+1)/2) / (n1*n0)

# concept-positive measures (higher = more deleterious)
def cp_measures(df):
    d = {f"TA:{p}": -df[p].to_numpy(float) for p in prom}      # transactivation: less = worse
    d["AlphaMissense"] = df["am_path"].to_numpy(float)
    # DMS under "higher fitness = more deleterious" hypothesis (to be checked):
    for c,lab in [("giac_wtnut","Giac:WT_Nutlin"),("giac_nullnut","Giac:Null_Nutlin"),
                  ("giac_nulletop","Giac:Null_Etop"),("kotler","Kotler")]:
        d[lab] = df[c].to_numpy(float)
    return d

# ---- (1) orientation check vs ClinVar (path vs benign, >=1 star) ----
lab = M[(M["clinvar_simple"].isin([0,1])) & (M["clinvar_stars"]>=1)].copy()
print("="*74)
print(f"(1) ORIENTATION CHECK vs ClinVar (path vs benign, >=1 star; n_labeled={len(lab)})")
print("    concept-positive measure should score pathogenic>benign  =>  AUC>0.5")
print("="*74)
meas = cp_measures(lab)
y = (lab["clinvar_simple"]==1).to_numpy().astype(int)
for name, mv in meas.items():
    ok = np.isfinite(mv)
    a = auc(mv[ok], y[ok])
    flag = "OK" if a>0.5 else "MIS-ORIENTED"
    print(f"  {name:16s} AUC={a:.3f}  (n={ok.sum():4d}, path={y[ok].sum()})   {flag}")
print("\n  (The 8 transactivation promoters + AlphaMissense should all read OK. The DMS null-selection")
print("   screens are the step-1 illustration: the anchor confirms/flags each screen's orientation.)")

# ---- (2) anchor validates the consensus ----
core = M[prom+["am_path"]].notna().all(axis=1).to_numpy()
D = M[core].copy()
cp = cp_measures(D)
core9 = [f"TA:{p}" for p in prom] + ["AlphaMissense"]
Rk = np.column_stack([rankdata(cp[n]) for n in core9])      # higher rank = more deleterious
consensus = Rk.mean(1)
D = D.assign(consensus=consensus)
cpct = (rankdata(consensus)-1)/(len(consensus)-1)
D = D.assign(cpct=cpct)

print("\n"+"="*74)
print(f"(2) ANCHOR VALIDATES THE CONSENSUS  (core panel N={core.sum()}; consensus = mean rank of 9)")
print("="*74)
# overall discrimination of the consensus
lab2 = D[(D["clinvar_simple"].isin([0,1])) & (D["clinvar_stars"]>=1)]
y2 = (lab2["clinvar_simple"]==1).to_numpy().astype(int)
print(f"  consensus AUC vs ClinVar path/benign: {auc(lab2['consensus'].to_numpy(),y2):.3f}"
      f"  (n={len(lab2)}, path={y2.sum()})")

# pathogenic fraction & VUS fraction across consensus tertiles
q1,q2 = np.quantile(cpct,[1/3,2/3])
D = D.assign(tert=np.where(cpct<=q1,"bottom",np.where(cpct>=q2,"top","middle")))
print("\n  Among ClinVar-labeled variants, by consensus tertile:")
print("   tertile   n_lab   %pathogenic   |   all core: %VUS(-1)   %unlabeled")
for t in ["bottom","middle","top"]:
    sub = D[D["tert"]==t]
    subl = sub[sub["clinvar_simple"].isin([0,1]) & (sub["clinvar_stars"]>=1)]
    pfrac = (subl["clinvar_simple"]==1).mean() if len(subl) else float("nan")
    vus = (sub["clinvar_simple"]==-1).mean()
    unl = sub["clinvar_simple"].isna().mean()
    print(f"   {t:7s}  {len(subl):5d}   {pfrac:10.2f}    |   {vus:12.2f}   {unl:9.2f}")
print("  (Expect: %pathogenic rises bottom->top; VUS/unlabeled concentrate in the middle tertile —")
print("   the concept decides the extremes and is silent in the middle.)")

# extremes vs middle agreement with the anchor
print("\n  Agreement (consensus AUC vs ClinVar) — extremes vs middle of the consensus axis:")
ext = (cpct<=q1)|(cpct>=q2); mid=~ext
for msk,lbl in [(ext,"extremes"),(mid,"middle")]:
    sub = D[msk & D["clinvar_simple"].isin([0,1]) & (D["clinvar_stars"]>=1)]
    ys = (sub["clinvar_simple"]==1).to_numpy().astype(int)
    a = auc(sub["consensus"].to_numpy(), ys) if ys.sum()>0 and (ys==0).sum()>0 else float("nan")
    print(f"    {lbl:8s}: consensus AUC {a:.3f}  (n={len(sub)}, path={ys.sum()}, benign={(ys==0).sum()})")
