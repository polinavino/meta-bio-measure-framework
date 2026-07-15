"""
TP53 domain, protocol step 5 (cluster measures into families by the order they induce).

Two candidate families of measures of one concept (deleteriousness of a TP53 missense variant),
all oriented concept-positive (higher = MORE deleterious):
  EXPERIMENTAL (8): Kato/Ishioka yeast transactivation promoters (negated %WT).
  COMPUTATIONAL (7): AlphaMissense, EVE, ESM1b (negated LLR), REVEL, CADD, PrimateAI, BayesDel.
Pairwise Spearman uses all shared variants per pair (pandas pairwise-complete).

Question: do the measures split into an experimental block and a computational block, each internally
coherent and the two diverging from one another? Reports within-block and between-block agreement.
Deterministic.
"""
import numpy as np, pandas as pd

# --- persist this script's console output to a tracked file (single source of truth for paper numbers) ---
import sys as _sys
_OUT = "/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/tp53/analysis/outputs/families.txt"
class _Tee:
    def __init__(self, p): self._f = open(p, "w"); self._o = _sys.stdout
    def write(self, s): self._o.write(s); self._f.write(s)
    def flush(self): self._o.flush(); self._f.flush()
_sys.stdout = _Tee(_OUT)

M = pd.read_csv("/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/tp53/data/tp53_matrix.csv")
prom = ["waf1","mdm2","bax","h1433s","aip1","gadd45","noxa","p53r2"]
# restrict to the core panel (all 8 promoters + AlphaMissense present)
core = M[prom+["am_path"]].notna().all(axis=1)
D = M[core]

# concept-positive columns (higher = more deleterious)
cp = pd.DataFrame(index=D.index)
for p in prom: cp[f"TA:{p}"] = -D[p]
cp["AlphaMissense"] = D["am_path"]
cp["EVE"]           = D["eve_score"]
cp["ESM1b"]         = -D["esm1b_llr"]          # LLR more negative = more damaging -> negate
cp["REVEL"]         = D["revel"]
cp["CADD"]          = D["cadd_phred"]
cp["PrimateAI"]     = D["primateai"]
cp["BayesDel"]      = D["bayesdel_addaf"]

exp_cols  = [f"TA:{p}" for p in prom]
comp_cols = ["AlphaMissense","EVE","ESM1b","REVEL","CADD","PrimateAI","BayesDel"]
names = exp_cols + comp_cols

S = cp[names].corr(method="spearman")          # pairwise-complete Spearman
print(f"=== TP53 measure families (core panel N={core.sum()}; pairwise-complete Spearman) ===\n")
print("Rank-correlation matrix (concept-positive; higher = more deleterious):")
print("            " + " ".join(f"{n[:7]:>7s}" for n in names))
for n in names:
    print(f"  {n:11s} " + " ".join(f"{S.loc[n,c]:+7.2f}" for c in names))

def block(rows, cols, diag):
    v = S.loc[rows, cols].to_numpy()
    v = v[np.triu_indices_from(v,1)] if diag else v.ravel()
    v = v[~np.isnan(v)]
    return v.min(), v.max(), v.mean()

ee = block(exp_cols,  exp_cols,  True)
cc = block(comp_cols, comp_cols, True)
ec = block(exp_cols,  comp_cols, False)
print(f"\nWithin EXPERIMENTAL (8 promoters, 28 pairs):   r = {ee[0]:.2f}..{ee[1]:.2f}, mean {ee[2]:.2f}")
print(f"Within COMPUTATIONAL (7 predictors, 21 pairs): r = {cc[0]:.2f}..{cc[1]:.2f}, mean {cc[2]:.2f}")
print(f"BETWEEN experimental and computational:        r = {ec[0]:.2f}..{ec[1]:.2f}, mean {ec[2]:.2f}")

# which computational predictor tracks the experimental block best / worst
mean_to_exp = {c: np.nanmean(S.loc[exp_cols, c].to_numpy()) for c in comp_cols}
order = sorted(mean_to_exp, key=mean_to_exp.get)
print("\nMean agreement of each computational predictor with the 8 experimental readouts:")
for c in order:
    print(f"  {c:12s} {mean_to_exp[c]:+.2f}")
print(f"\nClosest to experiment: {order[-1]} ({mean_to_exp[order[-1]]:+.2f}); "
      f"farthest: {order[0]} ({mean_to_exp[order[0]]:+.2f}).")
print("Reading: each block is internally coherent and the two blocks agree less with each other than")
print("within themselves (between-mean below both within-means) — a family split by measurement basis")
print("(experimental function vs computational prediction), the TP53 analog of kinase's family structure.")
