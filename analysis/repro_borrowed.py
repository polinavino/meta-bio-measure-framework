"""Reproduce/confirm the source-repo-carried numbers: serotonin S-score vs Gini, and the clock
figures (position-clock R2, coefficient cosine, cell-type |r|) from the clock repo's computed outputs."""
import pandas as pd, numpy as np
from scipy.stats import spearmanr

print("========== SEROTONIN: S-score vs Gini (cited -0.682) ==========")
df = pd.read_csv("/Users/polina/Documents/BioInfStuff/psychedelic-selectivity/data/selectivity_framework_results.csv")
print("rows:", len(df), "| n_receptors range:", df.n_receptors.min(), "-", df.n_receptors.max())
for label, sub in [("all rows", df),
                   ("n_receptors>=4", df[df.n_receptors>=4]),
                   ("n_receptors>=5", df[df.n_receptors>=5])]:
    if len(sub) > 3:
        rho = spearmanr(sub.s_score, sub.gini).correlation
        rp  = np.corrcoef(sub.s_score, sub.gini)[0,1]
        print(f"  {label:16s} n={len(sub):5d}  Spearman(s_score,gini)={rho:+.3f}  Pearson={rp:+.3f}")

print("\n========== CLOCKS: from the repo's computed analysis parquets ==========")
D = "/Users/polina/Documents/BioInfStuff/epigenetic-clock-desiderata/data/"
fr = pd.read_parquet(D+"functional_relationships.parquet")
print("\n-- functional_relationships (position-clock linear R2; rate-position log) --")
print(fr[["clock1","clock2","best_form","R2_linear","R2_log"]].to_string(index=False))
pos = fr[fr.clock1.isin(["Horvath","Hannum","PhenoAge"]) & fr.clock2.isin(["Horvath","Hannum","PhenoAge"])]
print(f"  position-clock pair R2_linear range: {pos.R2_linear.min():.3f}-{pos.R2_linear.max():.3f}  (cited 0.991-0.995)")

cs = pd.read_parquet(D+"cosine_similarity.parquet")
print("\n-- cosine_similarity (coefficient-vector cosines) --")
print(cs.to_string())
vals = cs.values.astype(float); off = vals[~np.eye(len(vals),dtype=bool)]
print(f"  off-diagonal cosine range: {off.min():.3f}-{off.max():.3f}  (cited 0.03-0.11; Horvath-DunedinPACE=0)")

d3 = pd.read_parquet(D+"d3_cell_type.parquet")
print("\n-- d3_cell_type (max correlation with cell-type proportions) --")
print(d3.to_string(index=False))
print(f"  max cell-type |r| across clocks: {d3.max_cell_cor.max():.3f}  (cited up to 0.372)")
