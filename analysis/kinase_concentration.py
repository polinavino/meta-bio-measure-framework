"""
Boundary-controlled concentration analysis for the kinase selectivity domain.

Referee challenge (empirical): the observed "disagreement concentrates in the middle"
may be (a) a mechanical rank-boundary compression artifact (ranks near 1 or n have less
room to spread), and (b) contradicted by zero-active compounds, which sit at an EXTREME
yet show maximal instability.

Design:
  1. Reproduce zero-active vs active instability (expect ~74 vs ~32).
  2. Resolve the extreme-counterexample principledly: zero-active compounds are exactly
     what the reliability gate D1/G1 excludes. Show the high-instability compounds ARE the
     ungated ones. Then test concentration only on the admissible (gated-in) set.
  3. Among gated-in compounds, test whether inter-MEASURE disagreement concentrates in the
     interior of the consensus selectivity coordinate, controlling for:
        - n_active and max activity (covariates),
        - mechanical rank-boundary compression, via an independence permutation null.
  4. Report honestly whether the interior peak survives the null + covariate control.

Data: klaeger_matrix.csv (222 x 343, pKd). Measure defs copied verbatim from the repo.
"""
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata

rng = np.random.default_rng(42)

M_df = pd.read_csv("/Users/polina/Documents/BioInfStuff/selectivity/klaeger_matrix.csv", index_col=0)
M = M_df.values.astype(float)
n_drugs, n_kin = M.shape
print(f"Klaeger matrix: {n_drugs} compounds x {n_kin} kinases")

# ---- measure definitions (verbatim from repo selectivity_analysis.py) ----
def s_score(P, threshold):        # Karaman fraction-hit; lower=selective
    return (P > threshold).astype(float).mean(axis=1)
def selectivity_entropy(P, baseline=5.0, eps=1e-10):  # higher=less selective
    sh = np.maximum(P - baseline, 0); rs = sh.sum(1, keepdims=True)
    rs = np.where(rs == 0, eps, rs); p = sh / rs
    lp = np.where(p > 0, np.log2(p + eps), 0); return -(p * lp).sum(1)
def gini_selectivity(P, baseline=5.0):  # higher=selective
    sh = np.maximum(P - baseline, 0); out = []
    for row in sh:
        rs = np.sort(row); n = len(rs); cs = np.cumsum(rs); tot = cs[-1]
        out.append(0.0 if tot == 0 else (2*np.sum(np.arange(1,n+1)*rs))/(n*tot)-(n+1)/n)
    return np.array(out)
def ratio_selectivity(P, top_n=1):  # higher=selective
    out = []
    for row in P:
        sd = np.sort(row)[::-1]; primary = sd[0]
        if len(sd) <= top_n: out.append(0.0); continue
        off = sd[top_n]; out.append(primary-5.0 if off <= 5.0 else primary-off)
    return np.array(out)

# ---- per-measure median ranking (rank 1 = most selective), matching repo ----
s_thr = np.arange(5.5, 8.25, 0.25)          # 11
base  = np.arange(5.0, 6.75, 0.25)          # 7
def to_rank(sel_score):  # higher sel_score -> rank 1
    return len(sel_score) - rankdata(sel_score, method="ordinal") + 1

meas_ranks = {}
meas_ranks["s_score"] = np.median([to_rank(-s_score(M,t)) for t in s_thr], axis=0)
meas_ranks["entropy"] = np.median([to_rank(-selectivity_entropy(M,b)) for b in base], axis=0)
meas_ranks["gini"]    = np.median([to_rank( gini_selectivity(M,b)) for b in base], axis=0)
meas_ranks["ratio"]   = np.median([to_rank( ratio_selectivity(M,n)) for n in range(1,6)], axis=0)
R = np.vstack([meas_ranks[k] for k in ["s_score","entropy","gini","ratio"]]).T  # (222,4)

# ---- covariates ----
n_active = (M > 6.0).sum(1)
max_act  = M.max(1)
gated_in = n_active > 0                      # D1/G1 reliability gate at tau*=6.0

# ---- disagreement per compound = std of its rank across the 4 measures ----
disagree = R.std(axis=1)

# ==== (1)+(2) zero-active vs active instability, and the gate resolution ====
print("\n[1/2] Reliability-gate resolution of the 'extreme counterexample'")
print(f"  zero-active (n_active==0): {np.sum(~gated_in)} compounds")
print(f"  inter-measure disagreement (std of rank across 4 measures):")
print(f"     ungated (zero-active): mean={disagree[~gated_in].mean():.1f}")
print(f"     gated-in   (active)  : mean={disagree[gated_in].mean():.1f}")
# also the repo's 30-config rank_std metric for the ~74 vs ~32 headline
all_scores = ([ -s_score(M,t) for t in s_thr] + [-selectivity_entropy(M,b) for b in base]
              + [gini_selectivity(M,b) for b in base] + [ratio_selectivity(M,n) for n in range(1,6)])
all_ranks = np.vstack([to_rank(s) for s in all_scores])   # (30,222)
rank_std30 = all_ranks.std(axis=0)
print(f"  repo 30-config rank_std: ungated={rank_std30[~gated_in].mean():.1f}, "
      f"gated-in={rank_std30[gated_in].mean():.1f}  (expect ~74 vs ~32)")

# ==== (3) concentration among gated-in, with null + covariate control ====
Rg = R[gated_in]; dis_g = disagree[gated_in]
na_g = n_active[gated_in]; mx_g = max_act[gated_in]
consensus = Rg.mean(axis=1)                       # consensus selectivity coordinate
coord = (rankdata(consensus)-1)/(len(consensus)-1)  # -> [0,1] percentile
middleness = 1 - 2*np.abs(coord - 0.5)            # 1 at center, 0 at extremes

print(f"\n[3] Concentration test on gated-in set (n={gated_in.sum()})")
# decile bins
bins = np.clip((coord*10).astype(int), 0, 9)
print("  decile | mean disagreement | n")
for b in range(10):
    m = bins == b
    print(f"    {b}    |   {dis_g[m].mean():6.2f}         | {m.sum()}")

# partial Spearman: disagreement ~ middleness, controlling n_active, max_act
def partial_spearman(x, y, controls):
    # residualize ranks of x and y on ranks of controls via OLS, then correlate
    from numpy.linalg import lstsq
    def resid(v):
        vr = rankdata(v)
        Z = np.column_stack([rankdata(c) for c in controls] + [np.ones_like(vr)])
        beta,_,_,_ = lstsq(Z, vr, rcond=None); return vr - Z@beta
    return spearmanr(resid(x), resid(y))
r_pc, p_pc = partial_spearman(middleness, dis_g, [na_g, mx_g])
r_raw, p_raw = spearmanr(middleness, dis_g)
print(f"\n  Spearman(disagreement, middleness): raw r={r_raw:.3f} (p={p_raw:.2e})")
print(f"  Partial Spearman controlling n_active,max_act: r={r_pc:.3f} (p={p_pc:.2e})")

# independence permutation null: destroy cross-measure agreement, keep each measure's marginal
# -> captures mechanical rank-boundary compression under NO concept structure
def disagree_vs_middle_null(nperm=1000):
    null_slope = []
    for _ in range(nperm):
        Rp = np.column_stack([rng.permutation(Rg[:,j]) for j in range(4)])
        dp = Rp.std(axis=1)
        cp = (rankdata(Rp.mean(1))-1)/(len(dp)-1); mp = 1-2*np.abs(cp-0.5)
        null_slope.append(spearmanr(mp, dp).correlation)
    return np.array(null_slope)
null = disagree_vs_middle_null()
obs = spearmanr(middleness, dis_g).correlation
print(f"\n  Interior-concentration slope (Spearman disagreement~middleness):")
print(f"     observed = {obs:.3f}")
print(f"     independence null: mean={null.mean():.3f}, 2.5-97.5% = [{np.percentile(null,2.5):.3f},{np.percentile(null,97.5):.3f}]")
print(f"     observed exceeds null: {obs > np.percentile(null,97.5)}")
print("\n  Interpretation: if observed >> null, interior concentration is structural,")
print("  not a boundary-compression artifact. If observed within null band, the")
print("  'middle' pattern is mechanical and the concentration claim must be downgraded.")
