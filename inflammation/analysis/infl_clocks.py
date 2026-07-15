"""
Topic 2 (inflammaging) — the thinnest instance, reported honestly.

Only TWO composite inflammatory-age clocks are publicly reconstructable on shared samples: SImAge
(per-sample values provided) and ipAGE (recomputed here as ElasticNet age-regression on the 46 cytokines,
fit on controls only, per Kalyakulina's method). iAge and IMM-AGE are not publicly obtainable. Two
measures are degenerate for the poset/near-tie machinery, so we widen the concept slightly to
"inflammatory burden" and add three canonical single-cytokine inflammaging markers used in the
literature: CXCL9 (the dominant iAge driver), IL6, TNF. All oriented concept-positive (higher = more
inflammatory burden / older). Anchor = chronological age. Data: SImAge cohort (n=343). Deterministic.

This is a smaller, more heterogeneous panel than the other two topics — treated as a qualitative check,
not a strong instance. See README limits.
"""
import numpy as np, pandas as pd
from scipy.stats import rankdata, spearmanr
from sklearn.linear_model import ElasticNetCV
from sklearn.preprocessing import StandardScaler

import sys as _sys
_OUT = "/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/inflammation/analysis/outputs/clocks.txt"
class _Tee:
    def __init__(self, p): self._f = open(p, "w"); self._o = _sys.stdout
    def write(self, s): self._o.write(s); self._f.write(s)
    def flush(self): self._o.flush(); self._f.flush()
_sys.stdout = _Tee(_OUT)

D = pd.read_excel("/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/inflammation/data/SImAge_data.xlsx")
BIO = list(D.columns[1:47])                                   # the 46 cytokine columns
age = D["Age"].to_numpy(float)
print(f"=== Inflammaging (SImAge cohort, n={len(D)}) ===")
print("Status:", D["Status"].value_counts().to_dict())

# recompute ipAGE: ElasticNet age-regression on 46 cytokines, fit on controls only
ctrl = D["Status"].astype(str).str.contains("Control", case=False)
sc = StandardScaler().fit(D.loc[ctrl, BIO])
en = ElasticNetCV(l1_ratio=[.1,.5,.7,.9,.95,.99,1], cv=5, max_iter=5000, random_state=0)
en.fit(sc.transform(D.loc[ctrl, BIO]), age[ctrl.to_numpy()])
ipage = en.predict(sc.transform(D[BIO]))
print(f"recomputed ipAGE: ElasticNet on {ctrl.sum()} controls, {np.sum(en.coef_!=0)}/46 nonzero coefs; "
      f"ipAGE vs age r={spearmanr(ipage,age).correlation:.3f}")

# concept-positive measures (higher = more inflammatory burden / older)
meas = {"SImAge": D["SImAge"].to_numpy(float), "ipAGE": ipage,
        "CXCL9": D["CXCL9"].to_numpy(float), "IL6": D["IL6"].to_numpy(float), "TNF": D["TNF"].to_numpy(float)}
names = list(meas); X = np.column_stack([meas[n] for n in names]); N=len(D)

print("\n(families) Spearman among the 5 measures:")
S = pd.DataFrame(X, columns=names).corr("spearman")
print("          " + " ".join(f"{n[:6]:>6s}" for n in names))
for n in names: print(f"  {n:8s} " + " ".join(f"{S.loc[n,c]:+.2f}" for c in names))
clk=["SImAge","ipAGE"]; cyt=["CXCL9","IL6","TNF"]
print(f"  clocks pair r={S.loc['SImAge','ipAGE']:+.2f}; within-cytokine mean "
      f"{S.loc[cyt,cyt].to_numpy()[np.triu_indices(3,1)].mean():+.2f}; "
      f"clocks-vs-cytokines mean {S.loc[clk,cyt].to_numpy().mean():+.2f}")

print("\n(anchor) each measure vs chronological age (Spearman):")
for n in names: print(f"  {n:8s} {spearmanr(meas[n],age).correlation:+.3f}")
cons = np.column_stack([rankdata(X[:,k]) for k in range(len(names))]).mean(1)
print(f"  consensus (mean rank) vs age: {spearmanr(cons,age).correlation:+.3f}")
cp=(rankdata(cons)-1)/(N-1); q1,q2=np.quantile(cp,[1/3,2/3])
print("  consensus vs age by tertile (mean age):",
      f"bottom {age[cp<=q1].mean():.1f} | middle {age[(cp>q1)&(cp<q2)].mean():.1f} | top {age[cp>=q2].mean():.1f}")

print("\n(near-tie) pairwise discordance among the 5 measures vs consensus separation:")
Rk=np.column_stack([rankdata(X[:,k]) for k in range(len(names))])
i,j=np.triu_indices(N,1); s0=np.sign(Rk[j,0]-Rk[i,0]); disc=np.zeros(len(i),bool)
for k in range(1,len(names)): disc|=(np.sign(Rk[j,k]-Rk[i,k])!=s0)
ax=cp; sep=np.abs(ax[i]-ax[j])
print(f"  overall discordance {disc.mean():.3f}")
for a,b in zip(np.quantile(sep,[0,.25,.5,.75])[:], np.quantile(sep,[.25,.5,.75,1.0])):
    m=(sep>=a)&(sep<=b); print(f"    sep [{a:.2f},{b:.2f}]: discordance {disc[m].mean():.3f}")

geq=np.ones((N,N),bool)
for k in range(len(names)): xk=X[:,k]; geq&=(xk[:,None]>=xk[None,:])
comp=geq|geq.T; np.fill_diagonal(comp,True)
print(f"\n(consensus poset) incomparable fraction = {1-(comp.sum()-N)/(N*(N-1)):.3f}  (5 measures, n={N})")
print("\nHonest reading: SImAge and ipAGE (both age-trained composites) agree strongly; the single")
print("cytokines induce different orders and agree less with the clocks — competing quantifications of")
print("inflammatory burden that disagree, as elsewhere. But this panel is small, heterogeneous, and only")
print("2 measures are true published clocks; treat as a qualitative check, not a strong instance.")
