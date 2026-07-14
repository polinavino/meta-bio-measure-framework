"""
Cleaner re-analysis addressing: (i) gate-cutoff sensitivity, (ii) is 'where disagreement
concentrates' just near-ties + density, not a concept-level middle/low-signal law?

Pairwise discordance design (avoids rank-boundary compression and per-object density noise):
  axis = n_active (raw data property, external-ish selectivity proxy)
  for each pair: separation = |Δn_active|; position = mean n_active percentile;
  discordant = the 4 measures do NOT unanimously agree on the pair's selectivity order.
Question: after controlling for separation (near-tie-ness), does POSITION still predict
discordance? If not -> no concept-level concentration; disagreement = near-ties, wherever they sit.
"""
import numpy as np, pandas as pd
from scipy.stats import rankdata
M = pd.read_csv("/Users/polina/Documents/BioInfStuff/selectivity/klaeger_matrix.csv", index_col=0).values.astype(float)

def s_score(P,t): return (P>t).mean(1)
def ent(P,b=5.0,e=1e-10):
    sh=np.maximum(P-b,0); rs=sh.sum(1,keepdims=True); rs=np.where(rs==0,e,rs); p=sh/rs
    return -(p*np.where(p>0,np.log2(p+e),0)).sum(1)
def gini(P,b=5.0):
    sh=np.maximum(P-b,0); o=[]
    for r in sh:
        rs=np.sort(r);n=len(rs);tot=rs.sum(); o.append(0.0 if tot==0 else (2*np.sum(np.arange(1,n+1)*rs))/(n*tot)-(n+1)/n)
    return np.array(o)
def ratio(P,k=1):
    o=[]
    for r in P:
        sd=np.sort(r)[::-1]; o.append(sd[0]-5.0 if (len(sd)<=k or sd[k]<=5.0) else sd[0]-sd[k])
    return np.array(o)
def to_rank(s): return len(s)-rankdata(s,method="ordinal")+1
s_thr=np.arange(5.5,8.25,.25); base=np.arange(5.0,6.75,.25)

# ---- (i) gate-cutoff sensitivity: does zero-active-vs-active instability depend on pKd cutoff? ----
print("=== (i) Gate-cutoff sensitivity (active defined as pKd > cutoff) ===")
R_full=np.vstack([np.median([to_rank(-s_score(M,t)) for t in s_thr],0),
   np.median([to_rank(-ent(M,b)) for b in base],0),
   np.median([to_rank(gini(M,b)) for b in base],0),
   np.median([to_rank(ratio(M,k)) for k in range(1,6)],0)]).T
dis_all=R_full.std(1)
for cut in [5.5,6.0,6.5,7.0]:
    na=(M>cut).sum(1); z=na==0
    print(f"  cutoff pKd>{cut}: #zero-active={z.sum():3d} | disagreement zero-active={dis_all[z].mean():5.1f} vs active={dis_all[~z].mean():5.1f}")

# ---- (ii) pairwise discordance vs separation and position (gate: n_active>0 at pKd>6) ----
na=(M>6.0).sum(1); g=na>0
R=R_full[g]; nag=na[g]; n=g.sum()
pos_pct=(rankdata(nag)-1)/(n-1)     # 0=most selective(fewest active),1=most promiscuous
i,j=np.triu_indices(n,1)
# discordant: not all 4 measures agree on order (lower rank = more selective)
di=np.sign(R[j]-R[i])               # (npairs,4): +1 if i more selective on that measure
disc=~(np.all(di==di[:,[0]],axis=1))   # True if the 4 measures disagree on direction
sep=np.abs(nag[i]-nag[j])
midpos=1-2*np.abs((pos_pct[i]+pos_pct[j])/2 - 0.5)   # 1=pair centered in middle, 0=at an end

print(f"\n=== (ii) Pairwise discordance (gated-in n={n}, {len(i)} pairs) ===")
print("  discordance rate by n_active separation (near-ties should dominate):")
sep_bins=[(0,0),(1,2),(3,5),(6,15),(16,999)]
for lo,hi in sep_bins:
    m=(sep>=lo)&(sep<=hi)
    print(f"    |Δn_active| in [{lo:3d},{hi:3d}]: discordance={disc[m].mean():.3f}  (n_pairs={m.sum()})")

print("\n  Among NEAR-TIE pairs (|Δn_active|<=2), discordance by middle-position of the pair:")
nt=sep<=2
for lo,hi,lab in [(0,.33,"ends"),(.33,.66,"mid-ish"),(.66,1.01,"center")]:
    m=nt&(midpos>=lo)&(midpos<hi)
    print(f"    {lab:8s} (midpos {lo:.2f}-{hi:.2f}): discordance={disc[m].mean():.3f}  (n_pairs={m.sum()})")

# logistic-style check: does midpos add beyond separation? compare discordance ~ sep alone vs +midpos
from numpy.linalg import lstsq
X0=np.column_stack([np.log1p(sep),np.ones_like(sep,dtype=float)])
X1=np.column_stack([np.log1p(sep),midpos,np.ones_like(sep,dtype=float)])
y=disc.astype(float)
b0,_,_,_=lstsq(X0,y,rcond=None); b1,_,_,_=lstsq(X1,y,rcond=None)
r0=y-X0@b0; r1=y-X1@b1
print(f"\n  Linear-prob model: R^2 sep-only={1-r0.var()/y.var():.4f}, sep+midpos={1-r1.var()/y.var():.4f}")
print(f"  midpos coefficient (added value beyond separation): {b1[1]:+.4f}")
print("\n  If discordance is dominated by small |Δn_active| and midpos adds ~0 beyond separation,")
print("  then disagreement = near-ties (wherever they fall), NOT a concept-level 'middle'.")
