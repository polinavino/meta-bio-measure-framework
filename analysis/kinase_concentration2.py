"""
Robustness check: use an INDEPENDENT selectivity coordinate (n_active, a raw data
property) instead of the mean-of-4-measures, to avoid coordinate/disagreement circularity.
Test whether inter-measure disagreement peaks at intermediate n_active beyond a null that
holds n_active fixed and permutes measure identities.
"""
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata
rng = np.random.default_rng(7)

M = pd.read_csv("/Users/polina/Documents/BioInfStuff/selectivity/klaeger_matrix.csv", index_col=0).values.astype(float)

def s_score(P,t): return (P>t).mean(1)
def ent(P,b=5.0,e=1e-10):
    sh=np.maximum(P-b,0); rs=sh.sum(1,keepdims=True); rs=np.where(rs==0,e,rs); p=sh/rs
    return -(p*np.where(p>0,np.log2(p+e),0)).sum(1)
def gini(P,b=5.0):
    sh=np.maximum(P-b,0); o=[]
    for r in sh:
        rs=np.sort(r);n=len(rs);tot=rs.sum()
        o.append(0.0 if tot==0 else (2*np.sum(np.arange(1,n+1)*rs))/(n*tot)-(n+1)/n)
    return np.array(o)
def ratio(P,k=1):
    o=[]
    for r in P:
        sd=np.sort(r)[::-1]
        o.append(sd[0]-5.0 if (len(sd)<=k or sd[k]<=5.0) else sd[0]-sd[k])
    return np.array(o)
def to_rank(s): return len(s)-rankdata(s,method="ordinal")+1

s_thr=np.arange(5.5,8.25,.25); base=np.arange(5.0,6.75,.25)
R=np.vstack([
   np.median([to_rank(-s_score(M,t)) for t in s_thr],0),
   np.median([to_rank(-ent(M,b)) for b in base],0),
   np.median([to_rank(gini(M,b)) for b in base],0),
   np.median([to_rank(ratio(M,k)) for k in range(1,6)],0)]).T
n_active=(M>6.0).sum(1)
g=n_active>0
Rg=R[g]; na=n_active[g]; dis=Rg.std(1)

# bin by n_active quantiles
print("Disagreement vs n_active (independent coordinate), gated-in n=",g.sum())
qs=np.quantile(na,[0,.2,.4,.6,.8,1.0])
print(" n_active range | mean disagreement | n")
for i in range(5):
    lo,hi=qs[i],qs[i+1]; m=(na>=lo)&(na<=hi) if i==4 else (na>=lo)&(na<hi)
    print(f"   [{lo:4.0f},{hi:4.0f}]   |     {dis[m].mean():6.2f}      | {m.sum()}")

# middle-ness on n_active axis (log scale since n_active skewed)
la=np.log1p(na); mid=1-2*np.abs((rankdata(la)-1)/(len(la)-1)-0.5)
obs=spearmanr(mid,dis).correlation
# null: permute each measure's ranks independently (n_active fixed), recompute disagreement
nn=[]
for _ in range(1000):
    Rp=np.column_stack([rng.permutation(Rg[:,j]) for j in range(4)])
    nn.append(spearmanr(mid,Rp.std(1)).correlation)
nn=np.array(nn)
print(f"\n interior-peak slope vs n_active: observed={obs:.3f}")
print(f" null (measures independent, n_active fixed): mean={nn.mean():.3f}, "
      f"95% CI=[{np.percentile(nn,2.5):.3f},{np.percentile(nn,97.5):.3f}]")
print(f" observed exceeds null 97.5%: {obs>np.percentile(nn,97.5)}")
# also: do measures agree MORE at the selective extreme than the promiscuous extreme?
sel_ext=na<=np.quantile(na,.1); prom_ext=na>=np.quantile(na,.9)
print(f"\n disagreement at selective extreme (few active): {dis[sel_ext].mean():.2f}")
print(f" disagreement at promiscuous extreme (many active): {dis[prom_ext].mean():.2f}")
