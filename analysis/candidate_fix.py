"""
Candidate FIX: pedestal-subtracted softplus hinge.

Problem (found): w_i = T*softplus((x_i-floor)/T) has a pedestal T*ln2 at x=floor that never
vanishes, so on floor-dominated panels (93.6% of Klaeger at floor) the inactive kinases dominate
the weight distribution and the measure inverts (selective compounds look uniform -> high entropy).

Fix: anchor the hinge at the floor by subtracting the pedestal:
     w_i = T * ( softplus((x_i - floor)/T) - softplus(0) )   [ = 0 at x_i = floor, smooth, monotone ]
For floored data (x_i >= floor) this is >= 0, so no clamp needed. Keeps D1 (gate), D2 (distributional),
D3 (smooth, no hard cutoff), D4 (sub-floor adds ~0). Verify all + un-inversion + panel convergence.
"""
import numpy as np, pandas as pd
from scipy.stats import spearmanr
M = pd.read_csv("/Users/polina/Documents/BioInfStuff/selectivity/klaeger_matrix.csv", index_col=0).values.astype(float)
n_drugs, n_kin = M.shape
FLOOR, TAU_STAR, T = 5.0, 6.0, 1.0
def to_ranks(s): return len(s) - s.argsort().argsort()

# ORIGINAL (broken) and FIXED candidates
def candidate_orig(P, floor=FLOOR, T=T, eps=1e-10):
    w = T*np.logaddexp(0.0,(P-floor)/T)
    rs=np.where(w.sum(1,keepdims=True)==0,eps,w.sum(1,keepdims=True)); p=w/rs
    return -(-(p*np.where(p>0,np.log2(p+eps),0)).sum(1))
def candidate_fixed(P, floor=FLOOR, T=T, eps=1e-10):
    w = T*(np.logaddexp(0.0,(P-floor)/T) - np.log(2.0))     # pedestal-subtracted; 0 at x=floor
    w = np.maximum(w, 0.0)                                    # safety for any x<floor
    rs=np.where(w.sum(1,keepdims=True)==0,eps,w.sum(1,keepdims=True)); p=w/rs
    return -(-(p*np.where(p>0,np.log2(p+eps),0)).sum(1))
def entropy(P,b=5.0,e=1e-10):
    s=np.maximum(P-b,0); rs=np.where(s.sum(1,keepdims=True)==0,e,s.sum(1,keepdims=True)); p=s/rs
    return -(-(p*np.where(p>0,np.log2(p+e),0)).sum(1))
def gini(P,b=5.0):
    s=np.maximum(P-b,0); o=[]
    for r in s:
        rs=np.sort(r);n=len(rs);t=rs.sum(); o.append(0.0 if t==0 else (2*np.sum(np.arange(1,n+1)*rs))/(n*t)-(n+1)/n)
    return np.array(o)
def s_score(P,thr=6.0): return -(P>thr).astype(float).mean(1)

print("=== (1) Un-inversion: Spearman(candidate, entropy/gini/s_score) ===")
for nm,fn in [("ORIGINAL",candidate_orig),("FIXED",candidate_fixed)]:
    c=to_ranks(fn(M))
    print(f"  {nm:8s}: vs entropy {spearmanr(c,to_ranks(entropy(M))).correlation:+.3f} | "
          f"vs gini {spearmanr(c,to_ranks(gini(M))).correlation:+.3f} | "
          f"vs s_score {spearmanr(c,to_ranks(s_score(M))).correlation:+.3f}")

print("\n=== (2) D3 baseline/floor robustness (worst-case pairwise Spearman, floor in [4.5,6.5]) ===")
betas=np.arange(4.5,6.6,0.25)
def worst(fn):
    rr=np.array([to_ranks(fn(M,floor=b)) for b in betas])
    return min(spearmanr(rr[i],rr[j]).correlation for i in range(len(betas)) for j in range(i+1,len(betas)))
def worst_hard():
    rr=np.array([to_ranks(entropy(M,b=b)) for b in betas])
    return min(spearmanr(rr[i],rr[j]).correlation for i in range(len(betas)) for j in range(i+1,len(betas)))
print(f"  FIXED candidate (smooth): {worst(candidate_fixed):+.3f}")
print(f"  hard entropy (baseline) : {worst_hard():+.3f}   (candidate should be >> this, ~0.99)")

print("\n=== (3) D4 monotonicity: add a sub-threshold off-target ===")
for nm,fn in [("ORIGINAL",candidate_orig),("FIXED",candidate_fixed)]:
    base=fn(M); new=fn(np.hstack([M,np.full((n_drugs,1),TAU_STAR-0.7)]))
    print(f"  {nm:8s}: non-increasing for {(new<=base+1e-9).mean()*100:.0f}% ; max increase {max(0.0,(new-base).max()):.4f}")

print("\n=== (4) Panel-size convergence p* (smallest panel, mean Spearman>0.90 vs full) ===")
panel=list(range(50,n_kin,30))+[n_kin]; np.random.seed(42); R=50
for nm,fn in [("FIXED cand",candidate_fixed),("entropy",entropy),("s_score",s_score),("gini",gini)]:
    ref=to_ranks(fn(M)); pstar=None
    for ps in panel:
        vals=[]
        for _ in range(R):
            idx=np.random.choice(n_kin,ps,replace=False)
            vals.append(spearmanr(ref,to_ranks(fn(M[:,idx]))).correlation)
        if np.mean(vals)>0.90: pstar=ps; break
    print(f"  {nm:10s} p* = {pstar}")

print("\n=== (5) FIXED candidate vs average-rank canonical (gated-in) — should now be a GOOD proxy ===")
# quick canonical via mean of the 4 established measures' ranks on gated-in (proxy for the sampled avg-rank)
na=(M>6.0).sum(1); g=na>0
from scipy.stats import rankdata
def rk(s): return len(s)-rankdata(s,method="ordinal")+1
cons=np.mean([rk(entropy(M))[g],rk(gini(M))[g],rk(s_score(M))[g]],axis=0)
print(f"  Spearman(FIXED candidate, consensus-of-3 distribution measures) = "
      f"{spearmanr(rk(candidate_fixed(M))[g], cons).correlation:+.3f}")
print(f"  (original was ~ -0.86)")
