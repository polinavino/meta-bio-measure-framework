"""
Bounded attempt at a candidate that satisfies ALL of D1-D4:
baseline-FREE (per-compound-anchored) concentration measures. A per-compound anchor is invariant to a
global additive shift of the profile, so D3 (baseline invariance) is satisfied by construction (there is
no floor parameter to vary). Question: can orientation (agreement with the consensus selectivity order)
be STRONG (>0.9) while keeping D4 monotonicity and fast panel convergence?
"""
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata
M = pd.read_csv("/Users/polina/Documents/BioInfStuff/selectivity/klaeger_matrix.csv", index_col=0).values.astype(float)
n_drugs,n_kin=M.shape

def s_score(P,t): return (P>t).mean(1)
def ent(P,b=5.0,e=1e-10):
    sh=np.maximum(P-b,0); rs=np.where(sh.sum(1,keepdims=True)==0,e,sh.sum(1,keepdims=True)); p=sh/rs
    return -(-(p*np.where(p>0,np.log2(p+e),0)).sum(1))
def gini(P,b=5.0):
    sh=np.maximum(P-b,0);o=[]
    for r in sh:
        rs=np.sort(r);n=len(rs);t=rs.sum();o.append(0.0 if t==0 else (2*np.sum(np.arange(1,n+1)*rs))/(n*t)-(n+1)/n)
    return np.array(o)
def rk(s): return len(s)-rankdata(s,method="ordinal")+1
# consensus selectivity order = mean rank of the 3 established distribution measures
na=(M>6.0).sum(1); g=na>0
cons_g = np.mean([rk(-s_score(M,6.0))[g], rk(ent(M))[g], rk(gini(M))[g]],axis=0)

def anchored_entropy(P, anchor="min", T=1.0, e=1e-10):
    if anchor=="min":  q=P.min(1,keepdims=True)
    elif anchor=="q10":q=np.quantile(P,0.10,axis=1,keepdims=True)
    elif anchor=="q25":q=np.quantile(P,0.25,axis=1,keepdims=True)
    elif anchor=="med":q=np.median(P,axis=1,keepdims=True)
    w=np.maximum(T*(np.logaddexp(0.0,(P-q)/T)-np.log(2.0)),0.0)
    rs=np.where(w.sum(1,keepdims=True)==0,e,w.sum(1,keepdims=True)); p=w/rs
    return -(-(p*np.where(p>0,np.log2(p+e),0)).sum(1))    # higher = more selective (neg entropy)

print("anchor  T   orient(vs consensus, gated)  D3(global shift)  D4(%non-incr)  p*")
panel=list(range(50,n_kin,30))+[n_kin]
for anchor in ["min","q10","q25","med"]:
    for T in [0.5,1.0,2.0]:
        val=anchored_entropy(M,anchor,T)
        orient=spearmanr(rk(val)[g], cons_g).correlation
        # D3 as invariance to global additive shift
        shifts=[spearmanr(rk(anchored_entropy(M+d,anchor,T)), rk(val)).correlation for d in [-1,-0.5,0.5,1]]
        d3=min(shifts)
        # D4: add a sub-threshold (weak) off-target at 5.3
        base=anchored_entropy(M,anchor,T); new=anchored_entropy(np.hstack([M,np.full((n_drugs,1),5.3)]),anchor,T)
        d4=(new<=base+1e-9).mean()*100
        # p*
        ref=rk(val); pstar=None; np.random.seed(0)
        for ps in panel:
            v=[spearmanr(ref,rk(anchored_entropy(M[:,np.random.choice(n_kin,ps,replace=False)],anchor,T))).correlation for _ in range(30)]
            if np.mean(v)>0.90: pstar=ps; break
        print(f"{anchor:4s}  {T:.1f}   {orient:+.3f}                     {d3:+.3f}          {d4:5.0f}       {pstar}")
print("\nGoal: orient > +0.9 AND D3 ~ +1 (shift-invariant) AND D4 = 100. If found -> repaired candidate.")
