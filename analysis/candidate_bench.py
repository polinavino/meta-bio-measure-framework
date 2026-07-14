"""
Benchmark the repo's EXACT softplus-entropy candidate against the other measures and the
average-rank canonical extension. Investigate the suspected floor-domination inversion.
Functions copied verbatim from selectivity/candidate_measure.py.
"""
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata
M = pd.read_csv("/Users/polina/Documents/BioInfStuff/selectivity/klaeger_matrix.csv", index_col=0).values.astype(float)
n_drugs, n_kin = M.shape
FLOOR, TAU_STAR, T = 5.0, 6.0, 1.0
def to_ranks(s): return len(s) - s.argsort().argsort()   # repo's exact ranker (rank1 = largest value)

# ---- verbatim repo functions ----
def candidate(P, floor=FLOOR, T=T, eps=1e-10):
    w = T*np.logaddexp(0.0,(P-floor)/T)
    rs = np.where(w.sum(1,keepdims=True)==0, eps, w.sum(1,keepdims=True)); p=w/rs
    return -(-(p*np.where(p>0,np.log2(p+eps),0)).sum(1))
def entropy(P,b=5.0,e=1e-10):
    s=np.maximum(P-b,0); rs=np.where(s.sum(1,keepdims=True)==0,e,s.sum(1,keepdims=True)); p=s/rs
    return -(-(p*np.where(p>0,np.log2(p+e),0)).sum(1))
def gini(P,b=5.0):
    s=np.maximum(P-b,0); out=[]
    for r in s:
        rs=np.sort(r);n=len(rs);t=rs.sum(); out.append(0.0 if t==0 else (2*np.sum(np.arange(1,n+1)*rs))/(n*t)-(n+1)/n)
    return np.array(out)
def s_score(P,thr=6.0): return -(P>thr).astype(float).mean(1)
def ratio(P,k=1): return np.array([np.sort(r)[::-1][0]-max(np.sort(r)[::-1][k] if len(r)>k else 5.0,5.0) for r in P])

# all measures higher = more selective; rank via to_ranks -> rank1 = most selective
funcs={'candidate':candidate,'entropy':entropy,'gini':gini,'s_score':s_score,'ratio':ratio}
ranks={k:to_ranks(f(M)) for k,f in funcs.items()}

print("=== Spearman among measure RANKINGS on full Klaeger (rank1=most selective) ===")
keys=list(funcs)
print("           "+" ".join(f"{k:>9s}" for k in keys))
for a in keys:
    print(f"{a:10s} "+" ".join(f"{spearmanr(ranks[a],ranks[b]).correlation:+9.3f}" for b in keys))

print("\n=== Is the candidate floor-dominated? candidate vs entropy at different floors ===")
for fl in [5.0,4.0,3.0,2.0,0.0]:
    c=to_ranks(candidate(M,floor=fl))
    print(f"  floor={fl:.1f}: Spearman(candidate, hard-entropy) = {spearmanr(c,ranks['entropy']).correlation:+.3f}")
# diagnostic: for a maximally selective compound (1 spike) what does candidate say?
frac_at_floor=(M<=5.0).mean()
print(f"\n  fraction of Klaeger entries at/below floor 5.0 = {frac_at_floor:.3f}  (pedestal weight per kinase = T*log(2) = {T*np.log(2):.3f})")

# gated-in comparison + vs a quick average-rank proxy (mean of the other 4 measure ranks = consensus)
na=(M>6.0).sum(1); g=na>0
cons=np.mean([ranks[k][g] for k in ['entropy','gini','s_score','ratio']],axis=0)  # consensus of the 4 (excl candidate)
print(f"\n=== On gated-in (n={g.sum()}): candidate vs consensus-of-4 and vs each ===")
cg=ranks['candidate'][g]
print(f"  Spearman(candidate, consensus-of-4) = {spearmanr(cg, cons).correlation:+.3f}")
for k in ['entropy','gini','s_score','ratio']:
    print(f"  Spearman(candidate, {k:8s}) = {spearmanr(cg, ranks[k][g]).correlation:+.3f}")
