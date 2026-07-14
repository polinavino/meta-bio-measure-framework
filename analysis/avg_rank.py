"""
Average-rank-over-linear-extensions canonical measure (Patil-Taillie; Bruggemann/De Loof),
demonstrated on the kinase (Klaeger) selectivity domain.

Consensus poset M-hat: a >= b iff all 4 canonical measures rank a at least as selective as b.
Canonical measure = mean rank of each compound over uniformly-sampled linear extensions of
that poset (weight-free, commits no further than the consensus order forces).

Sampler: Bubley-Dyer lazy adjacent-transposition Markov chain (uniform over linear extensions).
Convergence checked by two independent chains. Then benchmark the 4 existing measures and the
repo's reliability-gated softplus-entropy candidate against the canonical average-rank.
"""
import numpy as np, pandas as pd
from scipy.stats import rankdata, spearmanr

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
def softplus_entropy(P, floor=5.0, tau_star=6.0, T=1.0):
    # repo candidate: reliability gate + softplus hinge; higher = more selective (neg entropy)
    w = T*np.logaddexp(0.0,(P-floor)/T); rs=w.sum(1,keepdims=True); rs=np.where(rs==0,1e-10,rs)
    p=w/rs; H=-(p*np.where(p>0,np.log2(p+1e-12),0)).sum(1)
    val = -H                                  # higher = more selective
    val[P.max(1)<=tau_star] = np.nan          # D1 gate
    return val
def to_rank(s): return len(s)-rankdata(s,method="ordinal")+1  # 1 = most selective (highest score)

s_thr=np.arange(5.5,8.25,.25); base=np.arange(5.0,6.75,.25)
# gate to active compounds
na=(M>6.0).sum(1); g=na>0; Mg=M[g]; N=g.sum()
# 4 canonical measures -> median rank (1=most selective); higher selectivity score => rank 1
r_s = np.median([to_rank(-s_score(Mg,t)) for t in s_thr],0)
r_e = np.median([to_rank(-ent(Mg,b)) for b in base],0)
r_g = np.median([to_rank(gini(Mg,b)) for b in base],0)
r_r = np.median([to_rank(ratio(Mg,k)) for k in range(1,6)],0)
R = np.vstack([r_s,r_e,r_g,r_r]).T          # (N,4), lower rank = more selective

# consensus poset: dom[a,b]=True iff a at least as selective as b under ALL 4 (rank_a<=rank_b)
dom = np.all(R[:,None,:] <= R[None,:,:], axis=2)     # (N,N)
comparable = dom | dom.T
np.fill_diagonal(comparable, True)
frac_incomp = 1 - (comparable.sum()-N)/(N*(N-1))
print(f"Gated-in compounds N={N}; incomparable-pair fraction in consensus poset = {frac_incomp:.3f}")

def sample_avg_rank(seed, burn=200_000, nsamp=2000, gap=1000):
    rng=np.random.default_rng(seed)
    order=list(np.argsort(R.mean(1)))          # mean-rank sort = valid linear extension
    pos_of=np.empty(N,int);
    for p,x in enumerate(order): pos_of[x]=p
    rank_sum=np.zeros(N); ns=0; step=0; total=burn+nsamp*gap
    for step in range(total):
        i=rng.integers(N-1); a=order[i]; b=order[i+1]
        if not comparable[a,b]:                # incomparable -> lazy swap
            if rng.random()<0.5:
                order[i],order[i+1]=b,a
        if step>=burn and (step-burn)%gap==0:
            for p,x in enumerate(order): rank_sum[x]+=(p+1)   # 1=most selective
            ns+=1
    return rank_sum/ns

print("Running two independent Bubley-Dyer chains ...")
ar1=sample_avg_rank(1); ar2=sample_avg_rank(2)
conv=spearmanr(ar1,ar2).correlation
print(f"Chain convergence: Spearman(chain1, chain2) = {conv:.4f}  (want ~1.0)")
avg_rank=(ar1+ar2)/2                          # canonical measure (1=most selective)

# benchmark: Spearman of each measure's ranking vs canonical average-rank
print("\nSpearman of each existing measure's ranking vs canonical average-rank extension:")
for name,rk in [("S-score",r_s),("entropy",r_e),("Gini",r_g),("ratio",r_r)]:
    print(f"  {name:9s}: {spearmanr(rk, avg_rank).correlation:+.3f}")
# NOTE: candidate benchmarking is done authoritatively in candidate_bench.py / candidate_fix.py
# using the repo's EXACT candidate_measure.py. The line below used a quick reimplementation and is
# omitted here to avoid a second, non-authoritative number. See those scripts for the candidate result.

# how far is each measure from canonical, in mean |rank difference|?
print("\nMean |rank - canonical| (rank units, N compounds):")
for name,rk in [("S-score",r_s),("entropy",r_e),("Gini",r_g),("ratio",r_r)]:
    print(f"  {name:9s}: {np.abs(rk-avg_rank).mean():5.1f}")
print("\nInterpretation: the measure closest to the canonical average-rank extension is the")
print("best single proxy for the weight-free consensus aggregate; large gaps = commits beyond consensus.")
