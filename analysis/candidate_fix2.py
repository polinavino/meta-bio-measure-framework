"""
Confirm the tension: does the ORIGINAL candidate pass D3 only because pedestal-domination makes it
baseline-insensitive (and inverted)? Does ANY smoothing width give BOTH correct orientation AND D3?
Also test a scale-free variant (per-compound anchor).
"""
import numpy as np, pandas as pd
from scipy.stats import spearmanr
M = pd.read_csv("/Users/polina/Documents/BioInfStuff/selectivity/klaeger_matrix.csv", index_col=0).values.astype(float)
def to_ranks(s): return len(s)-s.argsort().argsort()
def entropy(P,b=5.0,e=1e-10):
    s=np.maximum(P-b,0); rs=np.where(s.sum(1,keepdims=True)==0,e,s.sum(1,keepdims=True)); p=s/rs
    return -(-(p*np.where(p>0,np.log2(p+e),0)).sum(1))
ENT_REF=to_ranks(entropy(M))
betas=np.arange(4.5,6.6,0.25)
def orient(fn,**k):    # + = correct (agrees w/ selectivity direction of entropy)
    return spearmanr(to_ranks(fn(M,**k)),ENT_REF).correlation
def d3(fn):            # worst-case pairwise Spearman as floor varies -> want ~+1
    rr=np.array([to_ranks(fn(M,floor=b)) for b in betas])
    return min(spearmanr(rr[i],rr[j]).correlation for i in range(len(betas)) for j in range(i+1,len(betas)))

def c_orig(P,floor=5.0,T=1.0,e=1e-10):
    w=T*np.logaddexp(0.0,(P-floor)/T); rs=np.where(w.sum(1,keepdims=True)==0,e,w.sum(1,keepdims=True)); p=w/rs
    return -(-(p*np.where(p>0,np.log2(p+e),0)).sum(1))
def c_fixed(P,floor=5.0,T=1.0,e=1e-10):
    w=np.maximum(T*(np.logaddexp(0.0,(P-floor)/T)-np.log(2.0)),0.0)
    rs=np.where(w.sum(1,keepdims=True)==0,e,w.sum(1,keepdims=True)); p=w/rs
    return -(-(p*np.where(p>0,np.log2(p+e),0)).sum(1))

print("T-sweep: orientation (vs entropy, + = correct) and D3 worst-case (want ~+1)")
print(f"{'T':>4} | {'orig orient':>11} {'orig D3':>8} | {'fixed orient':>12} {'fixed D3':>9}")
for T in [0.5,1.0,2.0,3.0,5.0]:
    print(f"{T:4.1f} | {orient(c_orig,T=T):>11.3f} {d3(lambda P,floor:c_orig(P,floor,T)):>8.3f} | "
          f"{orient(c_fixed,T=T):>12.3f} {d3(lambda P,floor:c_fixed(P,floor,T)):>9.3f}")

# scale-free variant: anchor per-compound (floor = per-row min), pedestal-subtracted, then entropy.
# baseline shift affects all kinases equally -> relative structure invariant -> should pass D3 AND orient.
def c_scalefree(P,floor=5.0,T=1.0,e=1e-10):     # 'floor' arg ignored except as global shift test
    rowmin=P.min(1,keepdims=True)
    w=np.maximum(T*(np.logaddexp(0.0,(P-rowmin)/T)-np.log(2.0)),0.0)
    rs=np.where(w.sum(1,keepdims=True)==0,e,w.sum(1,keepdims=True)); p=w/rs
    return -(-(p*np.where(p>0,np.log2(p+e),0)).sum(1))
def d3_shift(fn):   # D3 as GLOBAL additive shift of all values (true baseline invariance)
    rr=np.array([to_ranks(fn(M+delta)) for delta in np.arange(-1.0,1.01,0.25)])
    return min(spearmanr(rr[i],rr[j]).correlation for i in range(len(rr)) for j in range(i+1,len(rr)))
print("\nScale-free (per-compound-min anchor) variant:")
print(f"  orientation vs entropy: {spearmanr(to_ranks(c_scalefree(M)),ENT_REF).correlation:+.3f}")
print(f"  D3 under global shift [-1,+1]: worst pairwise = {d3_shift(c_scalefree):+.3f}")
print(f"  (compare fixed-candidate under global shift: {d3_shift(lambda P:c_fixed(P)):+.3f}; "
      f"orig: {d3_shift(lambda P:c_orig(P)):+.3f})")
