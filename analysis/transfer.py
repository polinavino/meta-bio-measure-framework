"""
Cross-domain transfer: does the (normalized separation -> inter-measure discordance) relationship
transfer across domains? Calibrate the curve in one domain, predict others out-of-domain.

Domains (different concepts, measures, object types):
  kinase   selectivity   4 measures (S-score,entropy,Gini,ratio)   compounds  (Klaeger 222)
  serotonin selectivity  4 measures (precomputed)                  compounds  (~13.5k)
  clocks   biological age 5 clocks (Horvath..DunedinPACE)          samples    (1385)
  smoking  exposure       4 signatures (AHRR,EpiSmoke,Joehanes,EpiTob) samples (~1153)

Per domain: sign-align measures to consensus direction; sample pairs; for each pair compute
  separation s = |percentile_i - percentile_j|  in [0,1]  (percentile on consensus rank)
  discordance d = min(a,k-a)/floor(k/2)  in [0,1]  (a = #measures ranking i above j; 0=unanimous,1=max split)
Bin s into deciles -> mean d per bin. Overlay; then cross-domain RMSE (fit A, predict B) vs within-domain.
"""
import numpy as np, pandas as pd
from scipy.stats import rankdata
rng = np.random.default_rng(0)

# ---------- load per-domain measure matrices (objects x measures) ----------
def kinase():
    M = pd.read_csv("/Users/polina/Documents/BioInfStuff/selectivity/klaeger_matrix.csv", index_col=0).values.astype(float)
    def s_score(P,t): return (P>t).mean(1)
    def ent(P,b=5.0,e=1e-10):
        sh=np.maximum(P-b,0); rs=np.where(sh.sum(1,keepdims=True)==0,e,sh.sum(1,keepdims=True)); p=sh/rs
        return -(-(p*np.where(p>0,np.log2(p+e),0)).sum(1))
    def gini(P,b=5.0):
        sh=np.maximum(P-b,0); o=[]
        for r in sh:
            rs=np.sort(r);n=len(rs);t=rs.sum(); o.append(0.0 if t==0 else (2*np.sum(np.arange(1,n+1)*rs))/(n*t)-(n+1)/n)
        return np.array(o)
    def ratio(P,k=1): return np.array([np.sort(r)[::-1][0]-max(np.sort(r)[::-1][k] if len(r)>k else 5.0,5.0) for r in P])
    return np.column_stack([-s_score(M,6.0), ent(M), gini(M), ratio(M)])   # orient roughly, sign-align fixes rest

def serotonin():
    df = pd.read_csv("/Users/polina/Documents/BioInfStuff/psychedelic-selectivity/data/selectivity_framework_results.csv")
    X = df[["s_score","entropy","gini","ratio"]].to_numpy(float)
    return X[~np.isnan(X).any(1)]

def clocks():
    df = pd.read_parquet("/Users/polina/Documents/BioInfStuff/epigenetic-clock-desiderata/data/clock_outputs.parquet")
    return df[["Horvath","Hannum","PhenoAge","GrimAge","DunedinPACE"]].to_numpy(float)

def smoking():
    a=pd.read_csv("/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/analysis/data/smoke50660.csv")
    b=pd.read_csv("/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/analysis/data/smoke42861.csv")
    X=pd.concat([a,b])[["ahrr","epismoke","joehanes","epitob"]].to_numpy(float)
    return X[~np.isnan(X).any(1)]

DOMAINS = {"kinase":kinase(),"serotonin":serotonin(),"clocks":clocks(),"smoking":smoking()}

def sign_align_ranks(X):
    # rank each measure; flip any measure negatively correlated with the mean rank so all point one way
    R = np.column_stack([rankdata(X[:,j]) for j in range(X.shape[1])])
    for _ in range(3):
        cons = R.mean(1)
        for j in range(R.shape[1]):
            if np.corrcoef(R[:,j], cons)[0,1] < 0:
                R[:,j] = (R.shape[0]+1) - R[:,j]
    return R

def pairs_curve(X, npairs=50000):
    N,k = X.shape
    R = sign_align_ranks(X)
    cons = R.mean(1)
    pct = (rankdata(cons)-1)/(N-1)
    # sample pairs (or all if small)
    allpairs = N*(N-1)//2
    if allpairs <= npairs:
        i,j = np.triu_indices(N,1)
    else:
        i = rng.integers(0,N,npairs); j = rng.integers(0,N,npairs)
        m=i!=j; i,j=i[m],j[m]
    a = (R[i] > R[j]).sum(1)                       # #measures ranking i above j
    disc = np.minimum(a, k-a) / (k//2)             # in [0,1]
    sep = np.abs(pct[i]-pct[j])                     # in [0,1]
    # decile-bin
    bins = np.clip((sep*10).astype(int),0,9)
    curve = np.array([disc[bins==b].mean() if (bins==b).any() else np.nan for b in range(10)])
    return curve, sep, disc, k

curves={}; raw={}
print(f"{'domain':10s} k  N     discordance by separation decile (0=near-tie .. 9=far)")
for name,X in DOMAINS.items():
    c,sep,disc,k = pairs_curve(X)
    curves[name]=c; raw[name]=(sep,disc)
    print(f"{name:10s} {k}  {X.shape[0]:5d} "+" ".join(f"{v:.2f}" for v in c))

# cross-domain transfer error: use domain A's decile curve to predict domain B's; RMSE over deciles
names=list(curves)
print("\nCross-domain prediction RMSE (rows=calibrate on, cols=predict) — low = transfers:")
print("           "+" ".join(f"{n[:8]:>8s}" for n in names))
for a in names:
    row=[]
    for b in names:
        rmse=np.sqrt(np.nanmean((curves[a]-curves[b])**2))
        row.append(rmse)
    print(f"{a:10s} "+" ".join(f"{v:8.3f}" for v in row))
within=np.mean([np.sqrt(np.nanmean((curves[a]-curves[a])**2)) for a in names])  # 0 by def
offdiag=[np.sqrt(np.nanmean((curves[a]-curves[b])**2)) for a in names for b in names if a!=b]
print(f"\nMean off-diagonal cross-domain RMSE = {np.mean(offdiag):.3f} (discordance is on a 0-1 scale)")
print(f"Between-domain spread of discordance at each decile (std across domains):")
allc=np.vstack([curves[n] for n in names])
print("  "+" ".join(f"{s:.2f}" for s in np.nanstd(allc,0)))
# monotone-decreasing check per domain (Spearman of curve vs decile)
from scipy.stats import spearmanr
print("\nMonotonic decrease (Spearman curve vs decile; want strongly negative):")
for n in names:
    print(f"  {n:10s}: {spearmanr(np.arange(10), curves[n]).correlation:+.3f}")
