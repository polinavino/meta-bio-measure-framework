"""Descriptive cross-domain separation/reversal curves.

Signs are fixed by declared conventions. Kinase eligibility is applied, serotonin
uses observed receptors and a log gap, and clocks are restricted to the position
type. The outcome is the fraction of indicator pairs in opposite strict directions;
ties are not reversals. Consensus axes share the scores under comparison. Distances
between binned curves do not establish held-out predictive transfer or a bound on
concept-level disagreement.
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
    M = M[(M > 6.0).any(1)]
    return np.column_stack([-s_score(M,6.0), ent(M), gini(M), ratio(M)])   # orient roughly, sign-align fixes rest

def serotonin():
    from serotonin_controls import load, score_rows
    matrix, measured = load(both_targets=True)
    X = score_rows(matrix)  # observed-only; top-two log gap, explicit max>6 gate
    return X[np.isfinite(X).all(1)]

def clocks():
    df = pd.read_parquet("/Users/polina/Documents/BioInfStuff/epigenetic-clock-desiderata/data/clock_outputs.parquet")
    # Position clocks only: rate/deviation outputs are different declared types.
    return df[["Horvath","Hannum"]].to_numpy(float)

def smoking():
    a=pd.read_csv("/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/analysis/data/smoke50660.csv")
    b=pd.read_csv("/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/analysis/data/smoke42861.csv")
    X=pd.concat([a,b])[["ahrr","epismoke","joehanes","epitob"]].to_numpy(float)
    return -X[~np.isnan(X).any(1)]

print("Descriptive consensus-axis curves; no held-out prediction or mathematical bound on constitutive disagreement.")
print("Discordance = fraction of indicator pairs with strict reversals; signs fixed by convention.")
DOMAINS = {"kinase":kinase(),"serotonin":serotonin(),"clocks":clocks(),"smoking":smoking()}

def sign_align_ranks(X):
    # X is oriented by declared conventions in its loader; disagreement never sets signs.
    return np.column_stack([rankdata(X[:,j]) for j in range(X.shape[1])])

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
    # Fraction of indicator pairs in opposite strict directions; ties get no reversal credit.
    signs = np.sign(R[i] - R[j])
    positive, negative = (signs > 0).sum(1), (signs < 0).sum(1)
    disc = positive * negative / (k * (k - 1) / 2)
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

# cross-domain transfer error: use domain A's descriptive decile curve with domain B's; RMSE over deciles
names=list(curves)
print("\nCross-domain descriptive curve RMSE (rows/columns=domains; curve distances only):")
print("           "+" ".join(f"{n[:8]:>8s}" for n in names))
for a in names:
    row=[]
    for b in names:
        rmse=np.sqrt(np.nanmean((curves[a]-curves[b])**2))
        row.append(rmse)
    print(f"{a:10s} "+" ".join(f"{v:8.3f}" for v in row))
within=np.mean([np.sqrt(np.nanmean((curves[a]-curves[a])**2)) for a in names])  # 0 by def
offdiag=[np.sqrt(np.nanmean((curves[a]-curves[b])**2)) for a in names for b in names if a!=b]
print(f"\nMean off-diagonal descriptive cross-domain RMSE = {np.mean(offdiag):.3f} (discordance is on a 0-1 scale)")
print(f"Between-domain spread of discordance at each decile (std across domains):")
allc=np.vstack([curves[n] for n in names])
print("  "+" ".join(f"{s:.2f}" for s in np.nanstd(allc,0)))
# monotone-decreasing check per domain (Spearman of curve vs decile)
from scipy.stats import spearmanr
print("\nMonotonic decrease (Spearman curve vs decile; want strongly negative):")
for n in names:
    print(f"  {n:10s}: {spearmanr(np.arange(10), curves[n]).correlation:+.3f}")
