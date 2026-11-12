"""Dependence-aware orientation checks and mean-rank anchor associations.

Directions are fixed by convention. Active count shares selectivity inputs, age
trains position clocks, and smoking labels are metadata. These checks do not
independently validate an entire consensus order. Clock summaries pool only the
position type; serotonin leaves unmeasured receptors unknown. The historical
ungated softplus is an explicitly obsolete negative control.
"""
import numpy as np, pandas as pd
from scipy.stats import rankdata, spearmanr

def sr(a, b): return spearmanr(a, b).correlation

# ---- measure functions (kinase conventions; raw: s_score,entropy higher=LESS selective) ----
def s_score(P,t=6.0): return (P>t).mean(1)
def entropy(P,b=5.0,e=1e-10):
    s=np.maximum(P-b,0); rs=np.where(s.sum(1,keepdims=True)==0,e,s.sum(1,keepdims=True)); p=s/rs
    return -(-(p*np.where(p>0,np.log2(p+e),0)).sum(1))          # returns -H (neg Shannon entropy): higher = MORE selective
def gini(P,b=5.0):
    s=np.maximum(P-b,0);o=[]
    for r in s:
        rs=np.sort(r);n=len(rs);t=rs.sum();o.append(0.0 if t==0 else (2*np.sum(np.arange(1,n+1)*rs))/(n*t)-(n+1)/n)
    return np.array(o)
def ratio(P,k=1): return np.array([np.sort(r)[::-1][0]-max(np.sort(r)[::-1][k] if len(r)>k else 5.0,5.0) for r in P])

def selectivity_measures(M, floor=5.0, active=6.0):
    # returned oriented CONCEPT-POSITIVE (higher = more selective).
    # NB: entropy() already returns -H (higher = more selective); s_score is fraction-hit (flip).
    return {"S-score": -s_score(M,active), "entropy": entropy(M,floor),
            "Gini": gini(M,floor), "ratio": ratio(M)}, (M>active).sum(1)

# ---- domain loaders: return (measures dict [concept-positive], anchor, anchor_name, expected_sign) ----
def load_kinase():
    M=pd.read_csv("/Users/polina/Documents/BioInfStuff/selectivity/klaeger_matrix.csv",index_col=0).values.astype(float)
    meas,na=selectivity_measures(M); return meas, na, "n_active", -1
def load_serotonin():
    from serotonin_controls import observed_anchor
    x, na, measured = observed_anchor()
    meas = dict(zip(["S-score", "entropy", "Gini", "ratio"], x.T))
    return meas, na, "observed n_active (shared-input proxy; variable coverage)", -1
def load_clocks():
    co=pd.read_parquet("/Users/polina/Documents/BioInfStuff/epigenetic-clock-desiderata/data/clock_outputs.parquet")
    ins=pd.read_parquet("/Users/polina/Documents/BioInfStuff/epigenetic-clock-desiderata/data/instability_scores.parquet")
    age=ins["age"].reindex(co.index).to_numpy(float); ok=~np.isnan(age)
    # A scalar consensus is formed only within the position-clock type.
    meas={c: co[c].to_numpy(float)[ok] for c in ["Horvath", "Hannum"]}
    return meas, age[ok], "chronological age", +1
def load_smoking():
    a=pd.read_csv("/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/analysis/data/smoke50660.csv")
    b=pd.read_csv("/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/analysis/data/smoke42861.csv")
    df=pd.concat([a,b],ignore_index=True)
    order={"0":0,"1":1,"2":2,"never":0,"former":1,"ex":1,"current":2}
    exp=df["smoking"].astype(str).map(order); ok=exp.notna()
    sigs=["ahrr","epismoke","joehanes","epitob"]
    meas={s: -df.loc[ok,s].to_numpy(float) for s in sigs}              # signatures lower=more exposure -> flip
    return meas, exp[ok].to_numpy(float), "never<former<current", +1

DOMAINS={"kinase":load_kinase,"serotonin":load_serotonin,"clocks":load_clocks,"smoking":load_smoking}

print("Anchors: active count shares inputs with selectivity scores; age is a clock training target.")
print("Smoking labels are external metadata. Correlation checks orientation, not causal validity.")
print("="*78)
print("(1) ORIENTATION CHECK — each concept-positive measure vs the external anchor")
print("    expected sign: kinase/serotonin NEG (more targets=less selective); clocks/smoking POS")
print("="*78)
consensus_store={}
for name,load in DOMAINS.items():
    meas,anchor,aname,exp=load(); N=len(anchor)
    print(f"\n{name}  (N={N}, anchor = {aname}; expected sign {'+' if exp>0 else '−'})")
    signs_ok=True
    for mname,mv in meas.items():
        r=sr(mv,anchor); ok = (r*exp)>0
        signs_ok &= ok
        print(f"   {mname:9s} vs anchor: {r:+.3f}   {'OK' if ok else 'MIS-ORIENTED'}")
    # consensus = mean rank of concept-positive measures (sign-aligned by known convention, not the anchor)
    Rk=np.column_stack([rankdata(mv) for mv in meas.values()]); cons=Rk.mean(1)
    consensus_store[name]=(cons,anchor,exp,aname)
    print(f"   → all measures correctly oriented: {signs_ok}")

# control: the un-gated kinase candidate (should FAIL orientation)
Mk=pd.read_csv("/Users/polina/Documents/BioInfStuff/selectivity/klaeger_matrix.csv",index_col=0).values.astype(float)
def cand_ungated(P,floor=5.0,T=1.0,e=1e-10):
    w=T*np.logaddexp(0.0,(P-floor)/T); rs=np.where(w.sum(1,keepdims=True)==0,e,w.sum(1,keepdims=True));p=w/rs
    return -(-(p*np.where(p>0,np.log2(p+e),0)).sum(1))                  # -H, intended higher=more selective
na_k=(Mk>6.0).sum(1)
print(f"\nCONTROL — un-gated kinase candidate (intended 'more selective', expected NEG vs n_active):")
print(f"   vs n_active: {sr(cand_ungated(Mk),na_k):+.3f}   -> POSITIVE = MIS-ORIENTED (the check catches it)")

print("\n"+"="*78)
print("(2) ANCHOR ASSOCIATIONS OF A MEAN-RANK SUMMARY — dependence must be declared")
print("="*78)
def fmt(x): return "n/a (constant)" if x is None or np.isnan(x) else f"{x:+.3f}"
for name,(cons,anchor,exp,aname) in consensus_store.items():
    c=cons*exp                                     # higher c ~ higher anchor expected
    cp=(rankdata(c)-1)/(len(c)-1)
    ar=rankdata(anchor)                            # anchor ranks (ties averaged) — robust to discrete anchors
    overall=sr(cp,anchor)
    # rank-tertiles (handle ties / discrete anchors like the 3-level smoking ordinal)
    q1,q2=np.quantile(ar,[1/3,2/3]); lo=ar<=q1; hi=ar>=q2; mid=(~lo)&(~hi)
    ext=lo|hi
    s_ext=sr(cp[ext],anchor[ext]); s_mid=sr(cp[mid],anchor[mid])
    pm=lambda mask: f"{cp[mask].mean():.2f}" if mask.sum() else "n/a"
    print(f"\n{name}: Spearman(consensus, anchor) overall = {overall:+.3f}")
    print(f"   mean consensus-pct by anchor tertile:  bottom {pm(lo)} | middle {pm(mid)} | top {pm(hi)}  (monotone = tracks)")
    print(f"   agreement Spearman:  extremes (n={ext.sum()}) {fmt(s_ext)} | middle (n={mid.sum()}) {fmt(s_mid)}")
    if mid.sum()==0: print("   (anchor too discrete/tied to form a middle tertile — see per-tertile bottom/top)")
    # Pooling bottom+top spans the whole anchor range while the middle spans one third,
    # and Spearman is attenuated by range restriction, so the pooled contrast above is
    # partly an artifact of the split. Each tertile alone carries the same restriction.
    _per = []
    for _nm, _m in [("bottom", lo), ("middle", mid), ("top", hi)]:
        _r = sr(cp[_m], anchor[_m]) if _m.sum() > 2 else float("nan")
        _per.append(f"{_nm} {fmt(_r)} (n={_m.sum()})")
    print(f"   per-tertile (comparable range restriction):  " + " | ".join(_per))
print("\n(Expect: strong overall agreement and a monotone consensus across anchor tertiles.")
print(" The pooled extremes-vs-middle gap is read against the per-tertile line: where one")
print(" extreme agrees no better than the middle, the gap is a floor or ceiling effect at the")
print(" other end rather than a symmetric extremes-versus-middle pattern.)")
