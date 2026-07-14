"""
Cross-domain demonstration of two protocol refinements:

(1) ORIENTATION CHECK (G2). Each measure, oriented to "more of the concept", should track a crude,
    theory-light EXTERNAL ANCHOR in the expected direction. This is the check that catches a measure
    which is internally consistent (stable/monotone) yet points backwards. Control: the un-gated kinase
    candidate, which FAILS it.
(2) EXTERNAL ANCHOR VALIDATES THE CONSENSUS (formal-spine §1.1). The consensus order of the measures
    should agree with the anchor especially at the EXTREMES (where both are unambiguous). We report the
    overall Spearman and the Spearman restricted to anchor-extreme vs anchor-middle objects.

Anchors:  kinase / serotonin selectivity -> n_active (fewer active targets = more selective; NEG expected)
          clocks biological age           -> chronological age                     (POS expected)
          smoking exposure                -> never<former<current ordinal          (POS expected)

Concept-positive orientation of each measure uses its KNOWN convention (not the anchor -> non-circular):
  selectivity: -s_score, -entropy, +gini, +ratio ;  clocks: +clock ;  smoking exposure: -signature.
Deterministic.
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
    df=pd.read_csv("/Users/polina/Documents/BioInfStuff/psychedelic-selectivity/data/serotonin_receptor_matrix.csv")
    R=df.iloc[:,1:].to_numpy(float); keep=(~np.isnan(R)).sum(1)>=4      # >=4 receptors measured -> N=949 (broader than §4.3's 297-compound 2A/2B-overlap set)
    R=R[keep]; Rf=np.where(np.isnan(R),5.0,R)                            # missing = non-binding floor (pKi 5)
    meas,na=selectivity_measures(Rf); return meas, na, "n_active(pKi>6)", -1
def load_clocks():
    co=pd.read_parquet("/Users/polina/Documents/BioInfStuff/epigenetic-clock-desiderata/data/clock_outputs.parquet")
    ins=pd.read_parquet("/Users/polina/Documents/BioInfStuff/epigenetic-clock-desiderata/data/instability_scores.parquet")
    age=ins["age"].reindex(co.index).to_numpy(float); ok=~np.isnan(age)
    meas={c: co[c].to_numpy(float)[ok] for c in co.columns}             # clocks: higher = older (concept-positive)
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
print("(2) EXTERNAL ANCHOR VALIDATES THE CONSENSUS — agreement, esp. at the extremes")
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
print("\n(Expect: strong overall agreement; monotone consensus across anchor tertiles;")
print(" and higher agreement at the extremes than in the middle — the concept decides the extremes.)")
