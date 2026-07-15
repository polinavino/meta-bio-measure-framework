"""
Build the unified TP53 variant x measure matrix from the acquired public sources.

Concept (latent) kappa = deleteriousness / functional impact of a TP53 missense variant.
Competing MEASURES (each claims to quantify kappa), left in their RAW native orientation here
(sign-alignment to "concept-positive = more deleterious" happens in the analysis scripts):

  8 Ishioka/Kato yeast transactivation promoters (%WT activity; HIGHER = more p53 function):
    WAF1, MDM2, BAX, h1433s(14-3-3s), AIP1(p53AIP1), GADD45, NOXA, P53R2      [NCI TP53 DB r21]
  7 computational predictors:
    AlphaMissense am_pathogenicity (0-1; HIGHER = worse)                       [Cheng 2023]
    EVE eve_score (0-1; HIGHER = worse), ESM1b LLR (more NEGATIVE = worse),
    REVEL, CADD phred, PrimateAI, BayesDel (all HIGHER = worse)                [EVE/ESM1b native; rest dbNSFP v4.1a]
  DMS (ProteinGym convention: HIGHER = higher proliferative "fitness"):
    Giacomelli 2018 {WT_Nutlin, Null_Nutlin, Null_Etoposide}, Kotler 2018

Field-PROPOSED consensus formula (for later benchmarking): median of the 8 promoters (Kato/Ishioka).

External ANCHOR (independent of the assay measures): ClinVar germline clinical significance + star rating.

All sources are on the canonical 393-aa P04637 / NM_000546 numbering. Join key = "R175H" style
(WT 1-letter)(codon)(ALT 1-letter). Deterministic; no randomness.

Output: tp53/data/tp53_matrix.csv  (one row per missense variant key; measure columns raw).
"""
import re
import numpy as np
import pandas as pd

# --- persist this script's console output to a tracked file (single source of truth for paper numbers) ---
import sys as _sys
_OUT = "/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/tp53/analysis/outputs/build_matrix.txt"
class _Tee:
    def __init__(self, p): self._f = open(p, "w"); self._o = _sys.stdout
    def write(self, s): self._o.write(s); self._f.write(s)
    def flush(self): self._o.flush(); self._f.flush()
_sys.stdout = _Tee(_OUT)

DATA = "/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/tp53/data"

AA3 = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G",
       "His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S",
       "Thr":"T","Trp":"W","Tyr":"Y","Val":"V","Ter":"*"}
AA1 = set("ACDEFGHIKLMNPQRSTVWY")

def is_missense_key(k):
    return isinstance(k,str) and len(k)>=3 and k[0] in AA1 and k[-1] in AA1 and k[0]!=k[-1] and k[1:-1].isdigit()

# ---------- 1. Ishioka 8-promoter transactivation (the field-native functional measures) ----------
ish = pd.read_csv(f"{DATA}/tp53_FunctionIshioka_r21.csv")
PROM = {"WAF1nWT":"waf1","MDM2nWT":"mdm2","BAXnWT":"bax","h1433snWT":"h1433s",
        "AIP1nWT":"aip1","GADD45nWT":"gadd45","NOXAnWT":"noxa","P53R2nWT":"p53r2"}
ish = ish.rename(columns=PROM)
ish["key"] = ish["AAchange"].astype(str).str.strip()
ish = ish[ish["key"].map(is_missense_key)].copy()
prom_cols = list(PROM.values())
# a few variants appear more than once -> average replicate rows
ish = ish.groupby("key", as_index=False)[prom_cols].mean()
ish["ish_median"] = ish[prom_cols].median(axis=1)   # field's proposed consensus formula (Kato/Ishioka)

# ---------- 2. AlphaMissense ----------
am = pd.read_csv(f"{DATA}/alphamissense_TP53_P04637.tsv", sep="\t")
am = am.rename(columns={"protein_variant":"key","am_pathogenicity":"am_path"})
am["key"] = am["key"].astype(str).str.strip()
am = am[am["key"].map(is_missense_key)][["key","am_path","am_class"]]

# ---------- 3. DMS assays (ProteinGym) ----------
def load_dms(fname, col):
    d = pd.read_csv(f"{DATA}/{fname}")[["mutant","DMS_score"]].rename(columns={"mutant":"key","DMS_score":col})
    d["key"] = d["key"].astype(str).str.strip()
    return d[d["key"].map(is_missense_key)]
giac_wt  = load_dms("P53_HUMAN_Giacomelli_2018_WT_Nutlin.csv","giac_wtnut")
giac_nn  = load_dms("P53_HUMAN_Giacomelli_2018_Null_Nutlin.csv","giac_nullnut")
giac_ne  = load_dms("P53_HUMAN_Giacomelli_2018_Null_Etoposide.csv","giac_nulletop")
kotler   = load_dms("P53_HUMAN_Kotler_2018.csv","kotler")

# ---------- 3b. additional computational predictors (all protein-keyed to the R175H key) ----------
def load_pred(fname, cols):
    d = pd.read_csv(f"{DATA}/{fname}"); d["key"] = d["key"].astype(str).str.strip()
    return d[["key"]+cols]
eve = load_pred("pred_eve_TP53.csv", ["eve_score"])            # 0-1, higher = more pathogenic
esm = load_pred("pred_esm1b_TP53.csv", ["esm1b_llr"])          # LLR, more NEGATIVE = more damaging
dbn = load_pred("pred_dbnsfp_TP53.csv", ["revel","cadd_phred","primateai","bayesdel_addaf"])  # all higher = worse

# ---------- 4. ClinVar external anchor (germline, GRCh38 only to dedupe) ----------
cv = pd.read_csv(f"{DATA}/clinvar_TP53_variant_summary.tsv", sep="\t", low_memory=False)
cv = cv[cv["Assembly"]=="GRCh38"].copy()
pat = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})")
def cv_key(name):
    m = pat.search(str(name))
    if not m: return None
    wt, pos, alt = AA3.get(m.group(1)), m.group(2), AA3.get(m.group(3))
    if wt is None or alt is None: return None
    k = f"{wt}{pos}{alt}"
    return k if is_missense_key(k) else None
cv["key"] = cv["Name"].map(cv_key)
cv = cv[cv["key"].notna()].copy()
STARS = {"practice guideline":4,"reviewed by expert panel":3,
         "criteria provided, multiple submitters, no conflicts":2,
         "criteria provided, single submitter":1,
         "criteria provided, conflicting classifications":1,
         "criteria provided, conflicting interpretations":1}
cv["stars"] = cv["ReviewStatus"].map(lambda s: STARS.get(str(s).strip(),0))
# collapse to one row/variant: keep the highest-star assertion
cv = cv.sort_values("stars", ascending=False).drop_duplicates("key")
def simple(sig):
    s=str(sig).lower()
    if "conflict" in s: return -1
    if "pathogenic" in s and "benign" not in s: return 1
    if "benign" in s and "pathogenic" not in s: return 0
    return -1   # uncertain / other
cv["clinvar_simple"] = cv["ClinicalSignificance"].map(simple)
cv = cv.rename(columns={"ClinicalSignificance":"clinvar_sig"})[
        ["key","clinvar_sig","clinvar_simple","stars"]].rename(columns={"stars":"clinvar_stars"})

# ---------- merge ----------
m = ish.merge(am, on="key", how="outer")
for d in (giac_wt, giac_nn, giac_ne, kotler, eve, esm, dbn):
    m = m.merge(d, on="key", how="outer")
m = m.merge(cv, on="key", how="left")
# derive wt/pos/alt
m["wt"]  = m["key"].str[0]
m["pos"] = m["key"].str.extract(r"(\d+)").astype(int)
m["alt"] = m["key"].str[-1]
m = m.sort_values("pos").reset_index(drop=True)

pred_cols = ["am_path","eve_score","esm1b_llr","revel","cadd_phred","primateai","bayesdel_addaf"]
front = ["key","wt","pos","alt"]
m = m[front + prom_cols + ["ish_median"] + pred_cols + ["am_class",
      "giac_wtnut","giac_nullnut","giac_nulletop","kotler",
      "clinvar_sig","clinvar_simple","clinvar_stars"]]
out = f"{DATA}/tp53_matrix.csv"
m.to_csv(out, index=False)

# ---------- report ----------
print(f"Wrote {out}: {len(m)} unique missense variant keys")
print("\nNon-null coverage per column:")
for c in prom_cols+["ish_median"]+pred_cols+["giac_wtnut","giac_nullnut","giac_nulletop","kotler",
                     "clinvar_simple"]:
    print(f"  {c:16s} {m[c].notna().sum():5d}")
print("\nComplete-panel counts (intersection):")
ish_ok = m[prom_cols].notna().all(axis=1)
print(f"  all 8 promoters present:            {ish_ok.sum()}")
print(f"  all 8 promoters + AlphaMissense:    {(ish_ok & m['am_path'].notna()).sum()}")
core = ish_ok & m["am_path"].notna()
for c in ["giac_wtnut","kotler"]+[p for p in pred_cols if p!="am_path"]:
    print(f"  core + {c:16s} {(core & m[c].notna()).sum()}")
print(f"  core + all 7 predictors:            {(core & m[[p for p in pred_cols]].notna().all(axis=1)).sum()}")
print("\nClinVar labels among core panel (all 8 promoters + AM):")
print("  ", m.loc[core,"clinvar_simple"].value_counts(dropna=False).to_dict(),
      "  (1=path,0=benign,-1=VUS/conflict,NaN=absent)")
print("  path/benign with >=1 star:",
      ((m["clinvar_simple"].isin([0,1])) & (m["clinvar_stars"]>=1) & core).sum())
