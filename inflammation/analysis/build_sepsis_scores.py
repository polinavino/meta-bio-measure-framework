"""
Topic 3 (sepsis transcriptomic signatures): score competing published signatures per sample on GSE65682.

Concept (latent) = inflammatory/immune dysregulation severity in sepsis. Competing MEASURES = per-sample
scores of published signatures, each a mean z-score over its genes (a "signature activity score"):
  Inflammatory   = HALLMARK_INFLAMMATORY_RESPONSE
  IFNg           = HALLMARK_INTERFERON_GAMMA_RESPONSE
  SRS7           = Davenport 7-gene Sepsis Response Signature set
  SRSq19         = SepstratifieR 19-gene set
  MARS8          = Scicluna MARS 4-endotype bi-signature genes
External anchor = 28-day mortality. Cohorts = MARS discovery vs validation (independent) for step 3.

NOTE (honest): mean-z over a gene SET is an approximation. It is the right object for co-directional
Hallmark sets; for SRS/SRSq/MARS (mixed-direction classifier genes) it is a proxy, not the official
PCA-based SRSq or the MARS classifier. We sign-align each score to the panel consensus in the analyses
and report the approximation. Deterministic.

Output: inflammation/data/sepsis_scores.csv (one row per sample).
"""
import json
import numpy as np, pandas as pd

DATA = "/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/inflammation/data"

expr = pd.read_parquet(f"{DATA}/sepsis_expr.parquet")               # probes x samples
# GEO platform table: skip the '#'/'!'-prefixed header block, read from the 'ID\t...' row
with open(f"{DATA}/GPL13667_probe_annotation.txt") as fh:
    hdr = next(i for i,l in enumerate(fh) if l.startswith("ID\t"))
ann = pd.read_csv(f"{DATA}/GPL13667_probe_annotation.txt", sep="\t", skiprows=hdr,
                  low_memory=False)
ann = ann[~ann["ID"].astype(str).str.startswith("!")]              # drop platform_table_end
sigs = json.load(open(f"{DATA}/signatures.json"))
pheno = pd.read_csv(f"{DATA}/sepsis_pheno.csv")

# probe -> symbol (old-platform symbols); collapse to symbol by mean across probes
idcol = [c for c in ann.columns if c.lower() in ("id","probe","probe_id")][0]
symcol = [c for c in ann.columns if "symbol" in c.lower()][0]
p2s = ann.set_index(idcol)[symcol].astype(str)
expr = expr.copy(); expr.index = expr.index.map(p2s)
expr = expr[expr.index.notna() & (expr.index != "nan")]
sym = expr.groupby(level=0).mean()                                  # genes x samples

# alias map: signature current symbols -> GPL13667 legacy symbols
ALIAS = {"ADGRE3":"EMR3", "ARL14EP":"C11orf46", "GLTSCR2":"GLTSCR2", "NOP53":"GLTSCR2"}
def resolve(g): return g if g in sym.index else ALIAS.get(g, g)

Z = sym.sub(sym.mean(1), axis=0).div(sym.std(1).replace(0, np.nan), axis=0)   # gene z across samples
def score(genes):
    present = [resolve(g) for g in genes if resolve(g) in Z.index]
    return Z.loc[present].mean(0), len(present), len(genes)

out = pd.DataFrame(index=sym.columns)
print("=== Sepsis signature scoring (GSE65682) ===")
print(f"expr {expr.shape[0]} probe-rows -> {sym.shape[0]} symbols x {sym.shape[1]} samples\n")
name_map = {"HALLMARK_INFLAMMATORY_RESPONSE":"Inflammatory","HALLMARK_INTERFERON_GAMMA_RESPONSE":"IFNg",
            "SRS_Davenport_7gene":"SRS7","SRSq_SepstratifieR_19gene":"SRSq19","MARS_endotype_bisignature":"MARS8"}
for key,short in name_map.items():
    if key not in sigs: print(f"  {short}: signature '{key}' not in signatures.json — skipped"); continue
    s,present,total = score(sigs[key]); out[short] = s
    print(f"  {short:12s} genes {present}/{total} present")

# attach phenotype (coerce string 'NA')
pheno = pheno.set_index("gsm")
out["mortality"] = pd.to_numeric(pheno["mortality_event_28days"].reindex(out.index).replace("NA",np.nan), errors="coerce")
out["endotype"]  = pheno["endotype_class"].reindex(out.index).replace("NA",np.nan)
out["cohort"]    = pheno["endotype_cohort"].reindex(out.index).replace("NA",np.nan)
out["age"]       = pd.to_numeric(pheno["age"].reindex(out.index).replace("NA",np.nan), errors="coerce")
out.index.name = "gsm"
out.to_csv(f"{DATA}/sepsis_scores.csv")

sc = list(name_map.values())
print(f"\nsamples: {len(out)}; with mortality: {out['mortality'].notna().sum()} "
      f"(deaths={int(out['mortality'].sum())}); with endotype: {out['endotype'].notna().sum()}")
print("cohorts:", out['cohort'].value_counts(dropna=False).to_dict())
print("\nSpearman among raw signature scores:")
print(out[sc].corr("spearman").round(2).to_string())
print(f"\nwrote sepsis_scores.csv")
