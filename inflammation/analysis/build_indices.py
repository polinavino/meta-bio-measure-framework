"""
Topic 1 (clinical inflammatory indices): compute the competing indices from raw NHANES blood counts.

Concept (latent) = systemic inflammatory burden. Competing MEASURES (all concept-positive: higher = more
inflammation), from routine CBC + CRP + albumin:
  NLR  = neutrophils / lymphocytes
  PLR  = platelets  / lymphocytes
  MLR  = monocytes  / lymphocytes
  SII  = platelets * neutrophils / lymphocytes      (systemic immune-inflammation index)
  SIRI = neutrophils * monocytes / lymphocytes      (systemic inflammation response index)
  CRP  = high-sensitivity C-reactive protein (mg/L)
  CAR  = CRP / albumin                              (CRP-to-albumin ratio)

Two independent NHANES cycles (2015-16, 2017-18) are kept separate = two cohorts for the reproducibility
test. External anchor = all-cause mortality (linked mortality file) + chronological age. Deterministic.

Output: inflammation/data/indices_2015.csv, indices_2017.csv (one row per adult with a usable CBC).
"""
import numpy as np, pandas as pd

DATA = "/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/inflammation/data"
IDX = ["NLR","PLR","MLR","SII","SIRI","CRP","CAR"]

def build(cycle):
    d = pd.read_csv(f"{DATA}/nhanes_{cycle}.csv")
    d = d[(d["lymph"] > 0)].copy()                      # ratios need lymphocytes > 0
    d["NLR"]  = d["neut"]     / d["lymph"]
    d["PLR"]  = d["platelet"] / d["lymph"]
    d["MLR"]  = d["mono"]     / d["lymph"]
    d["SII"]  = d["platelet"] * d["neut"] / d["lymph"]
    d["SIRI"] = d["neut"]     * d["mono"] / d["lymph"]
    d["CRP"]  = d["crp"]
    d["CAR"]  = d["crp"] / d["albumin"].where(d["albumin"] > 0)
    out = d[["seqn","age","sex","mortstat","permth"] + IDX]
    out.to_csv(f"{DATA}/indices_{cycle}.csv", index=False)
    return out

print("=== NHANES inflammatory indices ===")
for cyc in ["2015","2017"]:
    o = build(cyc)
    cc_idx  = o[IDX].notna().all(axis=1)
    cc_mort = cc_idx & o["mortstat"].notna()
    print(f"\ncycle {cyc}: {len(o)} rows with a CBC")
    print("  non-null per index:", {c:int(o[c].notna().sum()) for c in IDX})
    print(f"  complete 7-index panel: {cc_idx.sum()}")
    print(f"  complete panel + mortality: {cc_mort.sum()}"
          f"  (deaths={int(o.loc[cc_mort,'mortstat'].sum())})")
    print(f"  wrote indices_{cyc}.csv")
