"""
Assemble the NHANES cardiovascular risk-factor table for both cycles.

Joins, on SEQN, within each cycle:
  DEMO   age, sex, race/ethnicity
  BIOPRO total cholesterol (LBXSCH, mg/dL)      [reused from ../../inflammation/data]
  HDL    HDL cholesterol (LBDHDD, mg/dL)
  BPX    systolic BP, mean of the available readings
  BPQ    on antihypertensive medication (BPQ050A)
  SMQ    current smoker (SMQ040 in {1,2})
  DIQ    diagnosed diabetes (DIQ010 == 1)
  BMX    BMI (BMXBMI), for the office-based score that uses no laboratory value
  MORT   NHANES Linked Mortality File, fixed-width, incl. UCOD_LEADING

The anchor is cause-specific CARDIOVASCULAR death: UCOD_LEADING 001 (diseases of
heart) or 005 (cerebrovascular disease).  This is deliberately crude and is not
the outcome any of the scores was fitted to, which is what makes it usable as an
external check rather than a restatement of the measures.

Writes ../data/cvd_2015.csv and ../data/cvd_2017.csv.
"""
import numpy as np
import pandas as pd
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
INFL = HERE.parent.parent / "inflammation" / "data"

LEAD_CVD = {"001", "005"}  # diseases of heart, cerebrovascular disease


def read_mort(path):
    """NHANES Linked Mortality File, public-use fixed-width layout."""
    rows = []
    for line in open(path):
        line = line.rstrip("\n")
        if not line.strip():
            continue
        rows.append({
            "seqn": int(line[0:14].strip()),
            "eligstat": line[14:15].strip(),
            "mortstat": line[15:16].strip(),
            "ucod": line[16:19].strip(),
            "permth": pd.to_numeric(line[43:46].strip(), errors="coerce"),
        })
    d = pd.DataFrame(rows)
    d["mortstat"] = pd.to_numeric(d["mortstat"], errors="coerce")
    d["cvd_death"] = ((d["mortstat"] == 1) & (d["ucod"].isin(LEAD_CVD))).astype(float)
    d.loc[d["eligstat"] != "1", ["mortstat", "cvd_death"]] = np.nan
    return d[["seqn", "mortstat", "cvd_death", "permth"]]


def build(sfx, mort_file, out_name):
    def xpt(folder, name):
        return pd.read_sas(folder / f"{name}_{sfx}.XPT")

    demo = xpt(INFL, "DEMO")[["SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH3"]]
    bio = xpt(INFL, "BIOPRO")[["SEQN", "LBXSCH"]]
    hdl = xpt(DATA, "HDL")[["SEQN", "LBDHDD"]]
    bmx = xpt(DATA, "BMX")[["SEQN", "BMXBMI"]]

    bpx = xpt(DATA, "BPX")
    sbp_cols = [c for c in ["BPXSY1", "BPXSY2", "BPXSY3", "BPXSY4"] if c in bpx.columns]
    bpx = bpx[["SEQN"] + sbp_cols].copy()
    bpx["sbp"] = bpx[sbp_cols].mean(axis=1, skipna=True)
    bpx = bpx[["SEQN", "sbp"]]

    bpq = xpt(DATA, "BPQ")[["SEQN", "BPQ050A"]].copy()
    bpq["bp_treated"] = (bpq["BPQ050A"] == 1).astype(float)

    smq = xpt(DATA, "SMQ")[["SEQN", "SMQ040"]].copy()
    smq["smoker"] = smq["SMQ040"].isin([1, 2]).astype(float)

    diq = xpt(DATA, "DIQ")[["SEQN", "DIQ010"]].copy()
    diq["diabetes"] = (diq["DIQ010"] == 1).astype(float)

    d = demo
    for part in [bio, hdl, bmx, bpx, bpq[["SEQN", "bp_treated"]],
                 smq[["SEQN", "smoker"]], diq[["SEQN", "diabetes"]]]:
        d = d.merge(part, on="SEQN", how="left")

    d = d.rename(columns={"SEQN": "seqn", "RIDAGEYR": "age", "RIAGENDR": "sex",
                          "RIDRETH3": "race", "LBXSCH": "tc", "LBDHDD": "hdl",
                          "BMXBMI": "bmi"})
    d["seqn"] = d["seqn"].astype(int)
    d["female"] = (d["sex"] == 2).astype(float)
    # RIDRETH3 == 4 is Non-Hispanic Black; the Pooled Cohort Equations are
    # specified only for White and African-American adults, so every other
    # category is scored with the white coefficients, as is common practice.
    d["black"] = (d["race"] == 4).astype(float)

    d = d.merge(read_mort(INFL / mort_file), on="seqn", how="left")

    # The scores are specified for middle-aged and older adults free of prior
    # CVD. NHANES carries no reliable prevalent-CVD flag in these files, so the
    # cohort is defined by age alone and that limitation is stated in the README.
    need = ["age", "tc", "hdl", "sbp", "bp_treated", "smoker", "diabetes", "bmi", "female"]
    d = d[(d["age"] >= 40) & (d["age"] <= 79)]
    complete = d.dropna(subset=need + ["mortstat"]).copy()

    complete.to_csv(DATA / out_name, index=False)
    print(f"{out_name}: n={len(complete)} complete adults aged 40-79 "
          f"(deaths={int(complete['mortstat'].sum())}, "
          f"CVD deaths={int(complete['cvd_death'].sum())})")
    return complete


if __name__ == "__main__":
    build("I", "MORT_2015.dat", "cvd_2015.csv")
    build("J", "MORT_2017.dat", "cvd_2017.csv")
