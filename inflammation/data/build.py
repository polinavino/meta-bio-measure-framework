import pandas as pd

DATA = "/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/inflammation/data"

def xpt(name):
    return pd.read_sas(f"{DATA}/{name}", format="xport", encoding="latin1")

def read_mort(path):
    # Fixed-width, 1-indexed positions per SAS_ReadInProgramAllSurveys.sas
    colspecs = [(0, 6), (14, 15), (15, 16), (16, 19), (42, 45), (45, 48)]
    names = ["seqn", "eligstat", "mortstat", "ucod", "permth_int", "permth_exm"]
    df = pd.read_fwf(path, colspecs=colspecs, names=names,
                     na_values=[".", ".."], dtype={"ucod": str})
    df["seqn"] = df["seqn"].astype(int)
    return df[["seqn", "mortstat", "permth_int"]]

def build(cyc, suffix, mortfile, out):
    cbc = xpt(f"CBC_{suffix}.XPT")[["SEQN", "LBDLYMNO", "LBDNENO", "LBDMONO", "LBXPLTSI"]]
    crp = xpt(f"HSCRP_{suffix}.XPT")[["SEQN", "LBXHSCRP"]]
    bio = xpt(f"BIOPRO_{suffix}.XPT")[["SEQN", "LBXSAL"]]
    demo = xpt(f"DEMO_{suffix}.XPT")[["SEQN", "RIDAGEYR", "RIAGENDR"]]
    mort = read_mort(f"{DATA}/{mortfile}")

    df = demo.merge(cbc, on="SEQN", how="left") \
             .merge(crp, on="SEQN", how="left") \
             .merge(bio, on="SEQN", how="left")
    df["SEQN"] = df["SEQN"].astype(int)
    df = df.merge(mort, left_on="SEQN", right_on="seqn", how="left")

    out_df = pd.DataFrame({
        "seqn": df["SEQN"],
        "age": df["RIDAGEYR"],
        "sex": df["RIAGENDR"],
        "neut": df["LBDNENO"],
        "lymph": df["LBDLYMNO"],
        "mono": df["LBDMONO"],
        "platelet": df["LBXPLTSI"],
        "crp": df["LBXHSCRP"],
        "albumin": df["LBXSAL"],
        "mortstat": df["mortstat"],
        "permth": df["permth_int"],
    })
    path = f"{DATA}/{out}"
    out_df.to_csv(path, index=False)

    print(f"\n===== {cyc} -> {path} =====")
    print("rows:", len(out_df))
    print("non-null per column:")
    print(out_df.notna().sum().to_string())
    print("mortstat value counts:", out_df["mortstat"].value_counts(dropna=False).to_dict())
    print("permth range:", out_df["permth"].min(), "-", out_df["permth"].max())
    print("sample rows with mortality:")
    print(out_df[out_df["mortstat"].notna()].head(3).to_string())
    return out_df

build("2015-2016", "I", "MORT_2015.dat", "nhanes_2015.csv")
build("2017-2018", "J", "MORT_2017.dat", "nhanes_2017.csv")
