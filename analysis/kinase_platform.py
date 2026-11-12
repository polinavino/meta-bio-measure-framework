"""Matched-compound, harmonized-panel cross-assay comparison.

Reuses the kinase revision's curated target aliases and CID-resolved compound
overlap without importing its scripts (which write into the sibling repo).
Reports actual pair-state overlap, not only correlations of scalar summaries.
"""
import ast
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from framework_audit import legacy_ranks
from order_utils import pair_disagreement

SOURCE = Path(__file__).resolve().parents[2] / "selectivity"


def main():
    source = ast.parse((SOURCE/"platform_vs_definition.py").read_text())
    alias = next(ast.literal_eval(node.value) for node in source.body
                 if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "ALIAS" for t in node.targets))
    overlap = pd.read_csv(SOURCE/"davis_klaeger_overlap.csv")
    aff = pd.read_csv(SOURCE/"davis_affinity.csv")
    davis = aff.pivot(index="Drug_Index", columns="Protein_Index", values="Affinity").fillna(5.)
    klaeger = pd.read_csv(SOURCE/"klaeger_matrix.csv", index_col=0)
    names = pd.read_csv(SOURCE/"davis_proteins.csv").set_index("Protein_Index").Gene_Name.reindex(davis.columns)
    genes = names.str.replace(r"\(.*$", "", regex=True).str.replace(r"-(phosphorylated|nonphosphorylated)$", "", regex=True).map(lambda g: alias.get(g,g))
    dgene = davis.T.groupby(genes.values).max().T
    kcol = {sym: col for col in klaeger.columns for sym in col.split(";")}
    common = sorted(set(dgene.columns) & set(kcol))
    a = dgene.loc[overlap.Davis_Drug_Index, common].to_numpy(float)
    b = klaeger.loc[overlap.Klaeger_Name, [kcol[g] for g in common]].to_numpy(float)
    keep = (a > 6).any(1) & (b > 6).any(1)
    a, b = a[keep], b[keep]
    print(f"CID-resolved overlap={len(overlap)}; shared curated genes={len(common)}; eligible on both={len(a)}")
    ra, rb = legacy_ranks(a, "average"), legacy_ranks(b, "average")
    for k, name in enumerate(["S-score","entropy","Gini","ratio"]):
        print(f"  same-definition cross-platform rho, {name:8s} = {spearmanr(ra[:,k],rb[:,k]).statistic:.4f}")
    i, j = np.triu_indices(len(a), 1)
    da, ta = pair_disagreement(ra, i, j)
    db, tb = pair_disagreement(rb, i, j)
    intersection, union = da & db, da | db
    print(f"Pairs={len(i)}; incomparable Davis={da.mean():.4f}, Klaeger={db.mean():.4f}")
    jaccard = intersection.sum()/union.sum()
    print(f"Incomparable-set Jaccard={jaccard:.4f}; intersection={intersection.sum()}, union={union.sum()}")
    jointly_comparable = ~da & ~db
    dira = np.sign((ra[j]-ra[i]).sum(1))
    dirb = np.sign((rb[j]-rb[i]).sum(1))
    strict = jointly_comparable & (dira != 0) & (dirb != 0)
    print(f"Direction agreement on pairs comparable and strict on both={np.mean(dira[strict]==dirb[strict]):.4f} (pairs={strict.sum()})")
    rng = np.random.default_rng(2026); boot=[]
    for _ in range(2000):
        counts = np.bincount(rng.integers(len(a), size=len(a)), minlength=len(a))
        weight = counts[i]*counts[j]
        if weight[union].sum():
            boot.append(weight[intersection].sum()/weight[union].sum())
    print(f"Matched-compound bootstrap Jaccard 95% CI [{np.quantile(boot,.025):.4f}, {np.quantile(boot,.975):.4f}]")
    print("The compound and gene panels are harmonized; assay technologies, censoring and errors still differ.")
    print("This directly measures pair-set replication for this selected overlap, without attributing it uniquely to the concept.")


if __name__ == "__main__":
    main()
