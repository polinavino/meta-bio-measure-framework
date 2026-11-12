"""Serotonin panel coverage, missingness and implementation controls.

Untested receptors are unknown. Floor imputation is a sensitivity scenario,
never the default. The sibling 297-ligand analysis uses a quotient of pKi
values; this is explicitly compared with the kinase-style top-two difference.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr
from order_utils import pair_disagreement
from kinase_controls import neg_entropy, gini

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT.parent / "psychedelic-selectivity/data/serotonin_receptor_matrix.csv"


def score_rows(matrix, impute=False, ratio="gap"):
    rows = []
    for row in np.asarray(matrix, float):
        vals = np.where(np.isnan(row), 5., row) if impute else row[np.isfinite(row)]
        if len(vals) < 2 or vals.max() <= 6:
            rows.append([np.nan] * 4); continue
        p = vals[None, :]; ordered = np.sort(vals)[::-1]
        gap = ordered[0] - max(ordered[1], 5.)
        if ratio == "pki_quotient":
            gap = ordered[0] / ordered[1] if ordered[1] > 5 else np.nan
        rows.append([-(vals > 6).mean(), neg_entropy(np.maximum(p - 5, 0))[0], gini(p, 5)[0], gap])
    return np.asarray(rows)


def load(both_targets=False):
    df = pd.read_csv(PATH)
    matrix = df.iloc[:, 1:].to_numpy(float)
    measured = np.isfinite(matrix).sum(1)
    keep = measured >= 4
    if both_targets:
        keep &= df[["5HT2A", "5HT2B"]].notna().all(axis=1).to_numpy()
    return matrix[keep], measured[keep]


def observed_anchor():
    matrix, measured = load()
    x = score_rows(matrix)
    keep = np.isfinite(x).all(1)
    return x[keep], (matrix[keep] > 6).sum(1), measured[keep]


def describe(label, x, coverage):
    good = np.isfinite(x).all(1); x = x[good]; coverage = coverage[good]
    i, j = np.triu_indices(len(x), 1)
    disc, tied = pair_disagreement(x, i, j)
    print(f"{label}: eligible/defined={len(x)}; incomparable={disc.mean():.4f}; tie-only={tied.mean():.4f}")
    for k, name in enumerate(["S-score", "entropy", "Gini", "ratio"]):
        rho = spearmanr(x[:, k], coverage).statistic if len(np.unique(coverage)) > 1 else np.nan
        print(f"  {name:8s}: rho(score, measured receptors)={rho:+.3f}")
    print("  Spearman matrix (concept-positive):")
    print(pd.DataFrame(x, columns=["S", "entropy", "Gini", "ratio"]).corr("spearman").round(3).to_string())


def main():
    for both in [False, True]:
        m, measured = load(both)
        print(f"\n=== {'2A/2B covered subset' if both else 'broader panel'}: N={len(m)}, "
              f"tested receptors={measured.min()}..{measured.max()}, available columns={m.shape[1]} ===")
        observed = score_rows(m)
        imputed = score_rows(m, impute=True)
        quotient = score_rows(m, ratio="pki_quotient")
        describe("Observed-only, log gap", observed, measured)
        describe("Floor-imputed scenario, log gap", imputed, measured)
        describe("Observed-only, source pKi quotient", quotient, measured)
        common = np.isfinite(observed).all(1) & np.isfinite(imputed).all(1)
        print("Paired observed/floor-imputed ranking correlations:",
              [round(spearmanr(observed[common, k], imputed[common, k]).statistic, 4) for k in range(4)])
        # A fixed observed coverage count does not harmonize receptor identity.
        for count in [4, 5, 6]:
            keep = measured == count
            if keep.sum() >= 10:
                describe(f"Exactly {count} measured receptors (identities still vary)", observed[keep], measured[keep])
        patterns, counts = np.unique(np.isfinite(m), axis=0, return_counts=True)
        biggest = np.argmax(counts); same = np.all(np.isfinite(m) == patterns[biggest], axis=1)
        print(f"Largest identical observed receptor panel: n={same.sum()}; receptors={patterns[biggest].sum()}")
        if same.sum() >= 10:
            describe("Identical observed receptor panel", observed[same], measured[same])
    print("\nActive count is a shared-input proxy, not an external anchor. Sparse panels and selection preclude a general serotonin family claim.")


if __name__ == "__main__":
    main()
