"""Base-set, cohort-bootstrap and coverage controls across the demonstrations.

Bootstrap units are objects/participants, never the dependent pairs of objects
or pairs of indicators. Reproducibility does not identify concept-level truth.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr
from order_utils import pair_disagreement
from framework_audit import legacy_ranks

ROOT = Path(__file__).resolve().parents[1]


def base_sets(label, x, names):
    x = np.asarray(x, float)
    i, j = np.triu_indices(len(x), 1)
    full, ties = pair_disagreement(x, i, j)
    print(f"\n{label}: objects={len(x)}, full incomparable={full.mean():.4f}; tie-only={ties.mean():.4f}")
    for k, name in enumerate(names):
        keep = [a for a in range(len(names)) if a != k]
        disc, _ = pair_disagreement(x[:, keep], i, j)
        assert np.all(~disc | full)  # Removing a constraint cannot create a reversal.
        print(f"  omit {name:16s}: incomparable={disc.mean():.4f}")
    ranks = rankdata(x, axis=0)
    for a in range(len(names)):
        for b in range(a):
            if np.array_equal(ranks[:, a], ranks[:, b]):
                print(f"  EXACT same ordinal family: {names[a]} / {names[b]}")
    return full.mean()


def cohort_bootstrap(label, a, b, reps=400):
    iu = np.triu_indices(a.shape[1], 1)
    def vector(x):
        return pd.DataFrame(x).corr("spearman").to_numpy()[iu]
    va, vb = vector(a), vector(b)
    print(f"\n{label}: N={len(a)}/{len(b)}; indicator pairs={len(va)}; "
          f"matrix-vector rho={spearmanr(va,vb).statistic:.4f}; mean|delta r|={np.abs(va-vb).mean():.4f}")
    rng = np.random.default_rng(2026)
    rho, delta = [], []
    for _ in range(reps):
        aa = vector(a[rng.integers(len(a), size=len(a))])
        bb = vector(b[rng.integers(len(b), size=len(b))])
        rho.append(spearmanr(aa, bb).statistic)
        delta.append(np.abs(aa-bb).mean())
    print("  Independent participant bootstrap, 95% percentile intervals:")
    print(f"    matrix-vector rho [{np.quantile(rho,.025):.4f}, {np.quantile(rho,.975):.4f}]")
    print(f"    mean|delta r| [{np.quantile(delta,.025):.4f}, {np.quantile(delta,.975):.4f}]")
    print("  This quantifies cohort sampling uncertainty conditional on the chosen formulas; shared bias can replicate.")


def main():
    m = pd.read_csv(ROOT.parent / "selectivity/klaeger_matrix.csv", index_col=0).to_numpy(float)
    m = m[(m > 6).any(1)]
    base_sets("Kinase parameter-median ranks (gated population)", -legacy_ranks(m, "average"), ["S-score","entropy","Gini","ratio"])
    prom = ["waf1","mdm2","bax","h1433s","aip1","gadd45","noxa","p53r2"]
    d = pd.read_csv(ROOT / "tp53/data/tp53_matrix.csv").dropna(subset=prom+["am_path"])
    base_sets("TP53", np.column_stack([-d[prom].to_numpy(float), d.am_path.to_numpy(float)]), prom+["AlphaMissense"])
    print("TP53 experiment/prediction blocks are approximate correlation clusters, not asserted exact ordinal equivalence classes.")
    idx = ["NLR","PLR","MLR","SII","SIRI","CRP","CAR"]
    frames = [pd.read_csv(ROOT/f"inflammation/data/indices_{year}.csv").dropna(subset=idx) for year in [2015,2017]]
    for year, frame in zip([2015,2017], frames):
        print(f"NHANES {year}: complete-index population={len(frame)}; with mortality={frame.mortstat.notna().sum()}")
    sub = np.random.default_rng(0).choice(len(frames[0]), 600, replace=False)
    base_sets("NHANES complete-index subsample", frames[0].iloc[sub][idx].to_numpy(float), idx)
    corr = frames[0][["CRP","CAR"]].corr("spearman").iloc[0,1]
    print(f"CRP/CAR full-population rho={corr:.8f}; exactly identical rank vectors="
          f"{np.array_equal(rankdata(frames[0].CRP),rankdata(frames[0].CAR))}")
    cohort_bootstrap("NHANES index structure", frames[0][idx].to_numpy(float), frames[1][idx].to_numpy(float))
    sig = ["Inflammatory","IFNg","SRS7","SRSq19","MARS8"]
    s = pd.read_csv(ROOT/"inflammation/data/sepsis_scores.csv").dropna(subset=sig+["mortality"])
    discovery = s[s.cohort=="discovery"]
    signs = np.array([1 if spearmanr(discovery[k],discovery.mortality).statistic>=0 else -1 for k in sig])
    a = discovery[sig].to_numpy(float)*signs
    b = s[s.cohort=="validation"][sig].to_numpy(float)*signs
    print("\nSepsis: fit orientation on discovery and freeze on validation:",dict(zip(sig,signs)))
    cohort_bootstrap("Sepsis fixed-orientation structure", a, b)
    print("\nReproducible differences are descriptive. Constitutive attribution additionally needs independent evidence of admissibility and controls for systematic artifacts.")


if __name__ == "__main__":
    main()
