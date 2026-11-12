"""
Run the measurement protocol on competing 10-year cardiovascular risk scores.

Mirrors ../../tp53/analysis and ../../inflammation/analysis: families, consensus
poset and incomparability, near-tie law, orientation and anchor validation, and
cross-cohort reproducibility.

The near-tie law here is computed on the consensus axis only (see near_tie), so
its decline is descriptive and does not identify a constitutive component, as in inflammation/analysis/sep_pairwise.py.

Cohorts: NHANES 2015-16 and 2017-18, adults 40-79 with complete risk factors.
Anchor: cause-specific cardiovascular death (UCOD_LEADING 001 or 005) from the
NHANES Linked Mortality File, with all-cause death reported alongside because
the CVD-specific event count is small.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "analysis"))
from order_utils import pair_disagreement, average_rank
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr, rankdata
from itertools import combinations
from sklearn.metrics import roc_auc_score

from scores import MEASURES

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
OUT = HERE / "outputs"
OUT.mkdir(exist_ok=True)
NAMES = list(MEASURES)
# Shared numerical cutoffs for sensitivity only. PCE and Framingham predict
# different endpoints, so these are not interchangeable treatment thresholds.
THRESHOLDS = [0.075, 0.10]


def score_frame(path):
    d = pd.read_csv(path)
    for k, fn in MEASURES.items():
        d[k] = fn(d)
    return d.dropna(subset=NAMES).reset_index(drop=True)


def families(d, label, fh):
    print(f"\n=== Families: Spearman among {len(NAMES)} scores ({label}, N={len(d)}) ===", file=fh)
    S = d[NAMES].corr(method="spearman")
    print("          " + "".join(f"{n:>10s}" for n in NAMES), file=fh)
    for n in NAMES:
        print(f"{n:10s}" + "".join(f"{S.loc[n,m]:+10.3f}" for m in NAMES), file=fh)
    pairs = {f"{a}-{b}": S.loc[a, b] for a, b in combinations(NAMES, 2)}
    for k, v in sorted(pairs.items(), key=lambda x: -x[1]):
        print(f"  {k:22s} {v:+.3f}", file=fh)
    return pairs


def consensus(d, label, fh):
    """Pareto consensus order: a >= b iff every score ranks a at least as high."""
    X = d[NAMES].to_numpy()
    n = len(d)
    rng = np.random.default_rng(0)
    idx = rng.choice(n, size=min(n, 600), replace=False)  # pairs scale quadratically
    Y = X[idx]
    ge = (Y[:, None, :] >= Y[None, :, :]).all(2)
    comparable = ge | ge.T
    iu = np.triu_indices(len(Y), k=1)
    inc = 1.0 - comparable[iu].mean()
    print(f"\n=== Consensus poset ({label}, N={len(Y)} subsample, "
          f"{len(iu[0]):,} pairs) ===", file=fh)
    print(f"incomparable-pair fraction = {inc:.3f}", file=fh)
    mean_rank = np.mean([rankdata(X[:, j]) for j in range(X.shape[1])], axis=0)
    d = d.copy()
    d["consensus"] = mean_rank
    for n in NAMES:
        print(f"  Spearman({n:9s}, consensus mean-rank) = "
              f"{spearmanr(d[n], d['consensus']).correlation:+.3f}", file=fh)
    return d, inc


def near_tie(d, label, fh):
    """Discordance against separation on the consensus axis.

    The separation axis is the consensus mean-rank, which is built from the same
    three scores whose discordance is being measured, so near-tie-on-consensus
    and disagreement are linked by construction as well as by the concept. No
    strong independent concept-proxy axis exists here: the only candidates in the
    feature table (age, SBP) are themselves inputs to all three scores, so an
    "independent" axis would be no less circular. This is therefore the
    consensus-coordinate view only, and the decline below is descriptive; no mathematical bound on a constitutive
    component has been established.
    """
    X = d[NAMES].to_numpy()
    rng = np.random.default_rng(1)
    idx = rng.choice(len(d), size=min(len(d), 2000), replace=False)
    Y, c = X[idx], rankdata(d["consensus"].to_numpy()[idx]) / len(idx)
    i, j = np.triu_indices(len(Y), k=1)
    sep = np.abs(c[i] - c[j])
    R = np.array([rankdata(Y[:, k]) for k in range(Y.shape[1])])
    disc = np.zeros(len(i), dtype=float)
    for a, b in combinations(range(len(NAMES)), 2):
        disc += (R[a][i] - R[a][j]) * (R[b][i] - R[b][j]) < 0
    disc /= len(list(combinations(range(len(NAMES)), 2)))
    print(f"\n=== Near-tie law ({label}, {len(i):,} pairs) ===", file=fh)
    print("  [separation axis = consensus mean-rank, built from the same three scores;",
          file=fh)
    print("   no independent concept-proxy axis is available here; descriptive association only]",
          file=fh)
    print(f"overall discordance rate: {disc.mean():.3f}", file=fh)
    qs = np.quantile(sep, [0, .2, .4, .6, .8, 1.0])
    slope = []
    for lo, hi in zip(qs[:-1], qs[1:]):
        m = (sep >= lo) & (sep <= hi if hi == qs[-1] else sep < hi)
        print(f"  |dsep| in [{lo:.3f},{hi:.3f}]: discordance={disc[m].mean():.3f} "
              f"(n={m.sum():,})", file=fh)
        slope.append(disc[m].mean())
    return slope


def anchor(d, label, fh):
    print(f"\n=== Orientation and anchor ({label}) ===", file=fh)
    res = {}
    for ev, name in [("cvd_death", "CVD death"), ("mortstat", "all-cause death")]:
        y = d[ev].to_numpy()
        if np.nansum(y) < 5:
            print(f"  {name}: only {int(np.nansum(y))} events, not scored", file=fh)
            continue
        print(f"  {name} (events={int(np.nansum(y))}):", file=fh)
        for n in NAMES + ["consensus"]:
            a = roc_auc_score(y, d[n])
            res[(ev, n)] = a
            print(f"    {n:10s} AUC={a:.3f}  {'OK' if a > 0.5 else 'MIS-ORIENTED'}", file=fh)
    return res


def boundary(d, label, fh):
    """Shared-cutoff disagreement; no inferred treatment-eligibility decisions."""
    print(f"\n=== Decision-threshold effect ({label}) ===", file=fh)
    for t in THRESHOLDS:
        calls = {n: (d[n] >= t).astype(int) for n in NAMES}
        C = np.array([calls[n] for n in NAMES])
        disagree = (C.sum(0) > 0) & (C.sum(0) < len(NAMES))
        print(f"  cutoff {t:.3f}: {disagree.mean()*100:5.1f}% of adults get a "
              f"different above/below-cutoff classification "
              f"({int(disagree.sum())} of {len(d)})", file=fh)
        same_endpoint = calls["fram_lab"] != calls["fram_off"]
        print(f"    same-endpoint Framingham lab/office discrepancy: {same_endpoint.mean()*100:.1f}%", file=fh)
    print("  PCE predicts hard ASCVD and Framingham general CVD; shared-cutoff differences mix",
          "endpoint and formula differences and do not identify calibration error or treatment changes.", file=fh)
    return None


def reproducibility(p15, p17, fh):
    print("\n=== Cross-cohort reproducibility (2015-16 vs 2017-18) ===", file=fh)
    v15 = np.array([p15[k] for k in sorted(p15)])
    v17 = np.array([p17[k] for k in sorted(p17)])
    r = spearmanr(v15, v17).correlation
    print(f"inter-score correlation structure: Spearman(2015-vec, 2017-vec) = {r:.3f} "
          f"({len(v15)} pairs); mean|dr| = {np.abs(v15-v17).mean():.3f}", file=fh)
    for k, a, b in zip(sorted(p15), v15, v17):
        print(f"  {k:22s} 2015 {a:+.3f}   2017 {b:+.3f}   d={b-a:+.3f}", file=fh)
    return r


if __name__ == "__main__":
    with open(OUT / "cvd_protocol.txt", "w") as fh:
        frames, pair_sets = {}, {}
        for cyc, f in [("2015-16", "cvd_2015.csv"), ("2017-18", "cvd_2017.csv")]:
            d = score_frame(DATA / f)
            print(f"\n{'#'*72}\n# {cyc}: N={len(d)}  "
                  f"deaths={int(d.mortstat.sum())}  CVD deaths={int(d.cvd_death.sum())}\n"
                  f"{'#'*72}", file=fh)
            for n in NAMES:
                print(f"  {n:10s} median {d[n].median()*100:5.2f}%  "
                      f"IQR {d[n].quantile(.25)*100:5.2f}-{d[n].quantile(.75)*100:5.2f}%", file=fh)
            pair_sets[cyc] = families(d, cyc, fh)
            d, _ = consensus(d, cyc, fh)
            near_tie(d, cyc, fh)
            anchor(d, cyc, fh)
            boundary(d, cyc, fh)
            frames[cyc] = d
        reproducibility(pair_sets["2015-16"], pair_sets["2017-18"], fh)
    print(open(OUT / "cvd_protocol.txt").read())
