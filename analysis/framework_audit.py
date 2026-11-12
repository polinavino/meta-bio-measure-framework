"""Sensitivity checks prompted by the kinase paper's post-submission revisions.

Read-only with respect to source data and existing analyses. Run from any directory:
    python analysis/framework_audit.py --sampler-check > analysis/outputs/framework_audit.txt

These diagnostics do not replace the published analyses or establish a true
concept order. The legacy reconstruction deliberately retains its zero-weight
score conventions, so tie effects can be isolated from undefined-score handling.
"""
from pathlib import Path
import argparse
import hashlib
import subprocess

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parent / "selectivity"


def neg_entropy(weights):
    total = weights.sum(1, keepdims=True)
    p = weights / np.where(total == 0, 1.0, total)
    return (p * np.log2(np.where(p > 0, p, 1.0))).sum(1)


def gini(matrix, baseline):
    s = np.sort(np.maximum(matrix - baseline, 0), axis=1)
    total = s.sum(1)
    n = s.shape[1]
    value = 2 * (s * np.arange(1, n + 1)).sum(1) / (n * np.where(total == 0, 1, total)) - (n + 1) / n
    return np.where(total == 0, 0.0, value)


def legacy_ranks(matrix, method):
    """Reconstruct avg_rank.py's parameter-median ranks on its gated set."""
    n = len(matrix)

    def rank(score):
        return n - rankdata(score, method=method) + 1

    def legacy_entropy(baseline):
        # Preserve the original epsilon convention as well as the formula.
        w = np.maximum(matrix - baseline, 0)
        total = w.sum(1, keepdims=True)
        p = w / np.where(total == 0, 1e-10, total)
        return (p * np.where(p > 0, np.log2(p + 1e-10), 0)).sum(1)

    ordered = np.sort(matrix, axis=1)[:, ::-1]
    return np.column_stack([
        np.median([rank(-(matrix > t).mean(1)) for t in np.arange(5.5, 8.25, .25)], axis=0),
        np.median([rank(legacy_entropy(b)) for b in np.arange(5, 6.75, .25)], axis=0),
        np.median([rank(gini(matrix, b)) for b in np.arange(5, 6.75, .25)], axis=0),
        np.median([rank(ordered[:, 0] - np.maximum(ordered[:, k], 5)) for k in range(1, 6)], axis=0),
    ])


def pair_masks(scores):
    i, j = np.triu_indices(len(scores), 1)
    delta = scores[j] - scores[i]
    signs = np.sign(delta)
    strict = (delta.min(1) < 0) & (delta.max(1) > 0)
    sign_mismatch = (signs != signs[:, :1]).any(1)
    return i, j, strict, sign_mismatch


def describe(label, scores):
    _, _, strict, mismatch = pair_masks(scores)
    print(f"{label}: incomparable/strict reversal={strict.mean():.4f}; "
          f"sign-mismatch={mismatch.mean():.4f}; tie-only={np.mean(mismatch & ~strict):.4f}")


def sampler_check(x, start, saved, median8):
    """Use the existing TP53 chain budget from a different valid extension."""
    n = len(x)
    ge = np.ones((n, n), bool)
    for col in x.T:
        ge &= col[:, None] >= col[None, :]
    comparable = ge | ge.T
    order = list(np.argsort(start))
    rng = np.random.default_rng(1)
    rank_sum = np.zeros(n)
    positions = np.arange(1, n + 1)
    burn, nsamp, gap = 1_000_000, 2000, 2000
    for step in range(burn + nsamp * gap):
        i = rng.integers(n - 1)
        a, b = order[i], order[i + 1]
        if not comparable[a, b] and rng.random() < .5:
            order[i], order[i + 1] = b, a
        if step >= burn and (step - burn) % gap == 0:
            rank_sum[order] += positions
    result = rank_sum / nsamp
    print("Alternate-start chain: same 1,000,000 burn-in + 2,000 x 2,000 steps as TP53 script")
    print(f"rho(alternate-start mean ranks, stored canonical)={spearmanr(result, saved).statistic:.6f}")
    print(f"rho(alternate-start mean ranks, alternate start)={spearmanr(result, start).statistic:.6f}")
    print(f"rho(median8, alternate-start mean ranks)={spearmanr(-median8, result).statistic:.6f}")
    print(f"mean absolute difference from stored canonical={np.abs(result - saved).mean():.2f} ranks")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sampler-check", action="store_true", help="also run a full-budget TP53 chain from a different valid start")
    args = parser.parse_args()
    path = SOURCE / "klaeger_matrix.csv"
    print("SOURCE PROVENANCE")
    print("selectivity HEAD:", subprocess.check_output(["git", "-C", str(SOURCE), "rev-parse", "HEAD"], text=True).strip())
    print("comparison baseline: ad5cdd2 (2026-09-21, 'uploaded to wiley')")
    print("klaeger_matrix.csv SHA256:", hashlib.sha256(path.read_bytes()).hexdigest())
    matrix = pd.read_csv(path, index_col=0).to_numpy(float)
    active = (matrix > 6).sum(1)
    keep = active > 0
    m = matrix[keep]
    print(f"matrix={matrix.shape}; gated compounds={keep.sum()}")

    print("\n1. ROW ORDER AND TIES (same legacy parameter sweeps)")
    ordinal = legacy_ranks(m, "ordinal")
    average = legacy_ranks(m, "average")
    reversed_ordinal = legacy_ranks(m[::-1], "ordinal")[::-1]
    reversed_average = legacy_ranks(m[::-1], "average")[::-1]
    describe("legacy ordinal", ordinal)
    describe("tie-preserving average", average)
    describe("reversed-row ordinal, mapped back", reversed_ordinal)
    for label, a, b in [("ordinal", ordinal, reversed_ordinal), ("average", average, reversed_average)]:
        _, _, inc_a, _ = pair_masks(a)
        _, _, inc_b, _ = pair_masks(b)
        print(f"{label}: rows with changed measure ranks={np.any(a != b, axis=1).sum()}; "
              f"pairs changing comparability={np.count_nonzero(inc_a != inc_b)} / {len(inc_a)}")

    print("\n2. FORMULA IMPLEMENTATION (fixed default parameters, ties preserved)")
    s = -(m > 6).mean(1)
    h_log = neg_entropy(np.maximum(m - 5, 0))
    ka = 10.0 ** (m - m.max(1, keepdims=True))
    h_ka_excluded = neg_entropy(np.where(m > 5, ka, 0))
    h_ka_included = neg_entropy(ka)
    ordered = np.sort(m, axis=1)[:, ::-1]
    ratio = ordered[:, 0] - np.maximum(ordered[:, 1], 5)
    for label, h in [("baseline-subtracted log entropy", h_log),
                     ("Ka entropy, floor excluded", h_ka_excluded),
                     ("Ka entropy, floor included", h_ka_included)]:
        x = np.column_stack([s, h, gini(m, 5), ratio])
        describe(label, x)
        print(f"  rho(entropy, S)={spearmanr(h, s).statistic:+.4f}; "
              f"rho(entropy, ratio)={spearmanr(h, ratio).statistic:+.4f}")
    print(f"rho(n_active, S-score at tau=6)={spearmanr(active[keep], s).statistic:+.4f} "
          "(same information; not an independent anchor)")

    print("\n3. ACTIVE-COUNT SEPARATION, EXCLUDING THE S-SCORE")
    # Match kinase_pairwise.py: compute ranks over all compounds, then gate.
    r = legacy_ranks(matrix, "average")[keep]
    for label, scores in [("all four measures", r), ("entropy/Gini/ratio only", r[:, 1:])]:
        i, j, strict, mismatch = pair_masks(scores)
        separation = np.abs(active[keep][i] - active[keep][j])
        print(label)
        for low, high in [(0, 0), (1, 2), (3, 5), (6, 15), (16, 999)]:
            mask = (separation >= low) & (separation <= high)
            print(f"  separation {low}-{high}: strict={strict[mask].mean():.4f}; "
                  f"sign-mismatch={mismatch[mask].mean():.4f}; pairs={mask.sum()}")

    print("\n4. TIE-ONLY DISCORDANCE IN OTHER DOMAINS")
    prom = ["waf1", "mdm2", "bax", "h1433s", "aip1", "gadd45", "noxa", "p53r2"]
    tp = pd.read_csv(ROOT / "tp53/data/tp53_matrix.csv").dropna(subset=prom + ["am_path"]).reset_index(drop=True)
    describe("TP53 8-promoter panel", -tp[prom].to_numpy(float))
    idx = ["NLR", "PLR", "MLR", "SII", "SIRI", "CRP", "CAR"]
    inf = pd.read_csv(ROOT / "inflammation/data/indices_2015.csv").dropna(subset=idx).reset_index(drop=True)
    sub = np.random.default_rng(0).choice(len(inf), size=2000, replace=False)
    describe("NHANES 7-index panel (existing 2000-person subset)", inf.iloc[sub][idx].to_numpy(float))

    print("\n5. TP53 SAMPLER INITIALIZATION DIAGNOSTIC (not a convergence proof)")
    x = np.column_stack([-tp[prom].to_numpy(float), tp.am_path.to_numpy(float)])
    start = np.empty(len(tp), float)
    start[np.argsort(-x.mean(1))] = np.arange(1, len(tp) + 1)
    canonical = pd.read_csv(ROOT / "tp53/data/tp53_canonical.csv")
    aligned = tp[["key"]].merge(canonical[["key", "canonical_avg_rank"]], on="key", validate="one_to_one")
    saved = aligned.canonical_avg_rank.to_numpy(float)
    assert len(saved) == len(tp)
    assert np.isfinite(saved).all()
    median8 = -tp[prom].median(axis=1).to_numpy(float)
    print(f"rho(saved canonical, common starting order)={spearmanr(saved, start).statistic:.6f}")
    print(f"mean absolute displacement from starting order={np.abs(saved - start).mean():.2f} ranks")
    print(f"rho(median8, starting order)={spearmanr(-median8, start).statistic:.6f}")
    print(f"rho(median8, saved canonical)={spearmanr(-median8, saved).statistic:.6f}")
    # A monotone rescaling preserves every dominance edge but changes this initializer.
    scaled = x.copy()
    scaled[:, -1] *= 1000
    alt_start = rankdata(-scaled.mean(1), method="ordinal")
    print(f"rho(common start, start after AlphaMissense x1000)={spearmanr(start, alt_start).statistic:.6f}")
    print("Both starts respect the same dominance order. Re-run chains from dispersed starts, "
          "with monotone-rescaling checks and Monte Carlo errors, before treating agreement "
          "between same-start chains as convergence.")
    if args.sampler_check:
        sampler_check(x, alt_start, saved, median8)


if __name__ == "__main__":
    main()
