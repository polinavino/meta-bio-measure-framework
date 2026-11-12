"""Current kinase implementations, criteria checks, eligibility and panel recovery.

The candidate is the November-2026 gated exponential entropy. Historical hinge
experiments are not current validation. Formula sensitivity and random-panel
recovery are empirical, distinct from the construction guarantees in the source
paper (selectivity/paper/sections/methods.tex, Propositions 1 and 2).
"""
from pathlib import Path
import hashlib
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr
from order_utils import pair_disagreement

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parent / "selectivity"


def neg_entropy(w):
    total = w.sum(1, keepdims=True)
    p = w / np.where(total > 0, total, 1)
    out = (p * np.log2(np.where(p > 0, p, 1))).sum(1)
    return np.where(total[:, 0] > 0, out, np.nan)


def gini(m, floor):
    s = np.sort(np.maximum(m - floor, 0), axis=1)
    n = s.shape[1]; t = s.sum(1)
    out = 2 * (s * np.arange(1, n + 1)).sum(1) / (n * np.where(t > 0, t, 1)) - (n + 1) / n
    return np.where(t > 0, out, np.nan)


def candidate(m, floor, tau, scale=1):
    if scale < 1 or tau <= floor:
        raise ValueError("Weak-addition guarantee requires T >= tau - floor > 0")
    out = neg_entropy(np.where(m > floor, np.exp((m - floor) / (scale * (tau - floor))), 0))
    return np.where(m.max(1) > tau, out, np.nan)


def scores(m, floor, tau):
    ordered = np.sort(m, axis=1)[:, ::-1]
    ka = 10.0 ** (m - m.max(1, keepdims=True))
    return {
        "S-score": -(m > tau).mean(1),
        "log entropy": neg_entropy(np.maximum(m - floor, 0)),
        "Gini": gini(m, floor),
        "gap ratio": ordered[:, 0] - np.maximum(ordered[:, 1], floor),
        "Ka excluded": neg_entropy(np.where(m > floor, ka, 0)),
        "Ka included": neg_entropy(ka),
        "candidate": candidate(m, floor, tau),
    }


def load():
    aff = pd.read_csv(SOURCE / "davis_affinity.csv")
    davis = aff.pivot(index="Drug_Index", columns="Protein_Index", values="Affinity").fillna(5).to_numpy(float)
    return {"Davis": (davis, 5., 6.),
            "Klaeger": (pd.read_csv(SOURCE / "klaeger_matrix.csv", index_col=0).to_numpy(float), 5., 6.),
            "Anastassiadis": (pd.read_csv(SOURCE / "anastassiadis_matrix.csv", index_col=0).to_numpy(float), 20., 50.),
            "Metz": (pd.read_csv(SOURCE / "metz_matrix.csv", index_col=0).to_numpy(float), 4., 6.)}


def fraction_incomparable(x):
    x = x[np.isfinite(x).all(1)]
    i, j = np.triu_indices(len(x), 1)
    return pair_disagreement(x, i, j)[0].mean(), len(x)


def recovery(name, m, floor, tau, repetitions=80):
    full = scores(m, floor, tau)
    eligible = m.max(1) > tau
    ids = np.flatnonzero(eligible)
    rng = np.random.default_rng(42)
    i, j = np.triu_indices(len(ids), 1)
    if len(i) > 20000:
        sel = rng.choice(len(i), 20000, replace=False); i, j = i[sel], j[sel]
    i, j = ids[i], ids[j]
    names = ["log entropy", "Ka excluded", "candidate"]
    if name == "Anastassiadis":
        names.remove("Ka excluded")  # Association constants have no literal %inhibition interpretation.
    print("Random subpanels: mean eligibility retention; conditional Spearman; unconditional pair recovery")
    print("Pair denominator: sampled full-panel eligible pairs with a strict full-score order; "
          "lost eligibility=0 credit, new tie=0.5; intervals are panel 2.5/97.5% quantiles")
    for size in sorted(set([20, 50, min(100, m.shape[1]), m.shape[1]])):
        cov = []; rho = {k: [] for k in names}; pair = {k: [] for k in names}
        for _ in range(repetitions if size < m.shape[1] else 1):
            sub = m[:, rng.choice(m.shape[1], size, replace=False)]
            ok = sub.max(1) > tau
            both = eligible & ok
            cov.append(both.sum() / eligible.sum())
            sub_scores = scores(sub, floor, tau)
            for key in names:
                rho[key].append(spearmanr(full[key][both], sub_scores[key][both]).statistic)
                sf = np.sign(full[key][i] - full[key][j])
                ss = np.sign(sub_scores[key][i] - sub_scores[key][j])
                valid = sf != 0
                credit = np.where(ss == sf, 1., np.where(ss == 0, .5, 0.)) * (ok[i] & ok[j])
                pair[key].append(credit[valid].mean())
        for key in names:
            lo, hi = np.quantile(rho[key], [.025, .975])
            print(f"  targets={size:3d} {key:12s}: retention={np.mean(cov):.3f}; "
                  f"rho={np.mean(rho[key]):.3f} [{lo:.3f},{hi:.3f}]; pair={np.mean(pair[key]):.3f}")


def main():
    print("CURRENT KINASE CONTROLS — concept-positive scores, average ties, explicit eligibility")
    print("Source candidate SHA256:", hashlib.sha256((SOURCE / "candidate_measure.py").read_bytes()).hexdigest())
    print("Panel curation: Davis 433->363 constructs/reference targets; Klaeger 343->217; Anastassiadis 300->291; Metz 172")
    for name, (m, floor, tau) in load().items():
        eligible = m.max(1) > tau; gated = m[eligible]
        print(f"\n=== {name}: {m.shape[0]} compounds x {m.shape[1]} targets; eligible={eligible.sum()}; "
              f"floor={floor:g}; tau={tau:g}; T={tau-floor:g} ===")
        sc = scores(gated, floor, tau)
        print("Fixed-default implementation sensitivity; raw scores, ties retained")
        for ent in ["log entropy", "Ka excluded", "Ka included"]:
            if name == "Anastassiadis" and ent != "log entropy":
                continue
            x = np.column_stack([sc["S-score"], sc[ent], sc["Gini"], sc["gap ratio"]])
            frac, n = fraction_incomparable(x)
            print(f"  {ent:12s}: n={n}; incomparable={frac:.4f}; "
                  f"rho(ent,S)={spearmanr(sc[ent],sc['S-score']).statistic:+.3f}; "
                  f"rho(ent,ratio)={spearmanr(sc[ent],sc['gap ratio']).statistic:+.3f}")
        print("Eligibility-threshold sensitivity: scores compared on intersection, with T changed to tau-floor")
        offsets = [-.5, 0, .5, 1.] if name != "Anastassiadis" else [-10, 0, 10, 20]
        for off in offsets:
            t = tau + off; common = eligible & (m.max(1) > t)
            corr = spearmanr(candidate(m[common], floor, tau), candidate(m[common], floor, t)).statistic
            print(f"  tau={t:g}: eligible={(m.max(1)>t).sum()}; common={common.sum()}; candidate rho={corr:.4f}")
        # Empirical checks supplement, and do not replace, the source construction proofs.
        baseline = candidate(gated, floor, tau)
        weak = np.column_stack([gated, np.full(len(gated), tau)])
        violations = np.sum(candidate(weak, floor, tau) > baseline + 1e-10)
        synthetic = np.full((20, 30), floor)
        for row in range(20): synthetic[row, :row + 1] = tau + 1
        orientation = np.all(np.diff(candidate(synthetic, floor, tau)) <= 1e-10)
        # Choose two gates in an open interval containing no observed readout.
        above = gated[gated > floor].min()
        g1, g2 = floor + (above - floor) * np.array([.25, .75])
        T = tau - floor
        h1 = neg_entropy(np.where(gated > g1, np.exp((gated-g1)/T), 0))
        h2 = neg_entropy(np.where(gated > g2, np.exp((gated-g2)/T), 0))
        print(f"Criteria checks: D1 explicit max>tau eligibility; D2 no-crossing max score change={np.max(np.abs(h1-h2)):.2e}; "
              f"D3 observed weak-addition violations={violations}; D4 equal-strength orientation={orientation}")
        assert violations == 0 and orientation and np.allclose(h1, h2, atol=1e-10)
        recovery(name, m, floor, tau)
    print("\nConstruction guarantees are conditional on the stated domain and T bound; recovery is a separate empirical result.")


if __name__ == "__main__":
    main()
