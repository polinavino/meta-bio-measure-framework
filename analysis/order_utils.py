"""Shared weak-order semantics and checked average ranks of linear extensions.

All input columns are concept-positive. Strict reversals require both signs;
ties versus strict comparisons are reported separately and remain comparable.
Identical indicator vectors are collapsed to equivalence classes before sampling.
The aggregate ranks those classes and gives tied objects the same result.
"""
from concurrent.futures import ThreadPoolExecutor
import hashlib
import itertools
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.stats import rankdata


def pair_disagreement(x, i, j):
    delta = np.asarray(x)[j] - np.asarray(x)[i]
    strict = (delta.min(axis=1) < 0) & (delta.max(axis=1) > 0)
    tie_only = (delta == 0).any(axis=1) & (delta != 0).any(axis=1) & ~strict
    return strict, tie_only


def dominance(x):
    x = np.asarray(x, float)
    if x.ndim != 2 or not np.isfinite(x).all():
        raise ValueError("Expected a finite objects-by-indicators matrix")
    ge = np.ones((len(x), len(x)), bool)
    for col in x.T:
        ge &= col[:, None] >= col[None, :]
    return ge


def exact_rank_draws(before):
    n = len(before)
    a, b = np.where(before)
    draws = []
    for order in itertools.permutations(range(n)):
        ranks = np.empty(n, float)
        ranks[list(order)] = np.arange(1, n + 1)
        if np.all(ranks[a] < ranks[b]):
            draws.append(ranks)
    return np.asarray(draws)


def _executable():
    source = Path(__file__).with_name("_poset_sampler.cpp")
    digest = hashlib.sha256(source.read_bytes()).hexdigest()[:16]
    binary = Path(tempfile.gettempdir()) / f"meta_bio_poset_{digest}"
    if not binary.exists():
        subprocess.run(["clang++", "-O3", "-std=c++14", str(source), "-o", str(binary)], check=True)
    return binary


def _chain(binary, before, start, burn, draws, thin, seed):
    payload = np.asarray(before, np.uint8).tobytes() + np.asarray(start, np.int32).tobytes()
    out = subprocess.run([str(binary), str(len(before)), str(burn), str(draws), str(thin), str(seed)],
                         input=payload, capture_output=True, check=True)
    return np.frombuffer(out.stdout, np.float64).reshape(draws, len(before)).copy()


def diagnostics(draws):
    """Split R-hat plus batch-means Monte Carlo SE, in class-rank units."""
    c, d, n = draws.shape
    half = d // 2
    split = np.concatenate([draws[:, :half], draws[:, -half:]], axis=0)
    within = split.var(axis=1, ddof=1).mean(axis=0)
    between = half * split.mean(axis=1).var(axis=0, ddof=1)
    var_hat = (half - 1) / half * within + between / half
    rhat = np.sqrt(np.divide(var_hat, within, out=np.ones(n), where=within > 0))
    rhat[(within == 0) & (between > 0)] = np.inf
    batches = 16
    width = d // batches
    bm = draws[:, :batches * width].reshape(c, batches, width, n).mean(axis=2)
    mcse = bm.var(axis=1, ddof=1).mean(axis=0) ** .5 / np.sqrt(c * batches)
    # Include disagreement between independently initialized chains.
    mcse = np.maximum(mcse, draws.mean(axis=1).std(axis=0, ddof=1) / np.sqrt(c))
    return rhat, mcse


def average_rank(x, label="poset", seed=1729, max_rounds=6):
    """Uniform average class ranks; fail rather than publish unchecked estimates.

    Random insertion is an exact conditional Gibbs update of an extension.
    Diagnostics are empirical checks, not a proof of finite-time mixing.
    Initializers use ranks, so strictly increasing rescalings leave the entire
    computation invariant. Four dispersed starts include extreme priorities.
    """
    x = np.asarray(x, float)
    if not np.isfinite(x).all():
        raise ValueError("Non-finite indicators must be gated or excluded explicitly")
    # Rank encoding also makes equivalence classes and class order scale invariant.
    encoded = rankdata(x, axis=0, method="dense")
    unique, inverse = np.unique(encoded, axis=0, return_inverse=True)
    n = len(unique)
    ge = dominance(unique)
    before = ge & ~ge.T
    print(f"{label}: {len(x)} objects, {n} equivalence classes; tied objects share class rank")
    if n <= 8:
        draws = exact_rank_draws(before)
        print(f"Exact enumeration: {len(draws)} linear extensions")
        return draws.mean(axis=0)[inverse]
    binary = _executable()
    ranks = rankdata(unique, axis=0)
    rng = np.random.default_rng(seed)
    weights = [np.ones(unique.shape[1]), np.exp(rng.normal(0, 3, unique.shape[1]))]
    for col in [0, unique.shape[1] - 1]:
        w = np.full(unique.shape[1], 1e-6)
        w[col] = 1
        weights.append(w)
    starts = [np.argsort(-(ranks @ w), kind="stable") for w in weights]
    # Every positive weighted rank sum is a valid extension of strict dominance.
    for start in starts:
        p = np.empty(n, int); p[start] = np.arange(n)
        a, b = np.where(before)
        assert np.all(p[a] < p[b])
    burn, draws, thin = 512, 1024, 2
    for attempt in range(max_rounds):
        with ThreadPoolExecutor(max_workers=4) as pool:
            tasks = [pool.submit(_chain, binary, before, start, burn, draws, thin, seed + k)
                     for k, start in enumerate(starts)]
            samples = np.stack([task.result() for task in tasks])
        rh, se = diagnostics(samples)
        # Relative precision allows comparison across different panel sizes.
        passed = rh.max() < 1.01 and se.max() / n < .0025
        print(f"Insertion Gibbs: 4 dispersed chains; burn={burn} sweeps; draws={draws}; thin={thin}; "
              f"max split-Rhat={rh.max():.5f}; max MCSE={se.max():.3f} ranks "
              f"({se.max()/n:.6f} of classes); {'PASS' if passed else 'RETRY'}", flush=True)
        if passed:
            mean = samples.mean(axis=(0, 1))
            a, b = np.where(before)
            assert np.all(mean[a] < mean[b])
            return mean[inverse]
        burn *= 2
        draws *= 2
    raise RuntimeError(f"{label}: sampler diagnostics failed; suppress aggregate benchmarks")
