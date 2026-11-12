"""Post-remediation audit of numerical invariants, diagnostics and registered claims.

Exit nonzero on a failed check. Scientific limitations are reported separately
from software checks; passing this audit does not establish biological validity.
"""
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from framework_audit import legacy_ranks
from order_utils import dominance, pair_disagreement
from serotonin_controls import load as load_serotonin, score_rows

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pass_check(label, detail=""):
    print(f"PASS {label}" + (": " + detail if detail else ""))


def main():
    print("POST-REMEDIATION FRAMEWORK AUDIT")
    print("Checks concern the current local inputs and code; no source repositories are modified.\n")
    tests = subprocess.run([sys.executable, "-m", "unittest", "discover", "-s", "analysis", "-p", "test_*.py"],
                           cwd=ROOT, capture_output=True, text=True)
    if tests.returncode:
        raise AssertionError(tests.stdout + tests.stderr)
    pass_check("mathematical/native-kernel regression tests", re.search(r"Ran \d+ tests", tests.stderr).group())

    m = pd.read_csv(ROOT.parent/"selectivity/klaeger_matrix.csv", index_col=0).to_numpy(float)
    m = m[(m > 6).any(1)]
    r = legacy_ranks(m, "average")
    reversed_r = legacy_ranks(m[::-1], "average")[::-1]
    np.testing.assert_array_equal(r, reversed_r)
    pass_check("kinase row-order invariance", "all 800 parameter-median rank cells identical")
    ge = dominance(-r)
    i, j = np.triu_indices(len(r), 1)
    strict, tie_only = pair_disagreement(r, i, j)
    np.testing.assert_array_equal(strict, ~(ge | ge.T)[i,j])
    pass_check("strict-reversal/consensus agreement", f"{len(i)} pairs; tie-only cases remain comparable")
    # This is the sampler's complete poset/initializer input, before any random moves.
    encoded = rankdata(-r, axis=0, method="dense")
    transformed = np.column_stack([np.exp((-r[:,0])/100), -r[:,1]*1000, -r[:,2]+97, -r[:,3]])
    np.testing.assert_array_equal(encoded, rankdata(transformed, axis=0, method="dense"))
    pass_check("monotone-rescaling invariance", "identical encoded poset and rank-based sampler initialization")

    sampler_outputs = ["analysis/outputs/avg_rank.txt", "tp53/analysis/outputs/avg_rank.txt",
                       "inflammation/analysis/outputs/avg_rank.txt", "inflammation/analysis/outputs/sepsis_avg_rank.txt"]
    for rel in sampler_outputs:
        text = (ROOT/rel).read_text()
        lines = [line for line in text.splitlines() if line.startswith("Insertion Gibbs:")]
        assert lines and lines[-1].endswith("PASS"), rel
        last = lines[-1]
        rhat = float(re.search(r"max split-Rhat=([0-9.]+)", last).group(1))
        mcse = float(re.search(r"\(([0-9.]+) of classes\)", last).group(1))
        assert rhat < 1.01 and mcse < .0025, rel
        pass_check(f"sampler {rel}", f"max split-Rhat={rhat:g}; max relative MCSE={mcse:g}")

    prom = ["waf1","mdm2","bax","h1433s","aip1","gadd45","noxa","p53r2"]
    tp = pd.read_csv(ROOT/"tp53/data/tp53_matrix.csv").dropna(subset=prom+["am_path"])
    canonical = pd.read_csv(ROOT/"tp53/data/tp53_canonical.csv")
    joined = tp.merge(canonical[["key","canonical_avg_rank"]], on="key", validate="one_to_one")
    assert len(joined) == len(tp) and np.isfinite(joined.canonical_avg_rank).all()
    x = np.column_stack([-joined[prom].to_numpy(float), joined.am_path.to_numpy(float)])
    codes = pd.DataFrame(x).groupby(list(range(x.shape[1]))).ngroup()
    assert joined.canonical_avg_rank.groupby(codes).nunique().max() == 1
    pass_check("equivalent TP53 objects", "identical indicator vectors receive identical aggregate scores")

    matrix, measured = load_serotonin(both_targets=True)
    scores = score_rows(matrix)
    eligible = np.isfinite(scores).all(1)
    assert len(matrix) == 297 and eligible.sum() == 285
    assert np.isnan(matrix).any()
    assert not np.allclose(scores[eligible], score_rows(matrix, impute=True)[eligible])
    pass_check("serotonin missingness/eligibility", "observed-only default; 285/297 eligible; floor scenario differs")

    sys.path.insert(0, str(ROOT/"inflammation/analysis"))
    from sepsis_orientation import discovery_signs
    s = pd.read_csv(ROOT/"inflammation/data/sepsis_scores.csv")
    names = ["Inflammatory","IFNg","SRS7","SRSq19","MARS8"]
    expected = discovery_signs(s, names)
    flipped = s.copy()
    mask = flipped.cohort == "validation"
    flipped.loc[mask, "mortality"] = 1 - flipped.loc[mask, "mortality"]
    assert discovery_signs(flipped, names) == expected
    pass_check("sepsis held-out orientation", "changing validation outcomes cannot change discovery signs")
    print("\nREGISTERED MANUSCRIPT CLAIMS")
    spec = importlib.util.spec_from_file_location("check_claims", ROOT/"paper/check_claims.py")
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    mod.check()

    manifest = json.loads((ROOT/"analysis/outputs/INPUT_HASHES.json").read_text())
    for name, digest in manifest["files"].items():
        assert sha256(ROOT.parent/name) == digest, f"Stale input/code provenance: {name}"
    pass_check("source/code provenance", f"{len(manifest['files'])} current hashes; five repository revisions recorded")
    manuscript = (ROOT/"paper/manuscript.md").read_text()
    prohibited = ["no measure satisfies every requirement", "an axis independent of the four measures",
                  "rank-identical (1.00)", "the size of the constitutive component",
                  "different treat/no-treat", "while a gated smooth-hinge"]
    for phrase in prohibited:
        assert phrase not in manuscript, f"Unresolved audited claim: {phrase}"
    assert "not reconstructions of the published SRSq or MARS classifiers" in manuscript
    assert "cannot establish a calibration error" in manuscript
    pass_check("targeted manuscript claim corrections")
    print("\nSCIENTIFIC LIMITS THAT REMAIN")
    print("- Numerical diagnostics are empirical; they are not a finite-time mixing proof.")
    print("- Replicated scores can share inputs, training information, censoring and systematic bias.")
    print("- Sepsis results concern gene-set proxies, not the official SRSq/MARS classifiers.")
    print("- Sparse/selected receptor panels and heterogeneous clinical endpoints limit cross-domain interpretation.")
    print("- Registered numerical checks are not an exhaustive review of all biological claims or references.")
    print("\nRESULT: targeted numerical and reporting checks PASS; concept-level attribution remains conditional.")


if __name__ == "__main__":
    main()
