"""Hash source inputs and analysis code, including sibling-repository revisions."""
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    files = set()
    for directory in ["analysis", "tp53/analysis", "inflammation/analysis", "cvdrisk/analysis"]:
        files.update((ROOT/directory).glob("*.py"))
        files.update((ROOT/directory).glob("*.sh"))
        files.update((ROOT/directory).glob("*.R"))
    files.add(ROOT/"analysis/_poset_sampler.cpp")
    files.update((ROOT/"paper").glob("*.py"))
    files.add(ROOT/"paper/claim_registry.json")
    files.add(ROOT/"paper/manuscript.md")
    for directory in ["analysis/data", "tp53/data", "inflammation/data", "cvdrisk/data"]:
        files.update(p for p in (ROOT/directory).glob("*") if p.is_file())
    for sibling, names in {
        "selectivity": ["klaeger_matrix.csv", "davis_affinity.csv", "anastassiadis_matrix.csv", "metz_matrix.csv", "candidate_measure.py", "curate_panels.py", "davis_proteins.csv", "davis_klaeger_overlap.csv", "platform_vs_definition.py", "paper/sections/methods.tex"],
        "psychedelic-selectivity": ["data/serotonin_receptor_matrix.csv", "data/selectivity_framework_results.csv"],
        "epigenetic-clock-desiderata": ["data/clock_outputs.parquet", "data/instability_scores.parquet", "data/functional_relationships.parquet", "data/cosine_similarity.parquet", "data/d3_cell_type.parquet"],
        "methylation-biomarker-agreement": ["data/scores_GSE50660.rds", "data/scores_GSE42861.rds", "data/fasd_all_scores.rds"],
    }.items():
        files.update(ROOT.parent/sibling/name for name in names)
    manifest = {"files": {}, "repositories": {}}
    for path in sorted(files):
        key = str(path.relative_to(ROOT.parent))
        if not path.exists():
            raise FileNotFoundError(path)
        manifest["files"][key] = sha256(path)
    for path in [ROOT] + [ROOT.parent/name for name in ["selectivity", "psychedelic-selectivity", "epigenetic-clock-desiderata", "methylation-biomarker-agreement"]]:
        head = subprocess.check_output(["git", "-C", str(path), "rev-parse", "HEAD"], text=True).strip()
        status = subprocess.check_output(["git", "-C", str(path), "status", "--porcelain"], text=True)
        manifest["repositories"][path.name] = {"head": head, "dirty": bool(status)}
    target = ROOT/"analysis/outputs/INPUT_HASHES.json"
    target.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"Recorded {len(files)} source/code hashes and 5 repository revisions in {target}")


if __name__ == "__main__":
    main()
