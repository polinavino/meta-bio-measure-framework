"""Check registered manuscript quantities against regenerated analysis outputs.

This verifies the named numerical claims, not every sentence or its scientific
interpretation. The registry stores locations/patterns and rounding tolerance;
the expected values are always read from the current source outputs.
"""
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]


def extract(text, pattern, index=0):
    matches = list(re.finditer(pattern, text, flags=re.S))
    if len(matches) <= index:
        raise AssertionError(f"Missing registered pattern {pattern!r}")
    return float(matches[index].group(1).replace(",", "").replace("−", "-").rstrip("."))


def check():
    manuscript = (ROOT/"paper/manuscript.md").read_text()
    registry = json.loads((ROOT/"paper/claim_registry.json").read_text())
    for entry in registry:
        source = (ROOT/entry["source"]).read_text()
        expected = extract(source, entry["source_pattern"], entry.get("source_index", 0))
        expected *= entry.get("source_scale", 1)
        actual = extract(manuscript, entry["manuscript_pattern"], entry.get("manuscript_index", 0))
        if abs(actual - expected) > entry.get("tolerance", .00051):
            raise AssertionError(f"{entry['name']}: manuscript {actual}, source {expected}")
        print(f"PASS {entry['name']}: manuscript={actual:g}; source={expected:g}")
    print(f"Verified {len(registry)} registered quantitative claims (rounding-aware).")
    return len(registry)


if __name__ == "__main__":
    check()
