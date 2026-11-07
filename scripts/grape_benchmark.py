"""Run a small known-truth benchmark with GRAPE's external validation harness.

Requires a separately built GRAPE checkout and its validation dependencies.
No inference or scoring code is included in ngesh's runtime package.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import statistics
import subprocess
import sys
import tempfile

from ngesh import simulate


def benchmark(grape: Path, replicates: int) -> dict:
    sys.path.insert(0, str(grape / "validation"))
    import validate
    from ete3 import Tree

    executable = grape / "build/grape"
    scenarios = {
        "lexical-baseline": ("lexical", {}),
        "lexical-relaxed-clock": ("lexical", {"characters": {"clock_sigma": 0.7}}),
        "lexical-borrowing": ("lexical", {"horizontal": {"borrowing_rate": 0.02}}),
        "lexical-synonymy": ("lexical", {"characters": {"synonym_duration": 0.3}}),
        "morphology": ("morphology", {}),
        "genecontent": ("genecontent", {}),
        "stemma": ("stemma", {}),
    }
    records = []
    with tempfile.TemporaryDirectory() as tmp:
        for scenario, (preset, overrides) in scenarios.items():
            for replicate in range(replicates):
                config = {"tree": {"n_taxa": 20}, "characters": {"n_characters": 150}}
                for section, values in overrides.items():
                    config.setdefault(section, {}).update(values)
                out = Path(tmp) / f"{scenario}-{replicate}"
                simulate(
                    preset=preset,
                    seed=f"benchmark-v1:{scenario}:{replicate}",
                    config=config,
                    out=out,
                )
                truth = Tree((out / "truth.nwk").read_text(), format=1)
                source = out / "data.tsv"
                cognates = validate.grape.read_cognate_file(str(source))
                represented = {taxon for taxon, _ in cognates}
                truth_taxa = set(truth.get_leaf_names())
                if represented != truth_taxa:
                    raise ValueError(
                        f"{scenario}: TSV omits {sorted(truth_taxa - represented)} truth taxa"
                    )
                run = subprocess.run(
                    [
                        str(executable),
                        "--graph",
                        "adjusted",
                        "--community",
                        "louvain",
                        "--strategy",
                        "fixed",
                        "--seed",
                        "42",
                        str(source),
                    ],
                    capture_output=True,
                    text=True,
                    check=True,
                    timeout=120,
                )
                inferred = Tree(run.stdout.strip(), format=1)
                distances = validate.grape.compute_distance_matrix(
                    cognates, missing_data="max_dist"
                )
                taxa = sorted(represented)
                methods = {
                    "grape": inferred,
                    "nj": validate.dt_mod.nj(distances, taxa),
                    "upgma": validate.dt_mod.upgma(distances, taxa),
                }
                for method, tree in methods.items():
                    if set(tree.get_leaf_names()) != truth_taxa:
                        raise ValueError("inference changed the taxon roster")
                    scores = validate.clade_scores(tree, truth)
                    records.append(
                        {
                            "scenario": scenario,
                            "replicate": replicate,
                            "method": method,
                            "normalized_rf": validate.norm_rf(tree, truth),
                            "clade_recall": scores["recall"],
                            "clade_precision": scores["precision"],
                            "clade_f1": scores["f1"],
                        }
                    )
    summary = {}
    for scenario in scenarios:
        summary[scenario] = {}
        for method in ["grape", "nj", "upgma"]:
            rows = [r for r in records if r["scenario"] == scenario and r["method"] == method]
            summary[scenario][method] = {
                key: statistics.mean(r[key] for r in rows)
                for key in ["normalized_rf", "clade_recall", "clade_precision", "clade_f1"]
            }
    revision = subprocess.check_output(
        ["git", "-C", str(grape), "rev-parse", "HEAD"], text=True
    ).strip()
    return {
        "schema": "ngesh-grape-benchmark-v1",
        "grape_revision": revision,
        "replicates_per_scenario": replicates,
        "simulation": "20 taxa, 150 characters; parameters from presets and listed scenario overrides",
        "methods": "GRAPE native adjusted/louvain/fixed; NJ and UPGMA from GRAPE validation harness",
        "distance": "GRAPE legacy average synonym strategy and max_dist missing policy for baselines",
        "score": "GRAPE validation harness unrooted normalized RF and rooted clade precision/recall/F1",
        "limitations": [
            "small exploratory simulation; no confidence intervals or significance claims",
            "gene-content sparse TSV conflates known absence and missingness",
            "no reticulation network scoring; tree truth only",
        ],
        "scenarios": scenarios,
        "summary": summary,
        "records": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grape", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=5)
    args = parser.parse_args()
    result = benchmark(args.grape.resolve(), args.replicates)
    args.out.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result["summary"], indent=2))


if __name__ == "__main__":
    main()
