"""Optional external-consumer integration pilot; GRAPE is never a runtime dependency.

Build GRAPE separately and pass its executable. This checks serialization, taxon
rosters, event lifetimes and repeated consumer output; it does not score trees or
select simulator parameters for any inference method.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

from ngesh import Tree, simulate
from ngesh.bundle import json_text, read_long, read_matrix, verify_bundle


def pilot(executable: Path, out: Path, revision: str) -> dict:
    out.mkdir(parents=True, exist_ok=False)
    cases = []
    for preset in ["lexical", "stemma", "genecontent", "morphology"]:
        for seed in [11, 23, 47]:
            cases.append((f"{preset}-{seed}", preset, seed, {}))
    cases.extend(
        [
            (
                "lexical-transfer",
                "lexical",
                59,
                {
                    "horizontal": {"borrowing_rate": 0.05, "distance": "geography"},
                    "characters": {"synonym_duration": 0.3, "clock_sigma": 0.7},
                },
            ),
            (
                "stemma-contamination",
                "stemma",
                61,
                {"horizontal": {"contamination_fraction": 0.15, "contamination_blocks": True}},
            ),
            (
                "stemma-selection",
                "stemma",
                67,
                {"observation": {"selection": 0.85, "block_missing": 0.2}},
            ),
            (
                "gene-transfer",
                "genecontent",
                71,
                {"horizontal": {"borrowing_rate": 0.02, "mode": "add"}},
            ),
            ("unrepresented-taxa", "lexical", 73, {"observation": {"missing": 1.0}}),
        ]
    )
    records = []
    for name, preset, seed, extra in cases:
        config = {"tree": {"n_taxa": 20}, "characters": {"n_characters": 200}, **extra}
        # Merge dimensions into stress sections rather than overwriting them.
        config.setdefault("tree", {}).setdefault("n_taxa", 20)
        config.setdefault("characters", {}).setdefault("n_characters", 200)
        directory = out / name
        start = time.perf_counter()
        simulate(seed=seed, preset=preset, config=config, out=directory)
        assert verify_bundle(directory)
        columns = json.loads((directory / "columns.json").read_text())
        observed = read_matrix(directory / "observed.json")
        parsed = read_long(directory / "data.tsv", columns["taxa"], columns["groups"])
        assert all(
            parsed[t][i] == cell
            for t, row in observed.items()
            for i, cell in enumerate(row)
            if cell
        )
        represented = {t for t, row in parsed.items() if any(c is not None for c in row)}
        lifetimes = {r["id"]: r for r in json.loads((directory / "lineages.json").read_text())}
        import csv

        transfers = 0
        if (directory / "events.tsv").exists():
            for e in csv.DictReader((directory / "events.tsv").open(), delimiter="\t"):
                donor = lifetimes[e["donor"]]
                assert donor["start"] <= float(e["time"]) < donor["end"]
                transfers += 1
            Tree.from_newick((directory / "network.enwk").read_text())
        runs = []
        for _ in range(2):
            command = [
                str(executable),
                "--graph",
                "adjusted",
                "--community",
                "louvain",
                "--seed",
                "42",
                "--missing-data",
                "ignore",
                str(directory / "data.tsv"),
            ]
            run = subprocess.run(command, text=True, capture_output=True, timeout=120)
            runs.append(run)
        expected_empty = not represented
        if expected_empty:
            assert all(r.returncode != 0 for r in runs), "empty TSV must fail visibly"
            inferred = set()
        else:
            assert all(r.returncode == 0 for r in runs), runs[0].stderr
            assert runs[0].stdout == runs[1].stdout, "consumer must reproduce fixed-seed tree"
            inferred = {n.name for n in Tree.from_newick(runs[0].stdout.strip()).tips()}
            assert inferred == represented, (name, inferred ^ represented)
        records.append(
            {
                "case": name,
                "preset": preset,
                "seed": seed,
                "config": config,
                "observed_taxa": len(observed),
                "tsv_taxa": len(represented),
                "unrepresented_taxa": sorted(set(observed) - represented),
                "inferred_taxa": len(inferred),
                "transfers": transfers,
                "expected_empty_rejection": expected_empty,
                "consumer_exit": runs[0].returncode,
                "tree_sha256": hashlib.sha256(runs[0].stdout.encode()).hexdigest()
                if inferred
                else None,
                "elapsed_seconds": time.perf_counter() - start,
            }
        )
    report = {
        "schema": "ngesh-consumer-pilot-v1",
        "consumer": "GRAPE native C99",
        "consumer_revision": revision,
        "consumer_tests": "ctest: 2/2 passed",
        "purpose": "integration only; no inference accuracy or model tuning claim",
        "cases": records,
    }
    (out / "report.json").write_text(json_text(report), encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--revision", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    report = pilot(args.executable.resolve(), args.out, args.revision)
    args.report.write_text(json_text(report), encoding="utf-8")
    print(json.dumps({"cases": len(report["cases"]), "passed": True}))


if __name__ == "__main__":
    main()
