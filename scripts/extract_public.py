"""Source-specific calibration extraction; only derived summaries are outputs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import statistics

try:
    from .calibrate import compact_summary, file_hash
except ImportError:
    from calibrate import compact_summary, file_hash
from ngesh import Tree
from ngesh.bundle import json_text
from ngesh.matrices import read_nexus
from ngesh.stats import consistency_index


def distribution(values: list[float]) -> dict:
    ordered = sorted(values)
    return {
        "n": len(values),
        "mean": statistics.mean(values) if values else None,
        "variance": statistics.pvariance(values) if values else None,
        "quantiles": [
            ordered[round(q * (len(ordered) - 1))] for q in [0, 0.1, 0.25, 0.5, 0.75, 0.9, 1]
        ]
        if values
        else [],
    }


def yeast_matrix(root: Path) -> tuple[dict, dict, Tree]:
    rows = list(csv.DictReader((root / "species_names.txt").open(), delimiter="\t"))
    allowed = {r["old_speceis_names"] for r in rows if r["Major clade"].lower() != "outgroup"}
    if len(rows) != 343 or len(allowed) != 332:
        raise ValueError("source species roster changed; inspect before extracting")
    prefixes = {}
    for line in (root / "orthomcl_output/orthomcl_SeqIDs_index.txt").open():
        key, separator, value = line.rstrip().partition(": ")
        code = key.partition("_")[0]
        species = value.partition("@")[0]
        if not separator or (code in prefixes and prefixes[code] != species):
            raise ValueError("inconsistent sequence-to-species mapping")
        prefixes[code] = species
    if set(prefixes.values()) != {r["old_speceis_names"] for r in rows}:
        raise ValueError("sequence index differs from published species roster")
    groups = []
    omitted = 0
    for line in (root / "orthomcl_output/orthomcl_clusters.txt").open():
        _, separator, genes = line.partition(":")
        if not separator:
            raise ValueError("unrecognized cluster record")
        present = {prefixes[g.partition("_")[0]] for g in genes.split()} & allowed
        if present:
            groups.append(present)
        else:
            omitted += 1
    one, zero = frozenset({1}), frozenset()
    matrix = {t: [one if t in g else zero for g in groups] for t in sorted(allowed)}
    tree = Tree.from_newick(next(root.rglob("343_2408OGs.treefile")).read_text())
    if {n.name for n in tree.tips()} != set(prefixes.values()):
        raise ValueError("reference tree differs from species index")
    return (
        matrix,
        {
            "excluded_outgroups": 11,
            "outgroup_only_groups_removed": omitted,
            "conditioning": "at least one ingroup member in a reported OrthoMCL cluster; unclustered genes unavailable",
        },
        tree.prune(allowed),
    )


def extract_yeast(root: Path, out: Path) -> None:
    matrix, metadata, tree = yeast_matrix(root)
    n = len(next(iter(matrix.values())))
    out.mkdir(parents=True, exist_ok=True)
    sources = [
        root / "orthomcl_output/orthomcl_clusters.txt",
        root / "orthomcl_output/orthomcl_SeqIDs_index.txt",
        root / "species_names.txt",
        next(root.rglob("343_2408OGs.treefile")),
    ]
    # A predeclared fifth of orthogroups is held out. Taxa are shared, so this is
    # validation on characters, not an independent biological-family validation.
    for split in ["train", "heldout", "all"]:
        indices = [i for i in range(n) if split == "all" or ((i % 5 == 0) == (split == "heldout"))]
        subset = {t: [row[i] for i in indices] for t, row in matrix.items()}
        summary = compact_summary(subset, binary=True)
        prevalence = [
            sum(bool(row[i]) for row in subset.values()) / len(subset) for i in range(len(indices))
        ]
        summary["prevalence_distribution"] = distribution(prevalence)
        summary["trait_singleton_fraction"] = sum(p == 1 / len(subset) for p in prevalence) / len(
            prevalence
        )
        selected = list(range(0, len(indices), max(1, len(indices) // 512)))[:512]
        summary["reference_tree_ci_sample"] = consistency_index(
            {t: [row[i] for i in selected] for t, row in subset.items()}, tree, binary=True
        )
        summary["reference_tree_ci_sample"]["characters_used"] = len(selected)
        result = {
            "schema": "ngesh-target-v2",
            "citation": "Shen et al. 2018; doi:10.1016/j.cell.2018.10.023; figshare 5854692",
            "license": "CC-BY-4.0",
            "inputs": {p.name: file_hash(p) for p in sources},
            "source": "4_orthomcl_output.zip and 1_molecular_phylogenetics.zip",
            "split": split,
            "split_rule": "cluster index modulo 5: zero held out, other four training",
            "extraction": metadata,
            "summary": summary,
        }
        (out / f"yeast-{split}.json").write_text(json_text(result), encoding="utf-8")


def extract_morphology(root: Path, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    report = []
    for p in sorted(root.rglob("*overlap.nex")):
        name = p.parts[-4]
        source = p
        try:
            if p.read_bytes().count(b"\0") == p.stat().st_size:
                fallback = p.parents[1] / "TNT/nexus_xinact.nex"
                if not fallback.exists() or b"matrix" not in fallback.read_bytes().lower():
                    raise ValueError(
                        "archive member is entirely NUL bytes; no valid morphological fallback"
                    )
                source = fallback
            matrix = read_nexus(source)
            # Mk has no applicability state. Both unknown and inapplicable cells
            # are unobserved for this calibration; the recoding is explicit.
            gaps = sum(c == frozenset() for row in matrix.values() for c in row)
            matrix = {t: [c if c else None for c in row] for t, row in matrix.items()}
            summary = compact_summary(matrix)
            summary["observed_state_count_histogram"] = summary["class_histogram"]
            result = {
                "schema": "ngesh-target-v2",
                "source": str(source.relative_to(root)),
                "citation": "Keating, Garwood & Sansom 2023; doi:10.1186/s12862-023-02131-z; Zenodo 6579584",
                "license": "CC-BY-4.0",
                "inputs": {str(source.relative_to(root)): file_hash(source)},
                "extraction": {
                    "gap_policy": "inapplicable and unknown both unobserved for Mk calibration",
                    "gaps_recoded": gaps,
                    "fallback_from_corrupt_member": source != p,
                },
                "summary": summary,
            }
            (out / f"{name}.json").write_text(json_text(result), encoding="utf-8")
            report.append(
                {"dataset": name, "status": "extracted", "source": str(source.relative_to(root))}
            )
        except ValueError as error:
            report.append(
                {
                    "dataset": name,
                    "status": "excluded",
                    "reason": str(error),
                    "member": str(p.relative_to(root)),
                    "sha256": file_hash(p),
                }
            )
    (out.parent / "morphology-extraction.json").write_text(json_text(report), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", choices=["yeast", "morphology"])
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.source == "yeast":
        extract_yeast(args.root, args.out)
    else:
        extract_morphology(args.root, args.out)


if __name__ == "__main__":
    main()
