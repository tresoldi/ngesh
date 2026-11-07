"""Extract derived targets, fetch lexical CLDF sources, and fit preset grids.

Third-party inputs belong outside version control. This script commits nothing.
See docs/calibration.md for citations, licenses, and source-specific caveats.
"""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
from pathlib import Path
import statistics
import subprocess
from typing import Any

from ngesh.bundle import json_text
from ngesh.characters import Matrix
from ngesh.matrices import read_nexus
from ngesh.simulation import simulate
from ngesh.stats import summarize

LEXICAL = [
    "dravlex",
    "walworthpolynesian",
    "kitchensemitic",
    "savelyevturkic",
    "syrjaenenuralic",
    "utoaztecan",
    "chaconarawakan",
    "liusinitic",
    "wichmannmixezoquean",
    "constenlachibchan",
    "sidwellbahnaric",
    "leekoreanic",
    "mannburmish",
    "dunnaslian",
    "iecor",
    "tuled",
]


def compact_summary(matrix: Matrix, *, binary: bool = False) -> dict[str, Any]:
    summary = summarize(matrix, binary=binary, include_ci=False)
    distances = [
        v["distance"] for v in summary.pop("pairwise_distances") if v["distance"] is not None
    ]
    counts = summary.pop("classes_per_group")
    coverage = list(summary.pop("coverage").values())

    def distribution(values: list[float]) -> dict[str, Any]:
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

    ncols = len(next(iter(matrix.values()))) if matrix else 0
    constant = 0
    for i in range(ncols):
        cells = {row[i] for row in matrix.values() if row[i] is not None}
        constant += bool(cells) and len(cells) == 1
    summary["constant_fraction"] = constant / ncols if ncols else 0.0
    if binary:
        prevalence = [
            sum(bool(row[i]) for row in matrix.values()) / len(matrix) for i in range(ncols)
        ]
        summary["prevalence_distribution"] = distribution(prevalence)
        summary["trait_singleton_fraction"] = (
            sum(p == 1 / len(matrix) for p in prevalence) / ncols if ncols else 0
        )
    summary["distance_distribution"] = distribution(distances)
    summary["coverage_distribution"] = distribution(coverage)
    summary["class_distribution"] = distribution(counts)
    summary["class_histogram"] = {str(k): v for k, v in sorted(Counter(counts).items())}
    return summary


def cldf_matrix(directory: Path) -> Matrix:
    forms_path = directory / "forms.csv"
    cognates_path = directory / "cognates.csv"
    with forms_path.open(encoding="utf-8-sig", newline="") as f:
        forms = {row["ID"]: row for row in csv.DictReader(f)}
    cognates: dict[str, list[str]] = {}
    if cognates_path.exists():
        with cognates_path.open(encoding="utf-8-sig", newline="") as f:
            for row in csv.DictReader(f):
                value = row.get("Cognateset_ID", "")
                if value:
                    cognates.setdefault(row["Form_ID"], []).append(value)
    else:
        for fid, row in forms.items():
            value = row.get("Cognateset_ID", row.get("Cognacy", ""))
            if value:
                cognates[fid] = value.split()
    if not cognates:
        raise ValueError(
            "no cognate judgments found; unjudged forms are missing, not unique cognates"
        )
    groups = sorted({forms[fid]["Parameter_ID"] for fid in cognates if fid in forms})
    taxa = sorted({forms[fid]["Language_ID"] for fid in cognates if fid in forms})
    indices = {g: i for i, g in enumerate(groups)}
    result: Matrix = {t: [None] * len(groups) for t in taxa}
    labels: dict[str, dict[str, int]] = {g: {} for g in groups}
    for fid in sorted(cognates):
        if fid not in forms:
            raise ValueError(f"cognate refers to unknown form {fid}")
        form = forms[fid]
        taxon, group = form["Language_ID"], form["Parameter_ID"]
        i = indices[group]
        for label in cognates[fid]:
            state = labels[group].setdefault(label, len(labels[group]))
            result[taxon][i] = (result[taxon][i] or frozenset()) | {state}
    return result


def orthogroup_matrix(path: Path, separator: str) -> Matrix:
    """OrthoMCL groups: group: species|gene ... . Explicit separator required."""
    groups = []
    taxa = set()
    with path.open(encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            _, sep, genes = line.partition(":")
            if not sep:
                raise ValueError("expected OrthoMCL group: species|gene ...")
            present = set()
            for gene in genes.split():
                species, sep, _ = gene.partition(separator)
                if not sep:
                    raise ValueError(
                        "gene identifiers lack specified species separator; map source identifiers first"
                    )
                present.add(species)
            taxa.update(present)
            groups.append(present)
    one, zero = frozenset({1}), frozenset()
    return {t: [one if t in g else zero for g in groups] for t in sorted(taxa)}


def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def extract(
    source: Path,
    kind: str,
    citation: str,
    license_text: str,
    separator: str = "|",
    characters: list[int] | None = None,
) -> dict[str, Any]:
    if kind == "cldf":
        matrix = cldf_matrix(source)
        files = [f for f in source.glob("*.csv") if f.name in {"forms.csv", "cognates.csv"}]
    elif kind == "orthogroups":
        matrix = orthogroup_matrix(source, separator)
        files = [source]
    else:
        matrix = read_nexus(source, presence=kind == "binary-nexus", characters=characters)
        files = [source]
    return {
        "schema": "ngesh-target-v1",
        "source": str(source.name),
        "citation": citation,
        "license": license_text,
        "inputs": {f.name: file_hash(f) for f in files},
        "summary": compact_summary(matrix, binary=kind in {"orthogroups", "binary-nexus"}),
    }


def fetch_lexical(work: Path, out: Path) -> list[dict[str, Any]]:
    out.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)

    def fetch(name: str) -> dict[str, Any]:
        repo = ("tupian-language-resources/" if name == "tuled" else "lexibank/") + name
        dest = work / name
        try:
            if not dest.exists():
                subprocess.run(
                    [
                        "git",
                        "clone",
                        "--depth",
                        "1",
                        "--filter=blob:none",
                        "--sparse",
                        f"https://github.com/{repo}.git",
                        str(dest),
                    ],
                    check=True,
                    capture_output=True,
                )
                subprocess.run(
                    ["git", "-C", str(dest), "sparse-checkout", "set", "cldf"],
                    check=True,
                    capture_output=True,
                )
            revision = subprocess.check_output(
                ["git", "-C", str(dest), "rev-parse", "HEAD"], text=True
            ).strip()
            license_files = list(dest.glob("LICENSE*")) + list(dest.glob("COPYING*"))
            license_text = "\n".join(
                f.read_text(encoding="utf-8", errors="replace") for f in license_files
            )
            metadata = dest / "metadata.json"
            if metadata.exists():
                meta = json.loads(metadata.read_text())
                license_text = meta.get("license", license_text)
            result = extract(
                dest / "cldf",
                "cldf",
                f"https://github.com/{repo}/tree/{revision}",
                license_text or "See repository; no top-level license file located",
            )
            if metadata.exists():
                result["bibliographic_citation"] = meta.get("citation", meta.get("title", ""))
                result["source_url"] = meta.get("url", "")
            result["repository"] = repo
            result["revision"] = revision
            (out / (name + ".json")).write_text(json_text(result), encoding="utf-8")
            return {"name": name, "status": "ok"}
        except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as error:
            return {"name": name, "status": "failed", "error": str(error)}

    with ThreadPoolExecutor(max_workers=4) as executor:
        report = list(executor.map(fetch, LEXICAL))
    return report


def features(summary: dict[str, Any]) -> dict[str, float]:
    return {
        "distance": summary["distance_distribution"]["mean"],
        "singleton": summary["singleton_fraction"],
        "missing": summary["missing_fraction"],
        "classes_per_taxon": summary["class_distribution"]["mean"] / max(1, summary["taxa"]),
    }


def fit_legacy(targets: Path, out: Path, preset: str, replicates: int = 3) -> None:
    import itertools

    sources = [json.loads(p.read_text()) for p in sorted(targets.glob("*.json"))]
    if not sources:
        raise ValueError("no target summaries; extraction is required before fitting")
    target_features = [features(s["summary"]) for s in sources]
    target = {k: statistics.mean(f[k] for f in target_features) for k in target_features[0]}
    trials = []
    for rate, alpha, missing in itertools.product(
        [0.03, 0.1, 0.3, 0.6], [0.3, 0.7, 2.0], [0, 0.15, 0.4]
    ):
        config = {
            "tree": {"n_taxa": 50},
            "characters": {"n_characters": 200, "rate": rate, "alpha": alpha},
            "observation": {"missing": missing},
        }
        draws = []
        for replicate in range(replicates):
            # Include observation rather than fitting complete latent data to missing targets.
            from ngesh.observation import ObservationConfig, observe
            from ngesh.rng import RandomStream, derive_seed

            result = simulate(
                preset=preset,
                config=config,
                seed=derive_seed("calibration-v1", "replicate", replicate),
            )
            obs = observe(
                result,
                ObservationConfig(missing=missing),
                rng=RandomStream(derive_seed("calibration-v1", "missing", replicate)),
            )
            draws.append(features(compact_summary(obs.matrix)))
        simulated = {k: statistics.mean(d[k] for d in draws) for k in target}
        # Equal weights for bounded fractions; no hidden normalization tuning.
        score = sum((simulated[k] - target[k]) ** 2 for k in target)
        trials.append({"score": score, "config": config, "simulated_features": simulated})
    trials.sort(key=lambda t: t["score"])
    result = {
        "schema": "ngesh-fit-v1",
        "preset": preset,
        "target_files": [p.name for p in sorted(targets.glob("*.json"))],
        "target_features": target,
        "replicates_per_point": replicates,
        "procedure": "36-point grid, equal family weights, sum squared fractional features",
        "best": trials[0],
        "trials": trials,
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json_text(result), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    fetch = sub.add_parser("fetch-lexical")
    fetch.add_argument("--work", type=Path, required=True)
    fetch.add_argument("--out", type=Path, required=True)
    extraction = sub.add_parser("extract")
    extraction.add_argument("kind", choices=["cldf", "nexus", "binary-nexus", "orthogroups"])
    extraction.add_argument("source", type=Path)
    extraction.add_argument("--out", type=Path, required=True)
    extraction.add_argument("--citation", required=True)
    extraction.add_argument("--license", required=True)
    extraction.add_argument("--separator", default="|")
    extraction.add_argument("--characters", help="1-based inclusive morphology range, e.g. 1:152")
    fitting = sub.add_parser("fit")
    fitting.add_argument("--targets", type=Path, required=True)
    fitting.add_argument("--out", type=Path, required=True)
    fitting.add_argument(
        "--preset", choices=["lexical", "stemma", "genecontent", "morphology"], required=True
    )
    fitting.add_argument("--replicates", type=int, default=1)
    fitting.add_argument(
        "--legacy",
        action="store_true",
        help="reproduce the alpha3 preliminary grid without holdout",
    )
    args = parser.parse_args()
    if args.command == "fetch-lexical":
        print(json.dumps(fetch_lexical(args.work, args.out), indent=2))
    elif args.command == "fit":
        if args.legacy:
            fit_legacy(args.targets, args.out, args.preset, args.replicates)
        else:
            try:
                from .fit_presets import fit
            except ImportError:
                from fit_presets import fit
            fit(args.preset, args.targets, args.out, args.replicates)
    else:
        characters = None
        if args.characters:
            a, b = map(int, args.characters.split(":"))
            characters = list(range(a - 1, b))
        result = extract(
            args.source, args.kind, args.citation, args.license, args.separator, characters
        )
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json_text(result), encoding="utf-8")


if __name__ == "__main__":
    main()
