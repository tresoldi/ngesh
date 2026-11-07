"""Deterministic, checksummed dataset bundles and strict long-format readers."""

from __future__ import annotations

import csv
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
from typing import Any

from . import __version__
from .characters import Cell, Evolution, Matrix
from .models import ancestor_internal
from .stats import summarize
from .tree import Tree

HEADER = "Language\tConcept\tCognateset"


def json_text(value: Any) -> str:
    return json.dumps(value, sort_keys=True, ensure_ascii=False, indent=2, allow_nan=False) + "\n"


def matrix_json(matrix: Matrix) -> dict[str, list[list[int] | None]]:
    return {
        t: [sorted(c) if c is not None else None for c in row] for t, row in sorted(matrix.items())
    }


def read_matrix(path: str | Path) -> Matrix:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return {t: [frozenset(c) if c is not None else None for c in row] for t, row in data.items()}


def write_long(path: Path, matrix: Matrix, groups: list[str]) -> None:
    lines = [HEADER + "\n"]
    for taxon, row in sorted(matrix.items()):
        if len(row) != len(groups):
            raise ValueError("matrix width differs from groups")
        for group, cell in zip(groups, row):
            for token in (taxon, group):
                if not token or any(c in token for c in "\t\r\n"):
                    raise ValueError("TSV labels must be nonempty and contain no tabs/newlines")
            for state in sorted(cell or []):
                if not isinstance(state, int) or state < 0:
                    raise ValueError("state IDs must be nonnegative integers")
                line = f"{taxon}\t{group}\t{group}.{state}\n"
                if len(line.encode("utf-8")) >= 4096:
                    raise ValueError("long-format line exceeds GRAPE's 4095-byte limit")
                lines.append(line)
    path.write_text("".join(lines), encoding="utf-8")


def read_long(path: str | Path, taxa: list[str], groups: list[str]) -> Matrix:
    """Read sparse rows. No row means unknown here, never assumed binary zero.

    Use observed.json / matrices for the binary absent-vs-missing distinction.
    """
    result: Matrix = {t: [None] * len(groups) for t in taxa}
    indices = {g: i for i, g in enumerate(groups)}
    with Path(path).open(encoding="utf-8", newline="") as handle:
        if handle.readline().rstrip("\r\n") != HEADER:
            raise ValueError("unexpected data.tsv header")
        for line in handle:
            if len(line.encode("utf-8")) >= 4096:
                raise ValueError("overlong TSV line")
            values = line.rstrip("\r\n").split("\t")
            if len(values) != 3:
                raise ValueError("data.tsv requires exactly three columns")
            taxon, group, state = values
            prefix, sep, integer = state.rpartition(".")
            if (
                taxon not in result
                or group not in indices
                or not sep
                or not prefix
                or not integer.isascii()
                or not integer.isdigit()
            ):
                raise ValueError("invalid data.tsv record")
            i = indices[group]
            result[taxon][i] = (result[taxon][i] or frozenset()) | {int(integer)}
    return result


def _quote(text: str) -> str:
    return "'" + text.replace("'", "''") + "'"


def _nexus(
    path: Path,
    matrix: Matrix,
    groups: list[str],
    *,
    binary: bool,
    presence: bool = False,
    state_count: int | None = None,
) -> list[tuple[int, int]]:
    taxa = sorted(matrix)
    columns: list[tuple[int, int]] = []
    if binary:
        for i in range(len(groups)):
            states = (
                [1] if presence else sorted({s for row in matrix.values() for s in row[i] or []})
            )
            columns.extend((i, s) for s in states)
    else:
        columns = [(i, -1) for i in range(len(groups))]

    def token(cell: Cell, state: int) -> str:
        if cell is None:
            return "?"
        if binary:
            return "1" if state in cell else "0"
        if presence:
            return "1" if cell else "0"
        if not cell:
            return "-"
        return (
            str(next(iter(cell)))
            if len(cell) == 1
            else "{" + " ".join(map(str, sorted(cell))) + "}"
        )

    data = [f"{_quote(t)} " + " ".join(token(matrix[t][i], s) for i, s in columns) for t in taxa]
    # TOKENS permits integer states beyond the single-symbol NEXUS alphabet.
    states = sorted({s for row in matrix.values() for c in row for s in c or []})
    single_symbol = not binary and state_count is not None and max(states, default=0) < 10
    symbols = (
        "".join(str(i) for i in range(max(state_count, max(states, default=0) + 1)))
        if single_symbol
        else "01"
        if binary or presence
        else " ".join(map(str, states or [0]))
    )
    fmt = (
        "DATATYPE=STANDARD "
        + ("" if binary or single_symbol else "TOKENS ")
        + f'SYMBOLS="{symbols}" MISSING=? GAP=-;'
    )
    text = [
        "#NEXUS",
        "BEGIN DATA;",
        f"DIMENSIONS NTAX={len(taxa)} NCHAR={len(columns)};",
        "FORMAT " + fmt,
        "MATRIX",
        *data,
        ";",
        "END;",
    ]
    path.write_text("\n".join(text) + "\n", encoding="utf-8")
    return columns


def write_bundle(
    directory: str | Path,
    evolution: Evolution,
    *,
    params: dict[str, Any],
    observed: Matrix | None = None,
    truth_tree: Tree | None = None,
    observation_log: list[dict[str, Any]] | None = None,
    include_ci: bool = True,
    max_quartets: int = 1000,
    retained_groups: list[int] | None = None,
) -> dict[str, Any]:
    """Write files, then params.json last as the completion manifest.

    Files are deterministic; manifest does not hash itself (no recursive hash).
    Callers stage/rename a temporary directory when atomic publication is needed.
    """
    path = Path(directory)
    path.mkdir(parents=True, exist_ok=False)
    data = evolution.truth if observed is None else observed
    tree = truth_tree or evolution.tree
    groups = (
        evolution.groups
        if retained_groups is None
        else [evolution.groups[i] for i in retained_groups]
    )
    binary = evolution.config.model in {"dollo", "binary"}
    (path / "truth.nwk").write_text(tree.newick() + "\n", encoding="utf-8")
    if evolution.history_tree is not None:
        (path / "history_tree.nwk").write_text(
            evolution.history_tree.newick() + "\n", encoding="utf-8"
        )
    if params.get("tree_config", {}).get("ancestor_rate", 0) or any(
        n.attrs.get("sampled_ancestor") for n in tree.traverse()
    ):
        (path / "truth_sa_leaves.nwk").write_text(tree.newick() + "\n", encoding="utf-8")
        (path / "truth_sa_internal.nwk").write_text(
            ancestor_internal(tree).newick() + "\n", encoding="utf-8"
        )
    lineage_tree = evolution.history_tree or evolution.tree
    (path / "lineages.json").write_text(
        json_text(
            [
                {
                    "id": n.name,
                    "parent": n.parent.name if n.parent else None,
                    "start": n.parent.time if n.parent else 0.0,
                    "end": n.time,
                    "extinct": n.attrs.get("extinct", False),
                    "sampled_tip": n.attrs.get("sampled_tip", False),
                    "sampled_ancestor": n.attrs.get("sampled_ancestor", False),
                }
                for n in lineage_tree.traverse()
            ]
        ),
        encoding="utf-8",
    )
    write_long(path / "data.tsv", data, groups)
    (path / "truth.json").write_text(json_text(matrix_json(evolution.truth)), encoding="utf-8")
    (path / "observed.json").write_text(json_text(matrix_json(data)), encoding="utf-8")
    (path / "observations.json").write_text(json_text(observation_log or []), encoding="utf-8")
    columns = _nexus(path / "data.bin.nex", data, groups, binary=True, presence=binary)
    _nexus(
        path / "data.multi.nex",
        data,
        groups,
        binary=False,
        presence=binary,
        state_count=evolution.config.k if evolution.config.model == "mk" else None,
    )
    if evolution.config.model == "mk" and groups:
        with (path / "data.multi.nex").open("a", encoding="utf-8") as handle:
            handle.write(f"BEGIN SETS;\nCHARSET morphology = 1-{len(groups)};\nEND;\n")
            if evolution.config.ordered:
                handle.write(f"BEGIN ASSUMPTIONS;\nCTYPE ordered: 1-{len(groups)};\nEND;\n")
    # PHYLIP relaxed names are whitespace-free; map names to stable IDs explicitly.
    aliases = {t: f"T{i:06d}" for i, t in enumerate(sorted(data), 1)}
    rows = []
    for taxon in sorted(data):
        bits = "".join(
            "?" if data[taxon][i] is None else ("1" if s in data[taxon][i] else "0")
            for i, s in columns
        )
        rows.append(f"{aliases[taxon]} {bits}")
    (path / "data.bin.phy").write_text(
        f"{len(data)} {len(columns)}\n" + "\n".join(rows) + "\n", encoding="utf-8"
    )
    metadata = {
        "groups": groups,
        "taxa": sorted(data),
        "phylip_taxa": aliases,
        "binary_columns": [{"group": groups[i], "state": s} for i, s in columns],
        "model": evolution.config.model,
        "morphology_partition": list(range(1, len(groups) + 1))
        if evolution.config.model == "mk"
        else [],
        "ordered": getattr(evolution.config, "ordered", False),
        "retained_group_indices": retained_groups,
        "encoding": "TSV omits absent and missing cells; observed.json distinguishes [] and null. Binary present=1.",
    }
    (path / "columns.json").write_text(json_text(metadata), encoding="utf-8")
    with (path / "history.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["time", "kind", "donor", "recipient", "character", "old", "new"])
        for e in evolution.events:
            writer.writerow(
                [
                    format(e.time, ".17g"),
                    e.kind,
                    e.donor,
                    e.recipient,
                    evolution.groups[e.character],
                    json.dumps(e.old),
                    json.dumps(e.new),
                ]
            )
    horizontal_events = [e for e in evolution.events if e.donor]
    h = params.get("horizontal_config", {})
    if horizontal_events or h.get("borrowing_rate", 0) or h.get("contamination_fraction", 0):
        from .horizontal import displayed_network

        (path / "network.enwk").write_text(
            displayed_network(evolution.history_tree or evolution.tree, horizontal_events).newick()
            + "\n",
            encoding="utf-8",
        )
        with (path / "events.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["time", "kind", "donor", "recipient", "character", "old", "new"])
            for e in horizontal_events:
                writer.writerow(
                    [
                        format(e.time, ".17g"),
                        e.kind,
                        e.donor,
                        e.recipient,
                        evolution.groups[e.character],
                        json.dumps(e.old),
                        json.dumps(e.new),
                    ]
                )
    summary = summarize(data, tree, binary=binary, max_quartets=max_quartets, include_ci=include_ci)
    (path / "summary.json").write_text(json_text(summary), encoding="utf-8")
    manifest = dict(
        params,
        version=__version__,
        character_config=asdict(evolution.config),
        rates=evolution.rates,
        branch_rates=evolution.branch_rates,
        checksums={
            f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(path.iterdir())
        },
    )
    (path / "params.json").write_text(json_text(manifest), encoding="utf-8")
    return manifest


def verify_bundle(directory: str | Path) -> bool:
    try:
        path = Path(directory)
        manifest = json.loads((path / "params.json").read_text(encoding="utf-8"))
        files = manifest["checksums"]
        return bool(files) and all(
            Path(name).name == name
            and hashlib.sha256((path / name).read_bytes()).hexdigest() == digest
            for name, digest in files.items()
        )
    except (OSError, ValueError, KeyError, TypeError):
        return False
