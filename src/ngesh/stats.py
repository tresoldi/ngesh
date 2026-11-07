"""Descriptive benchmark statistics; no inference or external dependencies."""

from __future__ import annotations

from itertools import combinations
import math
import statistics
from typing import Any

from .characters import Matrix
from .tree import Tree


def pairwise_distances(matrix: Matrix) -> dict[tuple[str, str], float | None]:
    """Mean per-group Jaccard distance across jointly observed cells.

    Empty/empty contributes zero. In binary datasets this is Hamming distance
    across observed traits. Missing cells are pairwise deleted.
    """
    result = {}
    binary_states = {frozenset(), frozenset({1})}
    if all(cell is None or cell in binary_states for row in matrix.values() for cell in row):
        # Large orthogroup matrices otherwise require billions of Python cell
        # comparisons. Bitsets give the exact same Hamming distance in C.
        encoded = {}
        for taxon, row in matrix.items():
            present = bytearray((len(row) + 7) // 8)
            known = bytearray(len(present))
            for i, cell in enumerate(row):
                if cell is not None:
                    known[i // 8] |= 1 << (i % 8)
                    if cell:
                        present[i // 8] |= 1 << (i % 8)
            encoded[taxon] = (int.from_bytes(present, "little"), int.from_bytes(known, "little"))
        for a, b in combinations(sorted(matrix), 2):
            x, kx = encoded[a]
            y, ky = encoded[b]
            joint = kx & ky
            result[a, b] = ((x ^ y) & joint).bit_count() / joint.bit_count() if joint else None
        return result
    for a, b in combinations(sorted(matrix), 2):
        values = []
        for x, y in zip(matrix[a], matrix[b]):
            if x is None or y is None:
                continue
            union = x | y
            values.append(1 - len(x & y) / len(union) if union else 0.0)
        result[a, b] = statistics.mean(values) if values else None
    return result


def quartets(
    distances: dict[tuple[str, str], float | None], max_quartets: int = 1000
) -> dict[str, Any]:
    """Holland delta and Q-residual on deterministically spaced quartets.

    Q is the squared four-point residual divided by squared mean distance.
    All quartets are used if there are <=max_quartets. No RNG is consumed.
    """
    taxa = sorted({t for pair in distances for t in pair})
    total = math.comb(len(taxa), 4) if len(taxa) >= 4 else 0
    delta, residual = [], []
    target = 0
    next_index = 0
    limit = min(total, max_quartets)
    if limit < 0:
        raise ValueError("max_quartets must be nonnegative")

    # Unrank only the selected combinations; avoid iterating O(n^4) quartets.
    def unrank(rank: int) -> list[str]:
        chosen, start = [], 0
        for remaining in range(4, 0, -1):
            for j in range(start, len(taxa)):
                count = math.comb(len(taxa) - j - 1, remaining - 1)
                if rank < count:
                    chosen.append(taxa[j])
                    start = j + 1
                    break
                rank -= count
        return chosen

    while target < limit:
        a, b, c, d = unrank(next_index)
        ds = [
            distances.get(tuple(sorted(pair)))
            for pair in [(a, b), (c, d), (a, c), (b, d), (a, d), (b, c)]
        ]
        if all(x is not None for x in ds):
            lo, mid, hi = sorted([ds[0] + ds[1], ds[2] + ds[3], ds[4] + ds[5]])
            delta.append((hi - mid) / (hi - lo) if hi > lo else 0.0)
            residual.append((hi - mid) ** 2)
        target += 1
        next_index = target * total // limit if limit else 0
    mean = (
        statistics.mean(x for x in distances.values() if x is not None)
        if any(x is not None for x in distances.values())
        else 0
    )
    return {
        "delta": statistics.mean(delta) if delta else None,
        "q_residual": statistics.mean(residual) / mean**2
        if residual and mean
        else (0.0 if residual else None),
        "quartets_used": len(delta),
        "quartets_total": total,
        "quartet_sampling": "evenly_spaced_lexicographic",
    }


def consistency_index(matrix: Matrix, tree: Tree, *, binary: bool = False) -> dict[str, Any]:
    """Unordered Fitch CI, treating polymorphisms as ambiguous state sets.

    Missing cells admit any state; binary empty cells are observed zero.
    Lexical/textual empty cells admit absence as a separate state (-1).
    CI is undefined when no informative changes exist.
    """
    if not matrix:
        return {"ci": None, "minimum_steps": 0, "tree_steps": 0}
    if not binary and any(c is not None and len(c) > 1 for row in matrix.values() for c in row):
        columns = [
            (i, s)
            for i in range(len(next(iter(matrix.values()))))
            for s in sorted({s for row in matrix.values() for s in row[i] or []})
        ]
        recoded = {
            t: [
                None if row[i] is None else frozenset({1}) if s in row[i] else frozenset()
                for i, s in columns
            ]
            for t, row in matrix.items()
        }
        result = consistency_index(recoded, tree, binary=True)
        result["ci_encoding"] = "binary_membership_for_polymorphism"
        return result
    ncols = len(next(iter(matrix.values())))
    minimum = steps = 0
    nodes = list(tree.traverse("postorder"))
    for c in range(ncols):
        observed = set()
        for row in matrix.values():
            if row[c] is not None:
                observed.update(({1} if row[c] else {0}) if binary else (row[c] or {-1}))
        if not observed:
            continue
        minimum += len(observed) - 1
        # Sankoff handles multifurcations correctly (pairwise Fitch does not).
        alphabet = sorted(observed)
        costs = {}
        for node in nodes:
            own = matrix.get(node.name)
            allowed = (
                None
                if own is None or own[c] is None
                else (({1} if own[c] else {0}) if binary else (own[c] or {-1}))
            )
            vector = [0 if allowed is None or s in allowed else math.inf for s in alphabet]
            for child in node.children:
                child_cost = costs[child]
                best = min(child_cost)
                vector = [v + min(child_cost[i], best + 1) for i, v in enumerate(vector)]
            costs[node] = vector
        steps += min(costs[tree.root])
    return {"ci": minimum / steps if steps else None, "minimum_steps": minimum, "tree_steps": steps}


def summarize(
    matrix: Matrix,
    tree: Tree | None = None,
    *,
    binary: bool = False,
    max_quartets: int = 1000,
    include_ci: bool = True,
) -> dict[str, Any]:
    taxa = sorted(matrix)
    ncols = len(matrix[taxa[0]]) if taxa else 0
    counts, singleton = [], []
    for c in range(ncols):
        frequencies: dict[int, int] = {}
        for row in matrix.values():
            for state in row[c] or []:
                frequencies[state] = frequencies.get(state, 0) + 1
        counts.append(len(frequencies))
        singleton.extend(v == 1 for v in frequencies.values())
    missing = sum(cell is None for row in matrix.values() for cell in row)
    distances = pairwise_distances(matrix)
    result = {
        "taxa": len(taxa),
        "characters": ncols,
        "classes_per_group": counts,
        "singleton_fraction": statistics.mean(singleton) if singleton else 0,
        "missing_fraction": missing / (len(taxa) * ncols) if taxa and ncols else 0,
        "coverage": {
            t: sum(c is not None for c in matrix[t]) / ncols if ncols else 0 for t in taxa
        },
        "pairwise_distances": [{"a": a, "b": b, "distance": d} for (a, b), d in distances.items()],
        **quartets(distances, max_quartets),
    }
    if tree and include_ci:
        result.update(consistency_index(matrix, tree, binary=binary))
    return result
