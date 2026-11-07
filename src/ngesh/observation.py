"""Observation noise is applied to a copy, never to evolutionary truth."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

from .characters import Evolution, Matrix
from .rng import RandomStream
from .tree import Tree


@dataclass
class ObservationConfig:
    missing: float = 0.0
    coverage_a: float | None = None
    coverage_b: float | None = None
    block_missing: float = 0.0
    selection: float = 1.0
    lumping: float = 0.0
    splitting: float = 0.0
    remove_invariant: bool = False
    remove_singleton: bool = False
    taxon_sampling: float = 1.0

    def validate(self) -> None:
        for name in [
            "missing",
            "block_missing",
            "selection",
            "lumping",
            "splitting",
            "taxon_sampling",
        ]:
            v = getattr(self, name)
            if not math.isfinite(v) or not 0 <= v <= 1:
                raise ValueError(f"{name} must be in [0,1]")
        if (self.coverage_a is None) != (self.coverage_b is None):
            raise ValueError("coverage requires both beta shapes")
        for shape in [self.coverage_a, self.coverage_b]:
            if shape is not None and (not math.isfinite(shape) or shape <= 0):
                raise ValueError("coverage beta shapes must be finite and positive")


@dataclass
class Observed:
    matrix: Matrix
    tree: Tree
    log: list[dict[str, Any]]
    retained_groups: list[int]


def observe(
    evolution: Evolution, config: ObservationConfig | None = None, *, rng: RandomStream
) -> Observed:
    cfg = config or ObservationConfig()
    cfg.validate()
    data = {t: row[:] for t, row in evolution.truth.items()}
    log: list[dict[str, Any]] = []
    n = len(evolution.groups)
    # Selection masks are inherited down the tree, creating nested collections.
    masks: dict[str, set[int]] = {}
    for node in (evolution.history_tree or evolution.tree).traverse():
        parent_mask = masks[node.parent.name] if node.parent else set(range(n))
        if (
            node.parent is None
            or node.attrs.get("sampled_tip")
            or node.parent.attrs.get("sampled_ancestor")
        ):
            mask = set(parent_mask)
        else:
            mask = {
                i for i in sorted(parent_mask) if cfg.selection == 1 or rng.random() < cfg.selection
            }
        masks[node.name] = mask
    for taxon, row in sorted(data.items()):
        coverage = rng.beta(cfg.coverage_a, cfg.coverage_b) if cfg.coverage_a is not None else 1.0
        width = round(cfg.block_missing * n)
        start = rng.index(n - width + 1) if width else 0
        for i, cell in enumerate(row):
            if i not in masks[taxon]:
                row[i] = frozenset()
                log.append({"kind": "selection", "taxon": taxon, "character": i})
            if (
                (cfg.missing and rng.random() < cfg.missing)
                or (coverage < 1 and rng.random() > coverage)
                or (width and start <= i < start + width)
            ):
                row[i] = None
                log.append({"kind": "missing", "taxon": taxon, "character": i})
    for i in range(n):
        alphabet = sorted({s for row in data.values() for s in row[i] or []})
        next_id = max(alphabet, default=0) + 1
        mapping = {}
        for state in alphabet:
            earlier = [s for s in alphabet if s < state]
            mapping[state] = (
                mapping[earlier[rng.index(len(earlier))]]
                if earlier and cfg.lumping and rng.random() < cfg.lumping
                else state
            )
            if mapping[state] != state:
                log.append({"kind": "lumping", "character": i, "old": state, "new": mapping[state]})
        for taxon, row in sorted(data.items()):
            if row[i] is None:
                continue
            mapped = set()
            for state in sorted(row[i]):
                new = mapping[state]
                if cfg.splitting and rng.random() < cfg.splitting:
                    new, next_id = next_id, next_id + 1
                    log.append(
                        {
                            "kind": "splitting",
                            "taxon": taxon,
                            "character": i,
                            "old": state,
                            "new": new,
                        }
                    )
                mapped.add(new)
            row[i] = frozenset(mapped)
    if cfg.taxon_sampling < 1:
        keep = {t for t in sorted(data) if rng.random() < cfg.taxon_sampling}
        if not keep:
            raise ValueError("taxon subsampling removed every observation")
        for t in sorted(set(data) - keep):
            log.append({"kind": "taxon_sampling", "taxon": t})
        data = {t: row for t, row in data.items() if t in keep}
        tree = evolution.tree.prune(keep)
    else:
        tree = evolution.tree
    retained = []
    for i in range(n):
        frequency: dict[int, int] = {}
        cells = {row[i] for row in data.values() if row[i] is not None}
        for row in data.values():
            cell = row[i]
            if cell is None:
                continue
            states = (
                ({1} if cell else {0})
                if evolution.config.model in {"binary", "dollo"}
                else (cell or {-1})
            )
            for state in states:
                frequency[state] = frequency.get(state, 0) + 1
        invariant = len(cells) <= 1
        singleton = sum(v >= 2 for v in frequency.values()) < 2
        if (cfg.remove_invariant and invariant) or (cfg.remove_singleton and singleton):
            log.append(
                {
                    "kind": "ascertainment",
                    "character": i,
                    "invariant": invariant,
                    "singleton": singleton,
                }
            )
        else:
            retained.append(i)
    # Preserve unfiltered truth separately; only observed columns are filtered.
    data = {t: [row[i] for i in retained] for t, row in data.items()}
    return Observed(data, tree, log, retained)
