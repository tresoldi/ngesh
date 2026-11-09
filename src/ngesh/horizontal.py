"""Horizontal-process parameters and displayed primary-tree networks."""

from __future__ import annotations

from dataclasses import dataclass
import math

from .characters import Event
from .rng import positive
from .tree import Node, Tree


@dataclass
class HorizontalConfig:
    borrowing_rate: float = 0.0
    scale: float = 1.0
    distance: str = "patristic"
    diffusion: float = 1.0
    mode: str = "replace"
    contamination_fraction: float = 0.0
    contamination_blocks: bool = False
    # Contact episodes: sustained exchange between two contemporaneous lineages.
    episode_rate: float = 0.0  # onsets per lineage per unit time
    episode_duration: float = 0.2  # mean duration (exponential), time units
    episode_intensity: float = 0.5  # transfers per character per unit time during an episode
    episode_symmetric: bool = True  # each transfer goes either way, else partner -> initiator

    def validate(self) -> None:
        for name in ["borrowing_rate", "diffusion", "episode_rate", "episode_intensity"]:
            positive(getattr(self, name), name, zero=True)
        positive(self.episode_duration, "episode_duration")
        positive(self.scale, "scale")
        if self.distance not in {"patristic", "geography"}:
            raise ValueError("distance must be patristic or geography")
        if self.mode not in {"replace", "add"}:
            raise ValueError("mode must be replace or add")
        if (
            not math.isfinite(self.contamination_fraction)
            or not 0 <= self.contamination_fraction <= 1
        ):
            raise ValueError("contamination_fraction must be in [0,1]")


def displayed_network(tree: Tree, events: list[Event]) -> Tree:
    """Split recipient/donor edges at each event, then add zero-time hybrid edges.

    A displayed network records transfers, not inheritance probabilities. Multiple
    character transfers between the same lineages at the same time share one edge.
    """
    network = Tree.from_newick(tree.newick())
    originals = {n.name: n for n in network.traverse()}
    transfers = sorted({(e.time, e.donor, e.recipient) for e in events if e.donor})
    points: dict[str, list[float]] = {}
    for t, donor, recipient in transfers:
        points.setdefault(donor, []).append(t)
        points.setdefault(recipient, []).append(t)
    splits: dict[tuple[str, float], Node] = {}
    serial = 0
    for name, times in points.items():
        end = originals[name]
        parent = end.parent
        if parent is None:
            raise ValueError("horizontal event on root stem is not supported")
        start = parent.time
        for t in sorted(set(times)):
            if t < start - 1e-12 or t >= end.time:
                raise ValueError("horizontal event outside lineage lifetime")
            serial += 1
            split = Node(f"RET{serial:06d}", max(0, t - start), time=t)
            parent.children[parent.children.index(end)] = split
            split.parent = parent
            end.parent = None
            split.add(end)
            splits[name, t] = split
            parent, start = split, t
        end.length = end.time - start
    for t, donor, recipient in transfers:
        network.extra_edges.append((splits[donor, t], splits[recipient, t], 0.0))
    network.update_times()
    return network
