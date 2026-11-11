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


def episode_cells(
    lineage_tree: Tree,
    events: list[Event],
    episodes: list,
    data: dict[str, list],
    columns: list[int] | None = None,
) -> list[dict]:
    """Where each episode's copied readings sit in the observed data.

    A witness *carries* a reading from episode E at character i when the latest event at i
    along its lineage (itself and its ancestors) is an E transfer that changed the state, and
    its observed cell still holds a copied state (``new`` minus ``old``, or ``new`` if the
    transfer replaced a state with an overlapping one). *Holders* are the other observed
    witnesses whose cell holds a copied state; *donor_side* are the holders that descend from
    the episode founder on the donor's side (the founders' living descendants were the only
    side recorded before). Donor lineages are often lost, especially in manuscript traditions,
    so holders, not donor_side, locate the other side of an episode in the observed data.

    ``columns[j]`` is the simulated character behind observed column ``j`` (identity if None);
    the returned ``character`` is the observed column.
    """
    nodes = {n.name: n for n in lineage_tree.traverse()}
    founders = {ep.id: (ep.initiator, ep.partner) for ep in episodes}
    latest: dict[tuple[str, int], Event] = {}
    for e in events:
        key = (e.recipient, e.character)
        if key not in latest or e.time >= latest[key].time:
            latest[key] = e
    path = {t: [n.name for n in nodes[t].ancestors()] if t in nodes else [t] for t in data}
    found: dict[tuple[str, int], dict] = {}
    for t, cells in data.items():
        for i, cell in enumerate(cells):
            if not cell:
                continue
            char = columns[i] if columns is not None else i
            last = None
            for name in path[t]:
                e = latest.get((name, char))
                if e is not None and (last is None or e.time > last.time):
                    last = e
            if last is None or last.kind != "contact" or not last.episode or last.old == last.new:
                continue
            copied = set(last.new) - set(last.old) or set(last.new)
            if not copied & set(cell):
                continue
            row = found.setdefault(
                (last.episode, i), {"states": set(), "carriers": set(), "donors": set()}
            )
            row["states"] |= copied
            row["carriers"].add(t)
            row["donors"].add(last.donor)
    out = []
    for (ep, i), row in sorted(found.items()):
        holders = sorted(
            t
            for t, cells in data.items()
            if t not in row["carriers"] and cells[i] and row["states"] & set(cells[i])
        )
        side = set()
        for donor in row["donors"]:
            anc = {n.name for n in nodes[donor].ancestors()} if donor in nodes else {donor}
            side |= {f for f in founders.get(ep, ()) if f in anc}
        donor_side = [t for t in holders if side & set(path[t])]
        out.append(
            {
                "episode": ep,
                "character": i,
                "states": sorted(row["states"]),
                "carriers": sorted(row["carriers"]),
                "holders": holders,
                "donor_side": donor_side,
            }
        )
    return out
