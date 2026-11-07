"""Chronological CTMC scheduler shared by independent and horizontal evolution."""

from __future__ import annotations

import bisect
import heapq
import math
from typing import Any

from .characters import CharacterConfig, Event, Evolution, character_rates
from .horizontal import HorizontalConfig
from .rng import RandomStream
from .tree import Node, Tree


def run_process(
    tree: Tree,
    cfg: CharacterConfig,
    rng: RandomStream,
    horizontal: HorizontalConfig,
    hrng: RandomStream,
    grng: RandomStream,
) -> Evolution:
    if tree.extra_edges:
        raise ValueError(
            "evolution requires a primary tree; horizontal events generate its network"
        )
    cfg.validate()
    horizontal.validate()
    if horizontal.mode == "add" and cfg.model not in {"binary", "dollo"}:
        raise ValueError("add transfers require a binary presence/absence model")
    rates = character_rates(cfg, rng)
    count = cfg.n_characters
    # Names are stable lineage IDs; sample observations retain tip names.
    tip_names = [n.name for n in tree.tips()]
    if len(set(tip_names)) != len(tip_names):
        raise ValueError("tip names must be unique")
    reserved = set(tip_names)
    for i, node in enumerate(tree.traverse()):
        if not node.name or (node.children and node.name in reserved):
            node.name = f"I{i:06d}"
        if node.children:
            reserved.add(node.name)
    if len({n.name for n in tree.traverse()}) != len(list(tree.traverse())):
        raise ValueError("lineage names must be unique")
    nodes = list(tree.traverse())
    branch_rates = {n.name: rng.lognormal(cfg.clock_sigma) for n in nodes}
    counters = [cfg.pool_size] * count
    states: dict[Node, list[frozenset[int]]] = {}
    deadlines: dict[Node, list[dict[int, float]]] = {}
    active: set[Node] = set()
    positions: dict[Node, tuple[float, float, float]] = {}
    events: list[Event] = []
    heap: list[tuple[float, int, int, str, Node, int, Any]] = []
    serial = 0
    q = cfg.rate_matrix or [
        [
            (1 / (cfg.k - 1) if i != j and (not cfg.ordered or abs(i - j) == 1) else 0.0)
            for j in range(cfg.k)
        ]
        for i in range(cfg.k)
    ]
    qmax = max(map(sum, q))
    shared: dict[int, float] = {}
    if cfg.model == "mk" and cfg.block_correlation:
        for b in range(0, count, cfg.block_size):
            shared[b] = min(rates[b : b + cfg.block_size]) * qmax * cfg.block_correlation
    root_cells = []
    for i in range(count):
        if cfg.model in {"dollo", "binary"}:
            root_cells.append(frozenset({1}) if rng.random() < cfg.root_frequency else frozenset())
        else:
            root_cells.append(frozenset({0}))

    def queue(
        t: float, priority: int, kind: str, node: Node, char: int = -1, value: Any = None
    ) -> None:
        nonlocal serial
        serial += 1
        heapq.heappush(heap, (t, priority, serial, kind, node, char, value))

    def schedule_change(t: float, node: Node, c: int) -> None:
        base = rates[c]
        if cfg.model == "dollo":
            base *= cfg.loss_rate
        elif cfg.model == "binary":
            base *= max(cfg.loss_rate, cfg.gain_rate)
        elif cfg.model == "mk":
            base = base * qmax - shared.get((c // cfg.block_size) * cfg.block_size, 0)
        elif cfg.model == "textual":
            base *= 1 + cfg.omission_rate
        rate = base * branch_rates[node.name]
        if rate:
            when = t + rng.exponential(rate)
            if when < node.time:
                queue(when, 4, "change", node, c)

    def schedule_borrow(t: float, node: Node) -> None:
        rate = horizontal.borrowing_rate * count
        if rate:
            when = t + hrng.exponential(rate)
            if when < node.time:
                queue(when, 5, "borrow", node)

    def position(node: Node, t: float) -> tuple[float, float]:
        x, y, last = positions[node]
        if t > last:
            sd = horizontal.diffusion * math.sqrt(t - last)
            if sd:
                x += sd * grng.normal()
                y += sd * grng.normal()
            positions[node] = x, y, t
        return x, y

    def choose_donor(node: Node, t: float, copying: bool = False) -> Node | None:
        candidates = sorted(
            (
                n
                for n in active
                if n is not node
                and n.parent is not None
                and n.parent.time <= t < n.time
                and (not copying or (n.parent is not node.parent and n.parent.time < t))
            ),
            key=lambda n: n.name,
        )
        if not candidates:
            return None
        distances = []
        if horizontal.distance == "geography":
            x, y = position(node, t)
            distances = [math.dist((x, y), position(d, t)) for d in candidates]
        else:
            distances = [tree.distance_at(node, d, t) for d in candidates]
        # Subtract min distance before exp, avoiding all-zero underflow.
        minimum = min(distances)
        weights = [math.exp(-(d - minimum) / horizontal.scale) for d in distances]
        cumulative = []
        total = 0.0
        for weight in weights:
            total += weight
            cumulative.append(total)
        return candidates[
            min(bisect.bisect_right(cumulative, hrng.random() * total), len(candidates) - 1)
        ]

    if cfg.model == "dollo":
        lengths = [n.length for n in nodes]
        total_length = sum(lengths)
        for c in range(count):
            if root_cells[c] or not rates[c] or not total_length:
                continue
            draw = rng.random() * total_length
            for node in nodes:
                if draw < node.length:
                    t = (node.parent.time if node.parent else 0) + draw
                    queue(t, 4, "gain", node, c)
                    break
                draw -= node.length

    def schedule_block(t: float, node: Node, b: int) -> None:
        rate = shared[b] * branch_rates[node.name]
        if rate:
            when = t + rng.exponential(rate)
            if when < node.time:
                queue(when, 4, "block", node, b)

    def mk_transition(old: frozenset[int]) -> frozenset[int]:
        current = next(iter(old))
        draw = rng.random() * qmax
        for new, rate in enumerate(q[current]):
            if draw < rate:
                return frozenset({new})
            draw -= rate
        return old

    for node in nodes:
        start = node.parent.time if node.parent else 0.0
        # Starts at a split must follow parent completion, including zero edges.
        # Preorder insertion resolves chains of zero-length boundaries.
        queue(start, 1, "start", node)
        if node.time > start:
            queue(node.time, 0, "end", node)
    while heap:
        if len(events) > cfg.max_events:
            raise RuntimeError("max_events exceeded; reduce rates or increase explicit limit")
        t, _, _, kind, node, c, value = heapq.heappop(heap)
        if kind == "start":
            if node.parent is None:
                states[node] = root_cells[:]
                deadlines[node] = [{} for _ in rates]
                positions[node] = (0, 0, t)
            else:
                states[node] = states[node.parent][:]
                deadlines[node] = [dict(d) for d in deadlines[node.parent]]
                x, y = position(node.parent, t)
                positions[node] = (x, y, t)
            if t < node.time and not node.attrs.get("sampled_tip"):
                active.add(node)
                for i in range(count):
                    schedule_change(t, node, i)
                    for state, when in deadlines[node][i].items():
                        if when < node.time:
                            queue(when, 3, "expire", node, i, state)
                schedule_borrow(t, node)
                for b in shared:
                    schedule_block(t, node, b)
                if (
                    horizontal.contamination_fraction
                    and node.parent
                    and len(node.parent.children) > 1
                    and not node.parent.attrs.get("sampled_ancestor")
                ):
                    queue(t, 2, "contaminate", node)
            continue
        if kind == "end":
            active.discard(node)
            if horizontal.distance == "geography":
                position(node, t)
            continue
        if kind in {"borrow", "contaminate"}:
            donor = choose_donor(node, t, copying=kind == "contaminate")
            if donor is not None:
                if kind == "borrow":
                    chars = [hrng.index(count)]
                elif horizontal.contamination_blocks:
                    width = round(horizontal.contamination_fraction * count)
                    start = hrng.index(count - width + 1)
                    chars = list(range(start, start + width))
                else:
                    chars = [
                        i for i in range(count) if hrng.random() < horizontal.contamination_fraction
                    ]
                for i in chars:
                    old = states[node][i]
                    new = states[donor][i]
                    if horizontal.mode == "add" and kind == "borrow":
                        new = old | new
                    states[node][i] = new
                    deadlines[node][i] = {}  # incoming reading becomes established
                    events.append(
                        Event(
                            t,
                            node.name,
                            i,
                            tuple(sorted(old)),
                            tuple(sorted(new)),
                            kind,
                            donor.name,
                        )
                    )
            if kind == "borrow":
                schedule_borrow(t, node)
            continue
        if kind == "block":
            for i in range(c, min(c + cfg.block_size, count)):
                old = states[node][i]
                new = mk_transition(old)
                if new != old:
                    states[node][i] = new
                    events.append(
                        Event(
                            t,
                            node.name,
                            i,
                            tuple(sorted(old)),
                            tuple(sorted(new)),
                            "correlated_change",
                        )
                    )
            schedule_block(t, node, c)
            continue
        old = states[node][c]
        if kind == "gain":
            new = frozenset({1})
            states[node][c] = new
            events.append(Event(t, node.name, c, tuple(sorted(old)), (1,), "gain"))
            continue
        if kind == "change" and cfg.model in {"binary", "dollo", "mk"}:
            if cfg.model == "mk":
                new = mk_transition(old)
                event_kind = "mk_change"
            else:
                if cfg.model == "binary":
                    actual = cfg.loss_rate if old else cfg.gain_rate
                    if rng.random() * max(cfg.loss_rate, cfg.gain_rate) >= actual:
                        schedule_change(t, node, c)
                        continue
                if cfg.model == "dollo" and not old:
                    schedule_change(t, node, c)
                    continue
                new = frozenset() if old else frozenset({1})
                event_kind = "loss" if old else "gain"
            states[node][c] = new
            if new != old:
                events.append(
                    Event(t, node.name, c, tuple(sorted(old)), tuple(sorted(new)), event_kind)
                )
            schedule_change(t, node, c)
            continue
        if (
            kind == "change"
            and cfg.model == "textual"
            and rng.random() < cfg.omission_rate / (1 + cfg.omission_rate)
        ):
            states[node][c] = frozenset()
            deadlines[node][c].clear()
            events.append(Event(t, node.name, c, tuple(sorted(old)), (), "omission"))
            schedule_change(t, node, c)
            continue
        if kind == "expire":
            if deadlines[node][c].get(value) != t:
                continue
            new = old - {value}
            deadlines[node][c].pop(value, None)
            event_kind = "synonym_loss"
        else:
            if cfg.p_parallel and rng.random() < cfg.p_parallel:
                new_id = 1 + rng.index(cfg.pool_size)
            else:
                counters[c] += 1
                new_id = counters[c]
            if cfg.synonym_duration:
                for s in old - {new_id}:
                    if s not in deadlines[node][c]:
                        when = t + rng.exponential(1 / cfg.synonym_duration)
                        deadlines[node][c][s] = when
                        if when < node.time:
                            queue(when, 3, "expire", node, c, s)
                deadlines[node][c].pop(new_id, None)
                new = old | {new_id}
            else:
                new = frozenset({new_id})
                deadlines[node][c].clear()
            event_kind = "replacement"
            schedule_change(t, node, c)
        states[node][c] = frozenset(new)
        events.append(Event(t, node.name, c, tuple(sorted(old)), tuple(sorted(new)), event_kind))
    all_states = {n.name: states[n] for n in nodes}
    truth = {n.name: states[n][:] for n in tree.tips()}
    return Evolution(
        tree,
        cfg,
        [f"C{i + 1:04d}" for i in range(count)],
        truth,
        events,
        rates,
        branch_rates,
        all_states,
    )
