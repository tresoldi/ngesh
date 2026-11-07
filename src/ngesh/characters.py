"""Timed character processes. States are sets; empty means known absence."""
from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import TypeAlias

from .rng import RandomStream, gamma_categories, positive
from .tree import Tree

Cell: TypeAlias = frozenset[int] | None
Matrix: TypeAlias = dict[str, list[Cell]]


@dataclass
class CharacterConfig:
    model: str = "lexical"
    n_characters: int = 200
    rate: float = 0.1
    alpha: float | None = None
    categories: int = 0
    invariant: float = 0.0
    clock_sigma: float = 0.0
    p_parallel: float = 0.0
    pool_size: int = 3
    synonym_duration: float = 0.0

    def validate(self) -> None:
        if self.model != "lexical":
            raise ValueError("unknown character model")
        if not isinstance(self.n_characters, int) or self.n_characters < 1:
            raise ValueError("n_characters must be a positive integer")
        for name in ["rate", "clock_sigma", "synonym_duration"]:
            positive(getattr(self, name), name, zero=True)
        for name in ["invariant", "p_parallel"]:
            if not 0 <= getattr(self, name) <= 1:
                raise ValueError(f"{name} must be in [0,1]")
        if self.alpha is not None:
            positive(self.alpha, "alpha")
        if self.categories < 0 or not isinstance(self.categories, int):
            raise ValueError("categories must be a nonnegative integer")
        if self.categories and self.alpha is None:
            raise ValueError("categories require alpha")
        if self.pool_size < 1:
            raise ValueError("pool_size must be positive")


@dataclass
class Event:
    time: float
    recipient: str
    character: int
    old: tuple[int, ...]
    new: tuple[int, ...]
    kind: str = "replacement"
    donor: str = ""


@dataclass
class Evolution:
    tree: Tree
    config: CharacterConfig
    groups: list[str]
    truth: Matrix
    events: list[Event]
    rates: list[float]
    branch_rates: dict[str, float]
    all_states: Matrix = field(default_factory=dict, repr=False)


def character_rates(cfg: CharacterConfig, rng: RandomStream) -> list[float]:
    categories = gamma_categories(cfg.alpha, cfg.categories) if cfg.categories else None
    rates = []
    for _ in range(cfg.n_characters):
        if cfg.invariant and rng.random() < cfg.invariant:
            rates.append(0.0)
            continue
        multiplier = 1.0
        if categories:
            multiplier = categories[rng.index(len(categories))]
        elif cfg.alpha is not None:
            multiplier = rng.gamma(cfg.alpha, 1 / cfg.alpha)
        rates.append(cfg.rate * multiplier)
    return rates


def evolve(tree: Tree, config: CharacterConfig | None = None, *,
           rng: RandomStream | None = None, seed: str | int | bytes = 0) -> Evolution:
    """Independent Poisson replacements; branch state and synonym expiry inherited.

    Infinite-alleles identifiers are unique within each group. Parallel pools
    occupy 1..pool_size; innovations occupy larger integers. Root state is 0.
    """
    cfg = config or CharacterConfig()
    cfg.validate()
    rng = rng or RandomStream(seed)
    rates = character_rates(cfg, rng)
    counters = [cfg.pool_size] * cfg.n_characters
    states: dict[str, list[frozenset[int]]] = {}
    expiries: dict[str, list[dict[int, float]]] = {}
    branch_rates = {}
    events = []
    names = [n.name for n in tree.traverse()]
    if len(set(names)) != len(names):
        # Sampled observations share names with their parent; internal IDs differ.
        for i, node in enumerate(tree.traverse()):
            if node.children:
                node.name = f"I{i:05d}"
    for i, node in enumerate(tree.traverse()):
        if not node.name:
            node.name = f"I{i:05d}"
    for node in tree.traverse():
        if node.parent is None:
            cells = [frozenset({0}) for _ in rates]
            expires = [{} for _ in rates]
            start = 0.0
        else:
            cells = states[node.parent.name][:]
            expires = [dict(x) for x in expiries[node.parent.name]]
            start = node.parent.time
        multiplier = rng.lognormal(cfg.clock_sigma)
        branch_rates[node.name] = multiplier
        for c, rate in enumerate(rates):
            current = set(cells[c])
            deadlines = expires[c]
            t = start
            while True:
                replacement = t + rng.exponential(rate * multiplier) if rate else math.inf
                expiry = min(deadlines.values(), default=math.inf)
                when = min(replacement, expiry)
                if when >= node.time:
                    break
                old = tuple(sorted(current))
                if expiry <= replacement:
                    dead = [s for s, end in deadlines.items() if end <= when]
                    current.difference_update(dead)
                    for s in dead:
                        del deadlines[s]
                    kind = "synonym_loss"
                else:
                    if cfg.p_parallel and rng.random() < cfg.p_parallel:
                        new = 1 + rng.index(cfg.pool_size)
                    else:
                        counters[c] += 1
                        new = counters[c]
                    if cfg.synonym_duration:
                        for s in current - {new}:
                            if s not in deadlines:
                                deadlines[s] = when + rng.exponential(1 / cfg.synonym_duration)
                        deadlines.pop(new, None)
                        current.add(new)
                    else:
                        current = {new}
                        deadlines.clear()
                    kind = "replacement"
                events.append(Event(when, node.name, c, old, tuple(sorted(current)), kind))
                t = when
            cells[c] = frozenset(current)
        states[node.name] = cells
        expiries[node.name] = expires
    events.sort(key=lambda e: (e.time, e.recipient, e.character, e.kind))
    truth = {n.name: states[n.name][:] for n in tree.tips()}
    return Evolution(tree, cfg, [f"C{i + 1:04d}" for i in range(cfg.n_characters)],
                     truth, events, rates, branch_rates, states)
