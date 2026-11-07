"""Timed character processes. States are sets; empty means known absence."""
from __future__ import annotations

from dataclasses import dataclass, field
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
           rng: RandomStream | None = None, seed: str | int | bytes = 0,
           horizontal: object = None, horizontal_rng: RandomStream | None = None,
           geography_rng: RandomStream | None = None) -> Evolution:
    """Simulate timed CTMC events, optionally with contemporaneous donors."""
    from .horizontal import HorizontalConfig
    from .process import run_process
    from .rng import derive_seed
    cfg = config or CharacterConfig()
    hcfg = horizontal if isinstance(horizontal, HorizontalConfig) else HorizontalConfig()
    return run_process(tree,cfg,rng or RandomStream(seed),hcfg,
                       horizontal_rng or RandomStream(derive_seed(seed,'horizontal')),
                       geography_rng or RandomStream(derive_seed(seed,'geography')))
