"""Timed character processes. States are sets; empty means known absence."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, TypeAlias

if TYPE_CHECKING:
    from .horizontal import HorizontalConfig

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

    gain_rate: float = 0.1
    loss_rate: float = 0.1
    root_frequency: float = 0.5
    k: int = 2
    ordered: bool = False
    rate_matrix: list[list[float]] | None = None
    block_size: int = 1
    block_correlation: float = 0.0
    omission_rate: float = 0.0
    max_events: int = 1000000

    def validate(self) -> None:
        if self.model not in {"lexical", "textual", "dollo", "binary", "mk"}:
            raise ValueError("unknown character model")
        if not isinstance(self.n_characters, int) or self.n_characters < 1:
            raise ValueError("n_characters must be a positive integer")
        for name in [
            "rate",
            "clock_sigma",
            "synonym_duration",
            "gain_rate",
            "loss_rate",
            "omission_rate",
        ]:
            positive(getattr(self, name), name, zero=True)
        for name in ["invariant", "p_parallel", "root_frequency", "block_correlation"]:
            if not 0 <= getattr(self, name) <= 1:
                raise ValueError(f"{name} must be in [0,1]")
        if self.alpha is not None:
            positive(self.alpha, "alpha")
        if self.categories < 0 or not isinstance(self.categories, int):
            raise ValueError("categories must be a nonnegative integer")
        if self.categories and self.alpha is None:
            raise ValueError("categories require alpha")
        if not 2 <= self.k <= 6 or not isinstance(self.k, int):
            raise ValueError("k must be an integer from 2 to 6")
        if not isinstance(self.block_size, int) or self.block_size < 1:
            raise ValueError("block_size must be a positive integer")
        if self.max_events < 1:
            raise ValueError("max_events must be positive")
        if self.synonym_duration and self.model not in {"lexical", "textual"}:
            raise ValueError("synonymy is only defined for lexical/textual characters")
        if self.block_correlation and self.model != "mk":
            raise ValueError("correlated blocks are only defined for Mk")
        if self.p_parallel and self.model not in {"lexical", "textual"}:
            raise ValueError("parallel pools are only defined for lexical/textual models")
        if self.omission_rate and self.model != "textual":
            raise ValueError("omission_rate is only defined for textual models")
        if self.rate_matrix is not None:
            if (
                self.model != "mk"
                or len(self.rate_matrix) != self.k
                or any(len(row) != self.k for row in self.rate_matrix)
            ):
                raise ValueError("rate_matrix must be a k by k Mk matrix")
            for i, row in enumerate(self.rate_matrix):
                for j, value in enumerate(row):
                    if i != j and self.ordered and abs(i - j) > 1 and value != 0:
                        raise ValueError("ordered Mk rates must connect adjacent states")
                    if i != j:
                        positive(value, "off-diagonal rate", zero=True)
                    elif value != 0:
                        raise ValueError("supply zero diagonal; exit rates are computed")
            if max(map(sum, self.rate_matrix)) <= 0:
                raise ValueError("Mk rate matrix must allow changes")
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
    history_tree: Tree | None = field(default=None, repr=False)


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


def evolve(
    tree: Tree,
    config: CharacterConfig | None = None,
    *,
    rng: RandomStream | None = None,
    seed: str | int | bytes = 0,
    horizontal: HorizontalConfig | None = None,
    horizontal_rng: RandomStream | None = None,
    geography_rng: RandomStream | None = None,
) -> Evolution:
    """Simulate timed CTMC events, optionally with contemporaneous donors."""
    from .horizontal import HorizontalConfig
    from .process import run_process
    from .rng import derive_seed

    cfg = config or CharacterConfig()
    if horizontal is not None and not isinstance(horizontal, HorizontalConfig):
        raise TypeError("horizontal must be HorizontalConfig or None")
    hcfg = horizontal or HorizontalConfig()
    return run_process(
        tree,
        cfg,
        rng or RandomStream(seed),
        hcfg,
        horizontal_rng or RandomStream(derive_seed(seed, "horizontal")),
        geography_rng or RandomStream(derive_seed(seed, "geography")),
    )
