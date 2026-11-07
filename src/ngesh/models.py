"""Continuous-time tree models and deterministic control shapes."""
from __future__ import annotations

from dataclasses import dataclass
import math

from .rng import RandomStream, Seed, positive
from .tree import Node, Tree


@dataclass
class TreeConfig:
    birth: float = 1.0
    death: float = 0.0
    n_taxa: int | None = 50
    max_time: float | None = None
    polytomy: float = 0.0
    sampling: float = 1.0
    ancestor_rate: float = 0.0
    full: bool = False
    method: str = "fast"
    shape: str = "birth_death"
    newick: str | None = None
    height: float | None = None
    mean_height: float | None = None
    max_attempts: int = 1000
    max_lineages: int = 100000
    copy_sigma: float = 0.0
    burst_probability: float = 0.0
    burst_mean: float = 0.0

    def validate(self) -> None:
        for name in ["birth", "death", "polytomy", "ancestor_rate", "copy_sigma", "burst_mean"]:
            positive(getattr(self, name), name, zero=True)
        for name in ["sampling", "burst_probability"]:
            value = getattr(self, name)
            if not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError(f"{name} must be in [0,1]")
        if self.n_taxa is not None and (not isinstance(self.n_taxa, int) or self.n_taxa < 2):
            raise ValueError("n_taxa must be an integer >= 2, or None")
        if self.max_time is not None:
            positive(self.max_time, "max_time")
        if self.n_taxa is None and self.max_time is None and self.newick is None:
            raise ValueError("a stopping criterion is required")
        if self.method not in {"fast", "standard"}:
            raise ValueError("method must be fast or standard")
        if self.shape not in {"birth_death", "balanced", "caterpillar", "newick"}:
            raise ValueError("unknown tree shape")
        if self.height is not None and self.mean_height is not None:
            raise ValueError("choose height or mean_height")
        if self.max_attempts < 1 or self.max_lineages < 2:
            raise ValueError("simulation limits must be positive")


def control_tree(n: int, shape: str = "balanced") -> Tree:
    """Unit-height ultrametric controls; balanced splits differ by at most one."""
    if n < 2 or shape not in {"balanced", "caterpillar"}:
        raise ValueError("control tree needs n>=2 and balanced/caterpillar shape")
    root = Node("N0")
    stack = [(root, n, 0.0)]
    serial = 0
    while stack:
        node, count, start = stack.pop()
        left = count // 2 if shape == "balanced" else 1
        for size in [left, count - left]:
            serial += 1
            end = 1.0 if size == 1 else start + (1 - start) / 2
            child = node.add(Node(f"N{serial}", end - start))
            if size > 1:
                stack.append((child, size, end))
    tree = Tree(root)
    for i, tip in enumerate(tree.tips(), 1):
        tip.name = f"T{i:04d}"
    return tree


def sample_ancestors(tree: Tree, rate: float, rng: RandomStream) -> int:
    """Poisson sampling per lineage-time; retain continuation and zero-length sample."""
    positive(rate, "ancestor_rate", zero=True)
    if rate == 0:
        return 0
    count = 0
    for end in list(tree.traverse()):
        parent = end.parent
        if parent is None or end.attrs.get("sampled_tip"):
            continue
        start, finish = parent.time, end.time
        t = start + rng.exponential(rate)
        while t < finish:
            count += 1
            sample = Node(f"SA{count:04d}", t - start, attrs={"sampled_ancestor": True})
            parent.children[parent.children.index(end)] = sample
            sample.parent = parent
            end.parent = None
            sample.add(end)
            sample.add(Node(sample.name, 0, attrs={"sampled_tip": True}))
            parent, start = sample, t
            t += rng.exponential(rate)
        end.length = finish - start
    tree.update_times()
    return count


def ancestor_internal(tree: Tree) -> Tree:
    """Encode observations as labelled internal nodes instead of duplicate tips."""
    result = Tree.from_newick(tree.newick())
    source = list(tree.traverse())
    for old, new in zip(source, result.traverse()):
        new.attrs = dict(old.attrs)
    for node in list(result.traverse()):
        if node.attrs.get("sampled_tip") and node.parent:
            node.parent.children.remove(node)
    return result


def generate_tree(config: TreeConfig | None = None, *, seed: Seed = 0,
                  rng: RandomStream | None = None) -> Tree:
    """Birth–death process, conditioned on surviving and attaining the target.

    Exact counts stop at the first hitting time; with a deadline, successful
    trajectories must reach that count before the deadline. Time-only runs are
    conditioned on at least two extant tips. Extinction causes bounded retries.
    Both method names use the same event-driven kernel in 2.x.
    """
    cfg = config or TreeConfig()
    cfg.validate()
    rng = rng or RandomStream(seed)
    if cfg.newick is not None:
        tree = Tree.from_newick(cfg.newick)
        if tree.extra_edges:
            raise ValueError("input truth must be a tree")
    elif cfg.shape in {"balanced", "caterpillar"}:
        if cfg.n_taxa is None:
            raise ValueError("control shapes require n_taxa")
        tree = control_tree(cfg.n_taxa, cfg.shape)
    elif cfg.shape == "newick":
        raise ValueError("newick shape requires newick text")
    else:
        if cfg.birth == 0:
            raise ValueError("birth must be positive for a branching process")
        tree = _birth_death(cfg, rng)
    sample_ancestors(tree, cfg.ancestor_rate, rng)
    if not cfg.full or cfg.sampling < 1:
        keep = {n.name for n in tree.tips() if n.attrs.get("sampled_tip") or
                (not n.attrs.get("extinct", False) and
                 (cfg.sampling == 1 or rng.random() < cfg.sampling))}
        if cfg.full:
            keep.update(n.name for n in tree.tips() if n.attrs.get("extinct", False))
        if not keep:
            raise ValueError("sampling retained no observations; change seed or sampling probability")
        tree = tree.prune(keep)
    if cfg.height is not None:
        tree.rescale(cfg.height)
    if cfg.mean_height is not None:
        tree.rescale(cfg.mean_height, "mean")
    return tree


def _birth_death(cfg: TreeConfig, rng: RandomStream) -> Tree:
    for _ in range(cfg.max_attempts):
        root = Node("N0")
        active = [root]
        now, serial = 0.0, 0
        root.attrs["copy_rate"] = rng.lognormal(cfg.copy_sigma)
        success = False
        while active:
            rates = [cfg.birth * n.attrs["copy_rate"] + cfg.death for n in active]
            total = sum(rates)
            wait = rng.exponential(total)
            if cfg.max_time is not None and now + wait >= cfg.max_time:
                now = cfg.max_time
                success = cfg.n_taxa is None and len(active) >= 2
                break
            now += wait
            draw = rng.random() * total
            idx = 0
            while idx < len(rates) - 1 and draw >= rates[idx]:
                draw -= rates[idx]
                idx += 1
            node = active.pop(idx)
            node.length = now - node.time
            node.attrs["extinct"] = True
            birth = cfg.birth * node.attrs["copy_rate"]
            if rng.random() < birth / (birth + cfg.death):
                offspring = 2 + rng.poisson(cfg.polytomy)
                if cfg.burst_probability and rng.random() < cfg.burst_probability:
                    offspring += rng.poisson(cfg.burst_mean)
                if cfg.n_taxa is not None:
                    offspring = min(offspring, cfg.n_taxa - len(active))
                for _ in range(offspring):
                    serial += 1
                    child = node.add(Node(f"N{serial}", time=now,
                                         attrs={"extinct": False, "copy_rate": rng.lognormal(cfg.copy_sigma)}))
                    active.append(child)
            if cfg.n_taxa is not None and len(active) == cfg.n_taxa:
                success = True
                break
            if len(active) > cfg.max_lineages:
                raise RuntimeError("max_lineages exceeded; reduce rates or time horizon")
        if not success:
            continue
        for node in active:
            node.length = now - node.time
            node.attrs["extinct"] = False
        tree = Tree(root)
        for i, node in enumerate(tree.tips(), 1):
            node.name = f"T{i:04d}"
        return tree
    raise RuntimeError("no surviving tree met the stop criteria within max_attempts")
