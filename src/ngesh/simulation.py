"""Public orchestration and stable component/replicate seeds."""

from __future__ import annotations

from dataclasses import asdict, replace
from pathlib import Path
from typing import Any

from .bundle import write_bundle
from .characters import CharacterConfig, Evolution, evolve
from .models import TreeConfig, generate_tree
from .rng import RandomStream, Seed, derive_seed
from .horizontal import HorizontalConfig
from .observation import ObservationConfig, observe


PRESETS: dict[str, dict[str, Any]] = {
    "lexical": {
        "tree": {"crown_height": 1.0},
        "characters": {"model": "lexical", "rate": 0.3, "alpha": 0.7},
        "observation": {"missing": 0.15},
    },
    "genecontent": {
        "tree": {"crown_height": 1.0},
        "characters": {
            "model": "dollo",
            "rate": 1.0,
            "loss_rate": 1.0,
            "alpha": 0.3,
            "n_characters": 500,
            "root_frequency": 0.03,
        },
        "horizontal": {"mode": "add"},
    },
    "morphology": {
        "tree": {"crown_height": 1.0},
        "characters": {"model": "mk", "rate": 1.0, "k": 2, "n_characters": 100, "alpha": 2.0},
        "observation": {"coverage_a": 7.702921781245104, "coverage_b": 1.6842147553315705},
    },
    "stemma": {
        "tree": {
            "birth": 1.0,
            "death": 0.8,
            "ancestor_rate": 0.2,
            "copy_sigma": 0.6,
            "crown_height": 1.0,
        },
        "characters": {"model": "textual", "rate": 0.2, "p_parallel": 0.2, "omission_rate": 0.05},
    },
}


def simulate(
    *,
    seed: Seed = 0,
    preset: str = "lexical",
    config: dict[str, Any] | None = None,
    out: str | Path | None = None,
) -> Evolution:
    if preset not in PRESETS:
        raise ValueError(f"unknown preset: {preset}")
    config = config or {}
    unknown = set(config) - {"tree", "characters", "output", "horizontal", "observation"}
    if unknown:
        raise ValueError(f"unknown configuration sections: {sorted(unknown)}")
    base = PRESETS[preset]
    tree_options = dict(base.get("tree", {}))
    supplied_tree = config.get("tree", {})
    if any(
        supplied_tree.get(k) is not None
        for k in ("height", "mean_height", "crown_height", "origin_height")
    ):
        for k in ("height", "mean_height", "crown_height", "origin_height"):
            tree_options.pop(k, None)
    tree_cfg = TreeConfig(**(tree_options | supplied_tree))
    tree_cfg.validate()
    char_cfg = CharacterConfig(**(base.get("characters", {}) | config.get("characters", {})))
    seeds = {
        component: derive_seed(seed, component)
        for component in [
            "tree",
            "characters",
            "observation",
            "horizontal",
            "geography",
            "tree_sampling",
        ]
    }
    tree = generate_tree(
        replace(
            tree_cfg,
            full=True,
            sampling=1,
            height=None,
            mean_height=None,
            crown_height=None,
            origin_height=None,
        ),
        seed=seeds["tree"],
    )
    hcfg = HorizontalConfig(**(base.get("horizontal", {}) | config.get("horizontal", {})))
    ocfg = ObservationConfig(**(base.get("observation", {}) | config.get("observation", {})))
    if char_cfg.model in {"binary", "dollo"} and (ocfg.lumping or ocfg.splitting):
        raise ValueError("class lumping/splitting is not defined for binary presence/absence")
    srng = RandomStream(seeds["tree_sampling"])
    for sampling_attempts in range(1, tree_cfg.max_attempts + 1):
        selected = [
            n
            for n in tree.tips()
            if n.attrs.get("sampled_tip")
            or (n.attrs.get("extinct", False) and tree_cfg.full)
            or (
                not n.attrs.get("extinct", False)
                and (tree_cfg.sampling == 1 or srng.random() < tree_cfg.sampling)
            )
        ]
        if selected:
            break
    else:
        raise ValueError("tree sampling retained no observations within max_attempts")
    scaling = [
        ("crown", tree_cfg.crown_height),
        ("origin", tree_cfg.origin_height),
        ("origin", tree_cfg.height),
        ("mean", tree_cfg.mean_height),
    ]
    for mode, target in scaling:
        if target is None:
            continue
        from .rng import positive

        positive(target, "rescaling target")
        if mode == "crown":
            common = set(selected[0].ancestors())
            for tip in selected[1:]:
                common.intersection_update(tip.ancestors())
            crown = next(n for n in selected[0].ancestors() if n in common)
            base = max(n.time for n in selected) - crown.time
        elif mode == "origin":
            base = max(n.time for n in selected)
        else:
            base = sum(n.time for n in selected) / len(selected)
        if base <= 0:
            raise ValueError("cannot rescale a zero-height observation tree")
        for n in tree.traverse():
            n.length *= target / base
        tree.update_times()
    result = evolve(
        tree,
        char_cfg,
        seed=seeds["characters"],
        horizontal=hcfg,
        horizontal_rng=RandomStream(seeds["horizontal"]),
        geography_rng=RandomStream(seeds["geography"]),
    )
    keep = {n.name for n in selected}
    if keep != set(result.truth):
        result.history_tree = result.tree
        result.tree = result.tree.prune(keep)
        result.truth = {t: row for t, row in result.truth.items() if t in keep}
    observed = observe(result, ocfg, rng=RandomStream(seeds["observation"]))
    if out is not None:
        master = {
            "type": type(seed).__name__,
            "value": seed.hex() if isinstance(seed, bytes) else seed,
        }
        output = config.get("output", {})
        if set(output) - {"include_ci", "max_quartets"}:
            raise ValueError("unknown output configuration")
        write_bundle(
            out,
            result,
            params={
                "preset": preset,
                "master_seed": master,
                "derived_seeds": seeds,
                "rng_integer_seeds": {c: derive_seed(value) for c, value in seeds.items()},
                "tree_config": asdict(tree_cfg),
                "tree_sampling_attempts": sampling_attempts,
                "horizontal_config": asdict(hcfg),
                "observation_config": asdict(ocfg),
                "config": config,
            },
            observed=observed.matrix,
            truth_tree=observed.tree,
            observation_log=observed.log,
            retained_groups=observed.retained_groups,
            include_ci=output.get("include_ci", True),
            max_quartets=output.get("max_quartets", 1000),
        )
    return result
