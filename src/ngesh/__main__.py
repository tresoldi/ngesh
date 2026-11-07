"""ngesh tree, simulate and resumable factorial grid commands."""

from __future__ import annotations

import argparse
import copy
import itertools
import json
from pathlib import Path
import shutil
from typing import Any

from . import __version__
from .bundle import verify_bundle
from .config import load_config
from .models import TreeConfig, generate_tree
from .rng import derive_seed
from .simulation import PRESETS, simulate


def run_grid(config: dict[str, Any], out: Path) -> dict[str, int]:
    unknown = set(config) - {"seed", "preset", "replicates", "base", "grid"}
    if unknown:
        raise ValueError(f"unknown grid settings: {unknown}")
    seed = config.get("seed", 0)
    replicates = config.get("replicates", 1)
    if not isinstance(replicates, int) or replicates < 1:
        raise ValueError("replicates must be a positive integer")
    axes = config.get("grid", {})
    keys = sorted(axes)
    if any(not isinstance(axes[k], list) or not axes[k] for k in keys):
        raise ValueError("each grid axis must be a nonempty array")
    out.mkdir(parents=True, exist_ok=True)
    completed = skipped = 0
    for combination, values in enumerate(itertools.product(*(axes[k] for k in keys))):
        cfg = copy.deepcopy(config.get("base", {}))
        for key, value in zip(keys, values):
            section, sep, name = key.partition("__")
            if not sep:
                raise ValueError("grid axes must use section__parameter names")
            cfg.setdefault(section, {})[name] = value
        for replicate in range(replicates):
            run_seed = derive_seed(seed, f"grid:{combination}", replicate)
            dest = out / f"c{combination:04d}_r{replicate:05d}"
            if dest.exists():
                if verify_bundle(dest):
                    saved = json.loads((dest / "params.json").read_text())
                    if (
                        saved["master_seed"]["value"] == run_seed
                        and saved["config"] == cfg
                        and saved["preset"] == config.get("preset", "lexical")
                        and saved["version"] == __version__
                    ):
                        skipped += 1
                        continue
                raise ValueError(
                    f"{dest} exists but is incomplete, corrupt or from a different configuration"
                )
            temp = out / ("." + dest.name + ".partial")
            if temp.exists():
                shutil.rmtree(temp)
            simulate(seed=run_seed, preset=config.get("preset", "lexical"), config=cfg, out=temp)
            temp.rename(dest)
            completed += 1
    return {"completed": completed, "skipped": skipped}


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(prog="ngesh", description=__doc__)
    parser.add_argument("--version", action="version", version=__version__)
    subs = parser.add_subparsers(dest="command", required=True)
    tree = subs.add_parser("tree")
    tree.add_argument("--seed", default="0")
    tree.add_argument("--n-taxa", type=int, default=50)
    tree.add_argument("--birth", type=float, default=1)
    tree.add_argument("--death", type=float, default=0)
    tree.add_argument("--max-time", type=float)
    tree.add_argument("--time-only", action="store_true")
    tree.add_argument(
        "--shape", choices=["birth_death", "balanced", "caterpillar"], default="birth_death"
    )
    tree.add_argument("--out", type=Path)
    simulation = subs.add_parser("simulate")
    simulation.add_argument("--preset", choices=sorted(PRESETS), default="lexical")
    simulation.add_argument("--config", type=Path)
    simulation.add_argument("--seed", default="0")
    simulation.add_argument("--out", type=Path, required=True)
    grid = subs.add_parser("grid")
    grid.add_argument("config", type=Path)
    grid.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "tree":
            result = (
                generate_tree(
                    TreeConfig(
                        n_taxa=args.n_taxa,
                        birth=args.birth,
                        death=args.death,
                        max_time=args.max_time,
                        shape=args.shape,
                        stop="time" if args.time_only else "count",
                    ),
                    seed=args.seed,
                ).newick()
                + "\n"
            )
            if args.out:
                args.out.write_text(result, encoding="utf-8")
            else:
                print(result, end="")
        elif args.command == "simulate":
            simulate(
                seed=args.seed,
                preset=args.preset,
                config=load_config(args.config) if args.config else {},
                out=args.out,
            )
        else:
            print(json.dumps(run_grid(load_config(args.config), args.out), sort_keys=True))
    except (ValueError, TypeError, RuntimeError, OSError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
