"""Predeclared fitting/held-out validation of marginal preset statistics.

No tree-inference method enters selection. Family splits are fixed before fitting;
large datasets have equal weight. Unavailable sources cannot be synthesized.
"""

from __future__ import annotations

import argparse
import copy
import itertools
import json
from pathlib import Path
import statistics

try:
    from .calibrate import compact_summary
except ImportError:
    from calibrate import compact_summary
from ngesh import simulate
from ngesh.bundle import json_text
from ngesh.observation import ObservationConfig, observe
from ngesh.rng import RandomStream, derive_seed
from ngesh.simulation import PRESETS


# Frozen alpha3 defaults: rerunning a fit must not change its baseline after promotion.
BASELINES = {
    "lexical": {
        "characters": {"alpha": 0.7, "model": "lexical", "rate": 0.3},
        "observation": {"missing": 0.15},
        "tree": {},
    },
    "genecontent": {
        "characters": {
            "loss_rate": 0.1,
            "model": "dollo",
            "n_characters": 500,
            "rate": 1.0,
            "root_frequency": 0.5,
        },
        "horizontal": {"mode": "add"},
        "tree": {},
    },
    "morphology": {
        "characters": {"alpha": 0.7, "k": 4, "model": "mk", "n_characters": 100, "rate": 0.2},
        "tree": {},
    },
}

BUDGET = 150


def score_features(s: dict, preset: str) -> dict[str, float]:
    base = {"distance": s["distance_distribution"]["mean"], "delta": s["delta"] or 0.0}
    if preset == "genecontent":
        p = s["prevalence_distribution"]
        return base | {
            "prevalence_mean": p["mean"],
            "prevalence_q75": p["quantiles"][4],
            "prevalence_q90": p["quantiles"][5],
            "trait_singletons": s["trait_singleton_fraction"],
        }
    return base | {
        "singletons": s["singleton_fraction"],
        "classes": s["class_distribution"]["mean"] / (s["taxa"] if preset == "lexical" else 6),
        "constant": s["constant_fraction"],
    }


def coverage_fit(records: list[dict]) -> dict:
    mean = statistics.mean(r["summary"]["coverage_distribution"]["mean"] for r in records)
    variance = statistics.mean(
        r["summary"]["coverage_distribution"]["variance"] or 0 for r in records
    )
    # Match within-family taxon coverage variance, not between-family offsets.
    if variance <= 1e-9 or mean >= 0.999999:
        return {"missing": 1 - mean}
    variance = min(variance, mean * (1 - mean) * 0.99)
    concentration = mean * (1 - mean) / variance - 1
    return {
        "missing": 0.0,
        "coverage_a": mean * concentration,
        "coverage_b": (1 - mean) * concentration,
    }


def simulate_summary(preset: str, config: dict, target: dict, seed: int) -> dict:
    cfg = copy.deepcopy(BASELINES[preset])
    cfg.setdefault("tree", {}).update({"height": None, "mean_height": None})
    cfg.setdefault("observation", {}).update({"coverage_a": None, "coverage_b": None})
    if preset == "genecontent":
        cfg.setdefault("characters", {})["alpha"] = None
    for section, values in config.items():
        cfg.setdefault(section, {}).update(values)
    cfg.setdefault("tree", {})["n_taxa"] = target["summary"]["taxa"]
    cfg.setdefault("characters", {})["n_characters"] = min(BUDGET, target["summary"]["characters"])
    result = simulate(seed=seed, preset=preset, config=cfg)
    effective = PRESETS[preset].get("observation", {}) | cfg.get("observation", {})
    observed = observe(
        result, ObservationConfig(**effective), rng=RandomStream(derive_seed(seed, "observation"))
    )
    matrix = observed.matrix
    if preset == "genecontent":
        # The empirical clusters are conditioned on at least one ingroup gene.
        # Apply that same explicit projection only for calibration, not exports.
        indices = [i for i in range(len(result.groups)) if any(row[i] for row in matrix.values())]
        matrix = {t: [row[i] for i in indices] for t, row in matrix.items()}
    return compact_summary(matrix, binary=preset == "genecontent")


def fit(
    preset: str, targets: Path, out: Path, replicates: int = 1, validation_replicates: int = 3
) -> dict:
    files = sorted(targets.glob("*.json"))
    if preset == "genecontent":
        files = [targets / "yeast-train.json", targets / "yeast-heldout.json"]
    records = [dict(json.loads(p.read_text()), target_file=p.name) for p in files]
    if not records:
        raise ValueError("no empirical target summaries")
    if preset == "genecontent":
        train, heldout = records[:1], records[1:]
    else:
        train = [r for i, r in enumerate(records) if i % 4 != 0]
        heldout = [r for i, r in enumerate(records) if i % 4 == 0]
    if not train or not heldout:
        raise ValueError("fitting requires nonempty training and held-out sets")
    sampled = (
        [train[round(i * (len(train) - 1) / 7)] for i in range(8)] if len(train) > 8 else train
    )
    nuisance = coverage_fit(train) if preset != "genecontent" else {"missing": 0.0}
    if preset == "lexical":
        candidates = [
            {"characters": {"rate": r, "alpha": a}, "observation": nuisance}
            for r, a in itertools.product([0.1, 0.2, 0.3, 0.45, 0.6], [0.3, 0.7, 1.5, 3])
        ]
    elif preset == "morphology":
        candidates = [
            {
                "tree": {"crown_height": 1.0},
                "characters": {"rate": r, "alpha": a, "k": k},
                "observation": nuisance,
            }
            for r, a, k in itertools.product([0.1, 0.3, 0.6, 1, 2], [0.3, 0.7, 2], [2, 3, 4, 5, 6])
        ]
    elif preset == "genecontent":
        candidates = [
            {
                "tree": {"crown_height": 1.0},
                "characters": {"rate": 1.0, "loss_rate": loss, "root_frequency": root, "alpha": a},
                "observation": nuisance,
            }
            for loss, root, a in itertools.product(
                [0.1, 0.3, 1, 3], [0, 0.005, 0.01, 0.03], [0.3, 1, 3]
            )
        ]
    else:
        raise ValueError(
            "stemma fitting awaits the specified private traditions and reference histories"
        )

    def evaluate(config: dict, rows: list[dict], reps: int, stage: str) -> dict:
        scored = []
        for target in rows:
            actual = score_features(target["summary"], preset)
            draws = []
            summaries = []
            for rep in range(reps):
                seed = derive_seed(
                    "ngesh-fit-v2", stage + ":" + preset + ":" + target["target_file"], rep
                )
                summary = simulate_summary(preset, config, target, seed)
                draws.append(score_features(summary, preset))
                summaries.append(summary)
            average = {k: statistics.mean(d[k] for d in draws) for k in actual}
            loss = statistics.mean((average[k] - actual[k]) ** 2 for k in actual)
            scored.append(
                {
                    "target": target["target_file"],
                    "target_features": actual,
                    "simulated_features": average,
                    "score": loss,
                    "summaries": summaries,
                }
            )
        return {"score": statistics.mean(r["score"] for r in scored), "datasets": scored}

    trials = []
    for candidate in candidates:
        result = evaluate(candidate, sampled, replicates, "fit")
        trials.append({"config": candidate, "score": result["score"]})
    trials.sort(key=lambda t: t["score"])
    best = trials[0]
    # All holdout evaluations happen after selection; no test-score reselection.
    baseline = evaluate({}, heldout, validation_replicates, "validation")
    validation = evaluate(best["config"], heldout, validation_replicates, "validation")
    training_validation = evaluate(best["config"], train, validation_replicates, "validation")
    result = {
        "schema": "ngesh-fit-v2",
        "preset": preset,
        "selection_features": list(score_features(train[0]["summary"], preset)),
        "seed_master": "ngesh-fit-v2",
        "replicates_per_fit": replicates,
        "replicates_per_validation": validation_replicates,
        "maximum_simulated_characters": BUDGET,
        "training_targets": [r["target_file"] for r in train],
        "fitting_subset": [r["target_file"] for r in sampled],
        "heldout_targets": [r["target_file"] for r in heldout],
        "coverage_method": "mean and within-family variance matched on training targets only",
        "fit_kind": "character holdout on shared taxa"
        if preset == "genecontent"
        else "dataset holdout",
        "preexisting_preset": BASELINES[preset],
        "best": best,
        "trials": trials,
        "heldout_baseline": baseline,
        "heldout_selected": validation,
        "training_selected": training_validation,
        "limitations": [
            "coarse grid, finite Monte Carlo budget; no accuracy claim",
            "marginal summaries do not identify borrowing, contamination or absolute historical rates",
            "Beta coverage captures incomplete applicability/attestation, not its causal mechanism",
        ],
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json_text(result), encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preset", choices=["lexical", "morphology", "genecontent"], required=True)
    parser.add_argument("--targets", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=1)
    parser.add_argument("--validation-replicates", type=int, default=3)
    args = parser.parse_args()
    result = fit(args.preset, args.targets, args.out, args.replicates, args.validation_replicates)
    print(
        json.dumps(
            {
                "preset": args.preset,
                "best": result["best"],
                "heldout_baseline": result["heldout_baseline"]["score"],
                "heldout_selected": result["heldout_selected"]["score"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
