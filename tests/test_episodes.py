"""Contact episodes: sustained transfer between two contemporaneous lineages."""

import csv
import json
import tempfile
from pathlib import Path

import pytest

from ngesh import simulate
from ngesh.horizontal import HorizontalConfig

EPISODES = {"episode_rate": 1.0, "episode_duration": 0.3, "episode_intensity": 0.5}


def run(tmp, seed, horizontal, n_taxa=12, n_characters=60):
    out = Path(tmp) / str(seed)
    simulate(
        preset="lexical",
        seed=seed,
        config={
            "tree": {"n_taxa": n_taxa},
            "characters": {"n_characters": n_characters},
            "horizontal": horizontal,
        },
        out=out,
    )
    return out


def founder(lineages, name, episode):
    """The episode founder (initiator or partner) that lineage `name` descends from."""
    node = name
    while node is not None:
        if node in (episode["initiator"], episode["partner"]):
            return node
        node = lineages[node]["parent"]
    return None


def rows(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_transfers_stay_inside_their_episode():
    with tempfile.TemporaryDirectory() as tmp:
        for seed in range(5):
            out = run(tmp, seed, EPISODES)
            lineages = {x["id"]: x for x in json.loads((out / "lineages.json").read_text())}
            episodes = {e["episode"]: e for e in rows(out / "episodes.tsv")}
            assert episodes, "expected at least one episode"
            counts = dict.fromkeys(episodes, 0)
            for event in rows(out / "events.tsv"):
                if event["kind"] != "contact":
                    assert event["episode"] == ""
                    continue
                episode = episodes[event["episode"]]
                counts[event["episode"]] += 1
                t = float(event["time"])
                assert float(episode["start"]) <= t < float(episode["end"])
                sides = {
                    founder(lineages, event["donor"], episode),
                    founder(lineages, event["recipient"], episode),
                }
                assert sides == {episode["initiator"], episode["partner"]}
                for lineage in (event["donor"], event["recipient"]):
                    assert lineages[lineage]["start"] <= t < lineages[lineage]["end"]
            assert counts == {k: int(v["transfers"]) for k, v in episodes.items()}


def test_transfer_count_matches_intensity():
    # E[transfers] = intensity * characters * total episode duration.
    observed = expected = 0.0
    with tempfile.TemporaryDirectory() as tmp:
        for seed in range(20):
            out = run(tmp, f"rate-{seed}", EPISODES, n_characters=40)
            for episode in rows(out / "episodes.tsv"):
                observed += int(episode["transfers"])
                # an episode ongoing at the present only runs until the present
                present = max(
                    float(x["end"]) for x in json.loads((out / "lineages.json").read_text())
                )
                end = min(float(episode["end"]), present)
                expected += 0.5 * 40 * (end - float(episode["start"]))
    assert expected > 100
    assert abs(observed - expected) < 4 * expected**0.5


def test_one_way_episodes_only_import_into_the_initiator():
    with tempfile.TemporaryDirectory() as tmp:
        out = run(tmp, "oneway", dict(EPISODES, episode_symmetric=False))
        episodes = {e["episode"]: e for e in rows(out / "episodes.tsv")}
        contact = [e for e in rows(out / "events.tsv") if e["kind"] == "contact"]
        assert contact
        lineages = {x["id"]: x for x in json.loads((out / "lineages.json").read_text())}
        for event in contact:
            episode = episodes[event["episode"]]
            assert founder(lineages, event["recipient"], episode) == episode["initiator"]
            assert founder(lineages, event["donor"], episode) == episode["partner"]


def test_deterministic_and_zero_rate_is_identical():
    with tempfile.TemporaryDirectory() as tmp:
        a = run(tmp, "same", EPISODES)
        b = run(Path(tmp) / "again", "same", EPISODES)
        for name in ["events.tsv", "episodes.tsv", "data.tsv", "history.tsv"]:
            assert (a / name).read_text() == (b / name).read_text()
        plain = run(Path(tmp) / "p", "zero", {"borrowing_rate": 0.1})
        zero = run(Path(tmp) / "z", "zero", {"borrowing_rate": 0.1, "episode_rate": 0.0})
        for name in ["data.tsv", "history.tsv", "events.tsv", "truth.nwk"]:
            assert (plain / name).read_text() == (zero / name).read_text()
        assert not (zero / "episodes.tsv").exists()


def test_episodes_continue_through_splits():
    # Episodes are not cut at lineage splits. With pure-birth trees no side dies
    # out, so durations are exponential, right-censored only at the present:
    # the mean estimate is total exposure over completed episodes.
    exposure, completed, outlived = 0.0, 0, 0
    total = 0
    with tempfile.TemporaryDirectory() as tmp:
        for seed in range(20):
            out = run(tmp, f"split-{seed}", dict(EPISODES, episode_rate=0.5), n_taxa=40)
            lineages = {x["id"]: x for x in json.loads((out / "lineages.json").read_text())}
            present = max(x["end"] for x in lineages.values())
            for e in rows(out / "episodes.tsv"):
                start, end = float(e["start"]), float(e["end"])
                exposure += min(end, present) - start
                completed += end < present
                outlived += min(end, present) > lineages[e["initiator"]]["end"]
                total += 1
    assert completed > 40
    assert 0.25 < exposure / completed < 0.35
    assert outlived > total / 4


def test_validation():
    for bad in [{"episode_rate": -1}, {"episode_duration": 0}, {"episode_intensity": -0.1}]:
        with pytest.raises(ValueError):
            HorizontalConfig(**bad).validate()
