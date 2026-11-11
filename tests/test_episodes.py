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


def latest_events(path):
    """(lineage, character) -> latest history event at that lineage."""
    out = {}
    for r in rows(path / "history.tsv"):
        key = (r["recipient"], r["character"])
        if key not in out or float(r["time"]) >= float(out[key]["time"]):
            out[key] = r
    return out


def observed(path):
    cells = {}
    for r in rows(path / "data.tsv"):
        cells.setdefault(r["Language"], {}).setdefault(r["Concept"], set()).add(r["Cognateset"])
    return cells


def test_episode_cells_locate_copied_readings():
    with tempfile.TemporaryDirectory() as tmp:
        checked = 0
        for seed in range(4):
            out = run(tmp, f"cells-{seed}", EPISODES, n_taxa=16, n_characters=80)
            lineages = {x["id"]: x for x in json.loads((out / "lineages.json").read_text())}
            episodes = {r["episode"]: r for r in rows(out / "episodes.tsv")}
            cells = rows(out / "episode_cells.tsv")
            latest = latest_events(out)
            events = {
                (r["time"], r["recipient"], r["character"]): r for r in rows(out / "events.tsv")
            }
            data = observed(out)
            carriers = {}
            for row in cells:
                ep, char = row["episode"], row["character"]
                carr = row["carriers"].split(",")
                holders = [t for t in row["holders"].split(",") if t]
                donor_side = [t for t in row["donor_side"].split(",") if t]
                assert carr and not set(carr) & set(holders)
                assert set(donor_side) <= set(holders)
                for t in carr:
                    # the latest event at this character along t's lineage is an effective transfer of ep
                    node, last = t, None
                    while node is not None:
                        e = latest.get((node, char))
                        if e is not None and (
                            last is None or float(e["time"]) > float(last["time"])
                        ):
                            last = e
                        node = lineages[node]["parent"]
                    assert last["kind"] == "contact" and last["old"] != last["new"]
                    assert events[(last["time"], last["recipient"], char)]["episode"] == ep
                    checked += 1
                carriers.setdefault(ep, set()).update(carr)
            for ep, r in episodes.items():
                assert int(r["effective"]) <= int(r["transfers"])
                assert int(r["carriers"]) == len(carriers.get(ep, ()))
            assert all(set(r["holders"].split(",")) - {""} <= set(data) for r in cells)
        assert checked > 0


def test_lost_donor_lines_leave_holders_without_donor_side():
    """In manuscript traditions most exemplars are lost: readings copied from a lost line are
    still shared with surviving witnesses, which donor_side alone would miss."""
    with tempfile.TemporaryDirectory() as tmp:
        holders_only = 0
        for seed in range(3):
            out = Path(tmp) / f"stemma-{seed}"
            simulate(
                preset="stemma",
                seed=f"lost-donor-{seed}",
                config={
                    "tree": {"n_taxa": 30, "ancestor_rate": 0.0},
                    "characters": {"n_characters": 200},
                    "horizontal": {
                        **EPISODES,
                        "episode_rate": 0.6,
                        "episode_intensity": 3.0,
                        "scale": 100.0,
                    },
                },
                out=out,
            )
            for r in rows(out / "episode_cells.tsv"):
                holders_only += bool(r["holders"]) and not r["donor_side"]
        assert holders_only > 0
