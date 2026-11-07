import json
from pathlib import Path
import tempfile
import unittest

from ngesh import simulate
from ngesh.__main__ import run_grid
from ngesh.bundle import read_long, read_matrix, verify_bundle
from ngesh.config import toml_subset


class BundleTests(unittest.TestCase):
    def test_roundtrip_and_contract(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp) / "bundle"
            result = simulate(
                seed="roundtrip",
                config={
                    "tree": {"n_taxa": 8},
                    "characters": {"n_characters": 12, "synonym_duration": 2},
                    "observation": {"missing": 0},
                },
                out=dest,
            )
            self.assertEqual(
                read_long(dest / "data.tsv", sorted(result.truth), result.groups), result.truth
            )
            self.assertEqual(read_matrix(dest / "truth.json"), result.truth)
            self.assertTrue(verify_bundle(dest))
            for line in (dest / "data.tsv").read_text().splitlines()[1:]:
                taxon, group, state = line.split("\t")
                self.assertTrue(state.startswith(group + "."))
                self.assertTrue(state.rsplit(".", 1)[1].isdigit())
                self.assertLess(len(line.encode()) + 1, 4096)
            (dest / "data.tsv").write_text("corrupt")
            self.assertFalse(verify_bundle(dest))

    def test_grid_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = {
                "seed": "grid",
                "replicates": 2,
                "base": {"tree": {"n_taxa": 4}, "characters": {"n_characters": 3}},
                "grid": {"characters__rate": [0.1, 0.2]},
            }
            out = Path(tmp)
            self.assertEqual(run_grid(cfg, out), {"completed": 4, "skipped": 0})
            self.assertEqual(run_grid(cfg, out), {"completed": 0, "skipped": 4})
            seeds = [
                json.loads(p.read_text())["master_seed"]["value"] for p in out.glob("*/params.json")
            ]
            self.assertEqual(len(set(seeds)), 4)
            cfg["seed"] = "other"
            with self.assertRaises(ValueError):
                run_grid(cfg, out)

    def test_subset_toml(self):
        text = 'seed="a#b" # comment\n[tree]\nn_taxa=4\n[grid]\ncharacters__rate=[0.1,0.2]\n'
        import sys

        parsed = toml_subset(text)
        if sys.version_info >= (3, 11):
            import tomllib

            self.assertEqual(parsed, tomllib.loads(text))
        self.assertEqual(parsed["seed"], "a#b")

    def test_documented_preset_examples(self):
        from ngesh.config import load_config
        from ngesh.__main__ import main

        with tempfile.TemporaryDirectory() as tmp:
            for preset in ["lexical", "stemma", "genecontent", "morphology"]:
                out = Path(tmp) / preset
                main(
                    [
                        "simulate",
                        "--preset",
                        preset,
                        "--config",
                        f"examples/{preset}.toml",
                        "--seed",
                        "demo",
                        "--out",
                        str(out),
                    ]
                )
                self.assertTrue(verify_bundle(out))
                self.assertIn("tree", load_config(f"examples/{preset}.toml"))

    def test_zero_columns_and_sampled_ancestors(self):
        from ngesh.matrices import read_nexus
        from ngesh import Tree

        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "zero"
            simulate(
                seed=14,
                config={
                    "tree": {"n_taxa": 4, "ancestor_rate": 2},
                    "characters": {"rate": 0},
                    "observation": {"missing": 0, "remove_invariant": True},
                },
                out=out,
            )
            self.assertTrue(all(not row for row in read_nexus(out / "data.multi.nex").values()))
            leaves = Tree.from_newick((out / "truth_sa_leaves.nwk").read_text())
            internal = Tree.from_newick((out / "truth_sa_internal.nwk").read_text())
            leaf_samples = {
                n.name for n in leaves.tips() if n.length == 0 and n.name.startswith("SA")
            }
            internal_samples = {
                n.name for n in internal.traverse() if n.children and n.name.startswith("SA")
            }
            self.assertEqual(leaf_samples, internal_samples)

    def test_readme_doctest(self):
        import doctest

        result = doctest.testfile("../README.md", module_relative=True)
        self.assertEqual(result.failed, 0)

    def test_simulation_golden_and_hashseed(self):
        import os
        import subprocess
        import sys

        fixture = json.loads((Path(__file__).parent / "golden/simulation.json").read_text())
        with tempfile.TemporaryDirectory() as tmp:
            for value in ["1", "9999"]:
                out = Path(tmp) / value
                code = "from ngesh import simulate; import json; simulate(seed='simulation-golden-v1',config=json.loads(__import__('sys').argv[1]),out=__import__('sys').argv[2])"
                subprocess.run(
                    [sys.executable, "-c", code, json.dumps(fixture["config"]), str(out)],
                    check=True,
                    env=dict(os.environ, PYTHONHASHSEED=value, PYTHONPATH="src"),
                )
                self.assertEqual(
                    {name: (out / name).read_text() for name in fixture["outputs"]},
                    fixture["outputs"],
                )

    def test_subset_rejects_invalid_toml(self):
        from ngesh.config import toml_subset

        for text in [
            "x=[1,,2]",
            "x=01",
            "x=1__2",
            "[a..b]\nx=1",
            "x=1\n[x]\ny=2",
            "x=[1,2",
            'x="unterminated',
        ]:
            with self.subTest(text=text), self.assertRaises(ValueError):
                toml_subset(text)

    def test_sampled_ancestor_encodings_preserve_observations(self):
        from ngesh import Tree, simulate

        with tempfile.TemporaryDirectory() as tmp:
            for seed in range(20):
                out = Path(tmp) / str(seed)
                simulate(
                    preset="stemma",
                    seed=seed,
                    config={"tree": {"n_taxa": 10}, "characters": {"n_characters": 5}},
                    out=out,
                )
                truth = Tree.from_newick((out / "truth.nwk").read_text())
                leaves = Tree.from_newick((out / "truth_sa_leaves.nwk").read_text())
                internal = Tree.from_newick((out / "truth_sa_internal.nwk").read_text())
                expected = {n.name for n in truth.tips()}
                self.assertEqual(expected, {n.name for n in leaves.tips()})
                self.assertEqual(
                    expected, {n.name for n in internal.traverse() if n.name in expected}
                )

    def test_mk_nexus_declares_full_symbols_without_tokens(self):
        from ngesh import simulate
        from ngesh.matrices import read_nexus

        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "mk"
            simulate(
                preset="morphology",
                seed=7,
                config={
                    "tree": {"n_taxa": 5},
                    "characters": {"n_characters": 8, "k": 4, "rate": 0},
                },
                out=out,
            )
            text = (out / "data.multi.nex").read_text()
            self.assertNotIn("TOKENS", text)
            self.assertIn('SYMBOLS="0123"', text)
            self.assertEqual(read_nexus(out / "data.multi.nex"), read_matrix(out / "observed.json"))
