import math
from pathlib import Path
import statistics
import tempfile
import unittest

from ngesh import Tree, simulate
from ngesh.characters import CharacterConfig, evolve
from ngesh.horizontal import HorizontalConfig
from ngesh.models import control_tree
from ngesh.stats import consistency_index


class DomainTests(unittest.TestCase):
    def test_dollo_compatibility(self):
        result = evolve(
            control_tree(20),
            CharacterConfig(model="dollo", n_characters=100, rate=1, loss_rate=0.4),
            seed=6,
        )
        # Dollo permits repeated losses, which can make ordinary Mk/Fitch CI < 1.
        # Its defining compatibility criterion is a single origin per trait.
        self.assertTrue(
            all(e.old == (1,) and e.new == () for e in result.events if e.kind == "loss")
        )
        no_loss = evolve(
            control_tree(20),
            CharacterConfig(model="dollo", n_characters=100, rate=1, loss_rate=0),
            seed=6,
        )
        self.assertEqual(consistency_index(no_loss.truth, no_loss.tree, binary=True)["ci"], 1)
        gains = [e.character for e in result.events if e.kind == "gain"]
        self.assertEqual(len(gains), len(set(gains)))

    def test_binary_closed_form(self):
        # root absent, gain .4, loss .8, t=2: pi1*(1-exp(-(gain+loss)*t)).
        cfg = CharacterConfig(
            model="binary",
            n_characters=12000,
            rate=1,
            gain_rate=0.4,
            loss_rate=0.8,
            root_frequency=0,
        )
        result = evolve(Tree.from_newick("(A:2,B:2)R:0;"), cfg, seed=9)
        observed = statistics.mean(bool(c) for c in result.truth["A"])
        self.assertAlmostEqual(observed, (0.4 / 1.2) * (1 - math.exp(-2.4)), delta=0.015)

    def test_mk_closed_form(self):
        cfg = CharacterConfig(model="mk", n_characters=10000, rate=0.7, k=3)
        result = evolve(Tree.from_newick("(A:2,B:2)R:0;"), cfg, seed=2)
        observed = statistics.mean(c == {0} for c in result.truth["A"])
        expected = 1 / 3 + 2 / 3 * math.exp(-3 * 0.7 * 2 / 2)
        self.assertAlmostEqual(observed, expected, delta=0.015)

    def test_asymmetric_stationary(self):
        cfg = CharacterConfig(
            model="mk", n_characters=6000, rate=1, k=2, rate_matrix=[[0, 0.3], [0.7, 0]]
        )
        result = evolve(Tree.from_newick("(A:10,B:10)R:0;"), cfg, seed=2)
        self.assertAlmostEqual(
            statistics.mean(c == {1} for c in result.truth["A"]), 0.3, delta=0.02
        )

    def test_ordered_and_blocks(self):
        cfg = CharacterConfig(
            model="mk",
            n_characters=20,
            rate=3,
            k=4,
            ordered=True,
            block_size=5,
            block_correlation=0.8,
        )
        result = evolve(control_tree(6), cfg, seed=3)
        self.assertTrue(any(e.kind == "correlated_change" for e in result.events))
        self.assertTrue(all(abs(e.old[0] - e.new[0]) == 1 for e in result.events))

    def test_transfer_adds(self):
        result = evolve(
            control_tree(6),
            CharacterConfig(model="binary", n_characters=20, rate=1),
            seed=3,
            horizontal=HorizontalConfig(borrowing_rate=3, mode="add"),
        )
        self.assertTrue(any(e.donor for e in result.events))
        self.assertTrue(all(set(e.old) <= set(e.new) for e in result.events if e.donor))

    def test_all_presets(self):
        with tempfile.TemporaryDirectory() as tmp:
            for preset in ["lexical", "genecontent", "morphology", "stemma"]:
                result = simulate(
                    preset=preset,
                    seed=123,
                    config={"tree": {"n_taxa": 8}, "characters": {"n_characters": 20}},
                    out=Path(tmp) / preset,
                )
                self.assertGreaterEqual(len(result.truth), 8)
                self.assertTrue((Path(tmp) / preset / "data.multi.nex").exists())

    def test_dollo_gain_follows_branch_clock(self):
        from unittest.mock import patch
        from ngesh.rng import RandomStream

        tree = Tree.from_newick("(A:1,B:1)R:0;")
        cfg = CharacterConfig(
            model="dollo", n_characters=2000, rate=1, loss_rate=0, root_frequency=0, clock_sigma=1
        )
        with patch.object(RandomStream, "lognormal", side_effect=[1, 1, 9]):
            result = evolve(tree, cfg, seed=8)
        gains = [e for e in result.events if e.kind == "gain"]
        self.assertAlmostEqual(
            sum(e.recipient == "B" for e in gains) / len(gains), 0.9, delta=0.025
        )

    def test_absent_omissions_are_marked_silent(self):
        result = evolve(
            control_tree(3),
            CharacterConfig(model="textual", n_characters=8, rate=2, omission_rate=10),
            seed=4,
        )
        self.assertTrue(any(e.kind == "omission_silent" for e in result.events))
        self.assertTrue(
            all(not e.old and not e.new for e in result.events if e.kind == "omission_silent")
        )

    def test_sampled_ancestor_shares_clock(self):
        from ngesh.models import sample_ancestors
        from ngesh.rng import RandomStream

        tree = control_tree(5)
        sample_ancestors(tree, 3, RandomStream(3))
        result = evolve(tree, CharacterConfig(n_characters=5, clock_sigma=1), seed=2)
        for node in tree.traverse():
            if node.parent and node.parent.attrs.get("sampled_ancestor"):
                self.assertEqual(
                    result.branch_rates[node.name], result.branch_rates[node.parent.name]
                )
