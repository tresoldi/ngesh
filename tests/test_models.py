import math
import statistics
import unittest

from ngesh.models import (
    TreeConfig,
    ancestor_internal,
    control_tree,
    generate_tree,
    sample_ancestors,
)
from ngesh.rng import RandomStream


class TreeModelTests(unittest.TestCase):
    def test_exact_counts_both_methods(self):
        for method in ["fast", "standard"]:
            for polytomy in [0, 3]:
                tree = generate_tree(
                    TreeConfig(n_taxa=25, death=0.3, polytomy=polytomy, method=method), seed=7
                )
                self.assertEqual(len(tree.tips()), 25)
                times = [n.time for n in tree.tips()]
                self.assertLess(max(times) - min(times), 1e-12)

    def test_yule_height(self):
        n = 15
        heights = [
            max(t.time for t in generate_tree(TreeConfig(n_taxa=n), seed=i).tips())
            for i in range(1800)
        ]
        expected = sum(1 / k for k in range(1, n + 1))
        self.assertAlmostEqual(statistics.mean(heights), expected, delta=0.09)

    def test_sampling(self):
        counts = [
            len(
                generate_tree(TreeConfig(n_taxa=100, shape="balanced", sampling=0.4), seed=i).tips()
            )
            for i in range(400)
        ]
        self.assertAlmostEqual(statistics.mean(counts), 40, delta=1)

    def test_ancestors(self):
        counts = []
        for i in range(500):
            tree = control_tree(10)
            length = sum(n.length for n in tree.traverse())
            counts.append(sample_ancestors(tree, 2, RandomStream(i)))
            encoded = ancestor_internal(tree)
            self.assertEqual(len(encoded.tips()), 10)
            self.assertTrue(all(n.length == 0 for n in tree.tips() if n.attrs.get("sampled_tip")))
        self.assertAlmostEqual(statistics.mean(counts), 2 * length, delta=0.5)

    def test_time_and_rescale(self):
        tree = generate_tree(TreeConfig(n_taxa=None, max_time=3, full=True, death=0.5), seed=7)
        self.assertTrue(all(math.isclose(n.time, 3) for n in tree.tips(True)))
        tree = generate_tree(TreeConfig(n_taxa=10, max_time=10, height=2), seed=3)
        self.assertAlmostEqual(max(n.time for n in tree.tips()), 2)
        for shape in ["balanced", "caterpillar"]:
            tree = control_tree(15, shape)
            self.assertEqual(len(tree.tips()), 15)
            self.assertEqual({n.time for n in tree.tips()}, {1.0})

    def test_limits(self):
        for cfg in [TreeConfig(n_taxa=None), TreeConfig(birth=0), TreeConfig(sampling=-1)]:
            with self.assertRaises(ValueError):
                generate_tree(cfg)

    def test_time_stop_toml_compatible(self):
        tree = generate_tree(TreeConfig(stop="time", max_time=2), seed=91)
        self.assertTrue(all(abs(n.time - 2) < 1e-12 for n in tree.tips()))

    def test_simulate_mean_rescales_observed_tips(self):
        from ngesh import simulate

        result = simulate(
            preset="stemma",
            seed=54,
            config={
                "tree": {"n_taxa": 10, "sampling": 0.6, "mean_height": 2},
                "characters": {"n_characters": 5},
            },
        )
        self.assertAlmostEqual(statistics.mean(n.time for n in result.tree.tips()), 2)
        self.assertTrue(all(n.length >= 0 for n in result.tree.traverse()))

    def test_rescaling_reproducible_across_processes(self):
        import os
        import subprocess
        import sys

        code = "from ngesh import simulate; print(simulate(preset='stemma', seed=54, config={'tree':{'n_taxa':10,'sampling':.6,'mean_height':2},'characters':{'n_characters':5}}).tree.newick())"
        outputs = [
            subprocess.check_output(
                [sys.executable, "-c", code],
                env=dict(os.environ, PYTHONPATH="src", PYTHONHASHSEED=value),
            )
            for value in ["1", "999"]
        ]
        self.assertEqual(outputs[0], outputs[1])

    def test_yule_count_has_positive_terminal_lengths(self):
        for seed in range(100):
            tree = generate_tree(TreeConfig(n_taxa=20), seed=seed)
            self.assertTrue(all(tip.length > 0 for tip in tree.tips()))

    def test_birth_death_count_has_positive_terminal_lengths(self):
        for seed in range(100):
            tree = generate_tree(TreeConfig(n_taxa=20, death=0.5), seed=seed)
            self.assertTrue(all(tip.length > 0 for tip in tree.tips()))

    def test_crown_and_origin_height(self):
        from ngesh import simulate

        for seed in range(20):
            tree = simulate(
                seed=seed, config={"tree": {"n_taxa": 20}, "characters": {"n_characters": 2}}
            ).tree
            tips = tree.tips()
            common = set(tips[0].ancestors())
            for tip in tips[1:]:
                common.intersection_update(tip.ancestors())
            crown = next(n for n in tips[0].ancestors() if n in common)
            self.assertAlmostEqual(max(n.time for n in tips) - crown.time, 1)
            tree = simulate(
                seed=seed,
                config={
                    "tree": {"n_taxa": 20, "origin_height": 2},
                    "characters": {"n_characters": 2},
                },
            ).tree
            self.assertAlmostEqual(max(n.time for n in tree.tips()), 2)
