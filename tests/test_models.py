import math
import statistics
import unittest

from ngesh.models import TreeConfig, ancestor_internal, control_tree, generate_tree, sample_ancestors
from ngesh.rng import RandomStream


class TreeModelTests(unittest.TestCase):
    def test_exact_counts_both_methods(self):
        for method in ['fast', 'standard']:
            for polytomy in [0, 3]:
                tree = generate_tree(TreeConfig(n_taxa=25, death=.3, polytomy=polytomy, method=method), seed=7)
                self.assertEqual(len(tree.tips()), 25)
                times = [n.time for n in tree.tips()]
                self.assertLess(max(times) - min(times), 1e-12)

    def test_yule_height(self):
        n = 15
        heights = [max(t.time for t in generate_tree(TreeConfig(n_taxa=n), seed=i).tips()) for i in range(1800)]
        expected = sum(1/k for k in range(1, n))
        self.assertAlmostEqual(statistics.mean(heights), expected, delta=.09)

    def test_sampling(self):
        counts = [len(generate_tree(TreeConfig(n_taxa=100, shape='balanced', sampling=.4), seed=i).tips()) for i in range(400)]
        self.assertAlmostEqual(statistics.mean(counts), 40, delta=1)

    def test_ancestors(self):
        counts = []
        for i in range(500):
            tree = control_tree(10)
            length = sum(n.length for n in tree.traverse())
            counts.append(sample_ancestors(tree, 2, RandomStream(i)))
            encoded = ancestor_internal(tree)
            self.assertEqual(len(encoded.tips()), 10)
            self.assertTrue(all(n.length == 0 for n in tree.tips() if n.attrs.get('sampled_tip')))
        self.assertAlmostEqual(statistics.mean(counts), 2 * length, delta=.5)

    def test_time_and_rescale(self):
        tree = generate_tree(TreeConfig(n_taxa=None, max_time=3, full=True, death=.5), seed=7)
        self.assertTrue(all(math.isclose(n.time, 3) for n in tree.tips(True)))
        tree = generate_tree(TreeConfig(n_taxa=10, max_time=10, height=2), seed=3)
        self.assertAlmostEqual(max(n.time for n in tree.tips()), 2)
        for shape in ['balanced', 'caterpillar']:
            tree = control_tree(15, shape)
            self.assertEqual(len(tree.tips()), 15)
            self.assertEqual({n.time for n in tree.tips()}, {1.0})

    def test_limits(self):
        for cfg in [TreeConfig(n_taxa=None), TreeConfig(birth=0), TreeConfig(sampling=-1)]:
            with self.assertRaises(ValueError):
                generate_tree(cfg)
