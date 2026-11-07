import json
import os
from pathlib import Path
import statistics
import subprocess
import sys
import unittest

from ngesh import RandomStream, Tree, derive_seed
from ngesh.rng import gamma_categories


class CoreTests(unittest.TestCase):
    def test_golden(self):
        r = RandomStream("golden")
        result = {"seed": derive_seed("golden"), "random": [r.random() for _ in range(5)],
                  "gamma": [r.gamma(0.7) for _ in range(4)],
                  "poisson": [r.poisson(35) for _ in range(8)]}
        expected = json.loads((Path(__file__).parent / "golden/core.json").read_text())
        self.assertEqual(result, expected)

    def test_subprocess(self):
        outputs = []
        for value in ["1", "12345"]:
            env = dict(os.environ, PYTHONHASHSEED=value, PYTHONPATH="src")
            outputs.append(subprocess.check_output([sys.executable, "-c",
                "from ngesh import RandomStream; r=RandomStream('x'); print([r.gamma(2) for _ in range(8)])"], env=env))
        self.assertEqual(*outputs)

    def test_seed_types(self):
        self.assertEqual(len({derive_seed(x) for x in ["1", 1, b"1"]}), 3)
        for seed in [None, 1.2, True]:
            with self.assertRaises(TypeError):
                RandomStream(seed)

    def test_moments(self):
        r = RandomStream(521)
        samplers = [(lambda: r.exponential(2), 0.5, 0.25),
                    (lambda: r.gamma(0.7, 2), 1.4, 2.8),
                    (lambda: r.gamma(3, 2), 6, 12),
                    (lambda: r.poisson(35), 35, 35),
                    (lambda: r.beta(2, 3), 0.4, 0.04),
                    (lambda: r.lognormal(0.5), 1, 0.2840254)]
        for sample, mean, variance in samplers:
            draws = [sample() for _ in range(25000)]
            self.assertAlmostEqual(statistics.mean(draws), mean, delta=0.035 * max(1, mean))
            self.assertAlmostEqual(statistics.variance(draws), variance, delta=0.08 * variance)
        for a in [0.2, 1, 5]:
            rates = gamma_categories(a, 8)
            self.assertAlmostEqual(statistics.mean(rates), 1)
            self.assertEqual(sorted(rates), rates)

    def test_tree(self):
        tree = Tree.from_newick("(('a b':1,B:1)X:2,C:3)R:0;")
        a, b, c = tree.tips()
        self.assertEqual(tree.distance(a, c), 6)
        self.assertEqual(tree.distance(a, b), 2)
        self.assertEqual(len(tree.alive(1)), 2)
        self.assertEqual(len(tree.alive(2)), 3)
        self.assertEqual(Tree.from_newick(tree.newick()).newick(), tree.newick())
        pruned = tree.prune({'a b', 'C'})
        self.assertEqual([n.time for n in pruned.tips()], [3, 3])
        pruned.rescale(2)
        self.assertEqual([n.time for n in pruned.tips()], [2, 2])

    def test_network(self):
        tree = Tree.from_newick("((A:1,(B:0.5)X#H1:0.5)D:1,(#H1:0.5,C:1)E:1)R:0;")
        self.assertEqual(len(tree.extra_edges), 1)
        encoded = tree.newick()
        self.assertEqual(Tree.from_newick(encoded).newick(), encoded)

    def test_bad_newick(self):
        for text in ["(A:1,B:2)", "(A:-1,B:2);", "(A:1,B:2;;", "(A:nan,B:1);", "(A:1,B:1);trash"]:
            with self.assertRaises(ValueError):
                Tree.from_newick(text)
