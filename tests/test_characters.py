import statistics
import unittest

from ngesh import RandomStream, Tree
from ngesh.characters import CharacterConfig, evolve
from ngesh.models import control_tree
from ngesh.stats import consistency_index


class CharacterTests(unittest.TestCase):
    def test_poisson_counts(self):
        counts = []
        tree = Tree.from_newick('(A:3,B:3)R:0;')
        for i in range(400):
            result = evolve(tree, CharacterConfig(n_characters=10, rate=.7), seed=i)
            counts.append(len(result.events))
            self.assertTrue(all(0 <= e.time < 3 for e in result.events))
        self.assertAlmostEqual(statistics.mean(counts), 42, delta=1)
        self.assertAlmostEqual(statistics.variance(counts), 42, delta=7)

    def test_lexical_compatibility(self):
        result = evolve(control_tree(20), CharacterConfig(n_characters=30, rate=3), seed=4)
        self.assertEqual(consistency_index(result.truth, result.tree)['ci'], 1)
        self.assertTrue(any(sum(e.recipient == n.name and e.character == c for e in result.events) > 1 for n in result.tree.tips() for c in range(30)))

    def test_invariant_gamma(self):
        cfg = CharacterConfig(n_characters=1000, alpha=.5, categories=4, invariant=.2)
        result = evolve(control_tree(2), cfg, seed=12)
        self.assertGreater(sum(r == 0 for r in result.rates), 150)
        self.assertLessEqual(len(set(result.rates)), 5)

    def test_synonyms(self):
        result = evolve(control_tree(3), CharacterConfig(n_characters=10, rate=8, synonym_duration=.4), seed=7)
        self.assertTrue(any(e.kind == 'synonym_loss' for e in result.events))
        self.assertTrue(any(len(cell) > 1 for row in result.truth.values() for cell in row))

    def test_global_rng_untouched(self):
        import random
        random.seed(19)
        before = random.getstate()
        evolve(control_tree(3), seed=19)
        self.assertEqual(random.getstate(), before)
        self.assertEqual(RandomStream(19).random(), RandomStream(19).random())
