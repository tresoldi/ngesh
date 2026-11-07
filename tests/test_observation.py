import unittest

from ngesh.characters import CharacterConfig,evolve
from ngesh.models import control_tree
from ngesh.observation import ObservationConfig,observe
from ngesh.rng import RandomStream


class ObservationTests(unittest.TestCase):
    def test_truth_unchanged_and_missing(self):
        result = evolve(control_tree(10),CharacterConfig(n_characters=100),seed=2)
        original = {t:row[:] for t,row in result.truth.items()}
        obs = observe(result,ObservationConfig(missing=.2,coverage_a=3,coverage_b=1,block_missing=.1),rng=RandomStream(2))
        self.assertEqual(original,result.truth)
        self.assertGreater(sum(c is None for row in obs.matrix.values() for c in row),200)
        self.assertTrue(obs.log)

    def test_nested_selection(self):
        result=evolve(control_tree(10),CharacterConfig(n_characters=100,rate=0),seed=4)
        obs=observe(result,ObservationConfig(selection=.5),rng=RandomStream(4))
        self.assertTrue(any(not c for row in obs.matrix.values() for c in row))
        self.assertTrue(all(c is not None for row in obs.matrix.values() for c in row))

    def test_filters_and_taxa(self):
        result=evolve(control_tree(20),CharacterConfig(n_characters=30,rate=1),seed=5)
        obs=observe(result,ObservationConfig(remove_invariant=True,taxon_sampling=.5,lumping=.2,splitting=.1),rng=RandomStream(4))
        self.assertEqual(set(obs.matrix),{n.name for n in obs.tree.tips()})
        self.assertTrue(all(len(row)==len(obs.retained_groups) for row in obs.matrix.values()))
