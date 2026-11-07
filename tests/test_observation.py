import unittest

from ngesh.characters import CharacterConfig, evolve
from ngesh.models import control_tree
from ngesh.observation import ObservationConfig, observe
from ngesh.rng import RandomStream


class ObservationTests(unittest.TestCase):
    def test_truth_unchanged_and_missing(self):
        result = evolve(control_tree(10), CharacterConfig(n_characters=100), seed=2)
        original = {t: row[:] for t, row in result.truth.items()}
        obs = observe(
            result,
            ObservationConfig(missing=0.2, coverage_a=3, coverage_b=1, block_missing=0.1),
            rng=RandomStream(2),
        )
        self.assertEqual(original, result.truth)
        self.assertGreater(sum(c is None for row in obs.matrix.values() for c in row), 200)
        self.assertTrue(obs.log)

    def test_nested_selection(self):
        result = evolve(control_tree(10), CharacterConfig(n_characters=100, rate=0), seed=4)
        obs = observe(result, ObservationConfig(selection=0.5), rng=RandomStream(4))
        self.assertTrue(any(not c for row in obs.matrix.values() for c in row))
        self.assertTrue(all(c is not None for row in obs.matrix.values() for c in row))

    def test_filters_and_taxa(self):
        result = evolve(control_tree(20), CharacterConfig(n_characters=30, rate=1), seed=5)
        obs = observe(
            result,
            ObservationConfig(
                remove_invariant=True, taxon_sampling=0.5, lumping=0.2, splitting=0.1
            ),
            rng=RandomStream(4),
        )
        self.assertEqual(set(obs.matrix), {n.name for n in obs.tree.tips()})
        self.assertTrue(all(len(row) == len(obs.retained_groups) for row in obs.matrix.values()))

    def test_sampled_selection_is_nested(self):
        from ngesh.models import sample_ancestors

        tree = control_tree(8)
        sample_ancestors(tree, 4, RandomStream(5))
        result = evolve(tree, CharacterConfig(n_characters=100, rate=0), seed=4)
        obs = observe(result, ObservationConfig(selection=0.7), rng=RandomStream(4))
        for node in tree.traverse():
            if not node.attrs.get("sampled_ancestor"):
                continue
            witness = next(c for c in node.children if c.attrs.get("sampled_tip"))
            allowed = {i for i, c in enumerate(obs.matrix[witness.name]) if c}
            stack = [c for c in node.children if not c.attrs.get("sampled_tip")]
            while stack:
                child = stack.pop()
                if child.name in obs.matrix:
                    present = {i for i, c in enumerate(obs.matrix[child.name]) if c}
                    self.assertLessEqual(present, allowed)
                stack.extend(child.children)

    def test_filters_after_taxon_sampling(self):
        from ngesh import simulate

        for seed in range(5):
            result = simulate(
                preset="morphology",
                seed=seed,
                config={"tree": {"n_taxa": 30}, "characters": {"n_characters": 80}},
            )
            obs = observe(
                result,
                ObservationConfig(remove_invariant=True, remove_singleton=True, taxon_sampling=0.3),
                rng=RandomStream(seed),
            )
            for i in range(len(obs.retained_groups)):
                cells = [row[i] for row in obs.matrix.values() if row[i] is not None]
                self.assertGreater(len(set(cells)), 1)
                frequencies = {}
                for cell in cells:
                    for state in cell or {-1}:
                        frequencies[state] = frequencies.get(state, 0) + 1
                self.assertGreaterEqual(sum(v >= 2 for v in frequencies.values()), 2)

    def test_binary_absence_is_a_state(self):
        from ngesh.characters import Evolution, CharacterConfig

        tree = control_tree(4)
        names = [n.name for n in tree.tips()]
        truth = {name: [frozenset({1}) if i < 2 else frozenset()] for i, name in enumerate(names)}
        result = Evolution(
            tree, CharacterConfig(model="binary", n_characters=1), ["C1"], truth, [], [1], {}
        )
        obs = observe(result, ObservationConfig(remove_singleton=True), rng=RandomStream(2))
        self.assertEqual(obs.retained_groups, [0])
