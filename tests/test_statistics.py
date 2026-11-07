from pathlib import Path
import tempfile
import unittest

from ngesh import Tree, simulate
from ngesh.bundle import read_matrix
from ngesh.matrices import read_nexus, read_phylip
from ngesh.stats import consistency_index, pairwise_distances, quartets


class StatisticsTests(unittest.TestCase):
    def test_distances_missing(self):
        matrix = {
            "a": [frozenset({1}), None, frozenset()],
            "b": [frozenset({2}), frozenset({1}), frozenset()],
        }
        self.assertEqual(pairwise_distances(matrix)[("a", "b")], 0.5)

    def test_quartet_scores(self):
        tree = Tree.from_newick("((A:1,B:1):1,(C:1,D:1):1):0;")
        tips = tree.tips()
        distances = {
            (a.name, b.name): tree.distance(a, b) for i, a in enumerate(tips) for b in tips[i + 1 :]
        }
        self.assertEqual(quartets(distances)["delta"], 0)
        distances["A", "B"] = 5
        distances["C", "D"] = 5
        self.assertEqual(quartets(distances)["delta"], 1)
        self.assertGreater(quartets(distances)["q_residual"], 0)

    def test_polytomy_ci(self):
        tree = Tree.from_newick("(A:1,B:1,C:1,D:1):0;")
        data = {
            "A": [frozenset({0})],
            "B": [frozenset({0})],
            "C": [frozenset({1})],
            "D": [frozenset({1})],
        }
        self.assertEqual(consistency_index(data, tree)["tree_steps"], 2)

    def test_matrix_roundtrips(self):
        with tempfile.TemporaryDirectory() as tmp:
            for preset in ["lexical", "genecontent", "morphology"]:
                out = Path(tmp) / preset
                simulate(
                    preset=preset,
                    seed=4,
                    config={
                        "tree": {"n_taxa": 6},
                        "characters": {"n_characters": 12},
                        "observation": {"missing": 0.2},
                    },
                    out=out,
                )
                expected = read_matrix(out / "observed.json")
                self.assertEqual(
                    read_nexus(out / "data.multi.nex", presence=preset == "genecontent"), expected
                )
                binary = read_nexus(out / "data.bin.nex", presence=True)
                import json

                aliases = json.loads((out / "columns.json").read_text())["phylip_taxa"]
                self.assertEqual(
                    read_phylip(out / "data.bin.phy", names={v: k for k, v in aliases.items()}),
                    binary,
                )

    def test_binary_bitsets_match_cell_distances(self):
        from ngesh import RandomStream

        rng = RandomStream("binary-comparison")
        matrix = {
            str(t): [
                None if (u := rng.random()) < 0.2 else frozenset({1}) if u < 0.6 else frozenset()
                for _ in range(137)
            ]
            for t in range(8)
        }
        for (a, b), actual in pairwise_distances(matrix).items():
            expected = [
                int(x != y) for x, y in zip(matrix[a], matrix[b]) if x is not None and y is not None
            ]
            self.assertAlmostEqual(actual, sum(expected) / len(expected))
