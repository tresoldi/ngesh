import statistics
import unittest

from ngesh import Tree, simulate
from ngesh.characters import CharacterConfig, evolve
from ngesh.horizontal import HorizontalConfig, displayed_network
from ngesh.models import control_tree


class HorizontalTests(unittest.TestCase):
    def test_zero_identical(self):
        a = evolve(control_tree(12),seed='identical')
        b = evolve(control_tree(12),seed='identical',horizontal=HorizontalConfig(borrowing_rate=0))
        self.assertEqual(a.truth,b.truth)
        self.assertEqual(a.events,b.events)

    def test_contemporaneous_state(self):
        tree = control_tree(12)
        result = evolve(tree,CharacterConfig(n_characters=8,rate=2),seed=9,
                        horizontal=HorizontalConfig(borrowing_rate=3))
        lookup = {n.name:n for n in tree.traverse()}
        transfers = [e for e in result.events if e.donor]
        self.assertGreater(len(transfers),10)
        # Replay all chronological events and inherited states independently.
        cells = {}
        timeline = [(n.parent.time if n.parent else 0,0,i,n) for i,n in enumerate(tree.traverse())]
        timeline += [(e.time,1,i,e) for i,e in enumerate(result.events)]
        for _,kind,_,item in sorted(timeline,key=lambda x:x[:3]):
            if kind == 0:
                cells[item.name] = cells[item.parent.name][:] if item.parent else [frozenset({0})]*8
            else:
                e = item
                if e.donor:
                    donor = lookup[e.donor]
                    self.assertIn(donor,tree.alive(e.time))
                    self.assertEqual(tuple(sorted(cells[e.donor][e.character])),e.new)
                cells[e.recipient][e.character] = frozenset(e.new)
        network = displayed_network(tree,transfers)
        parsed = Tree.from_newick(network.newick())
        self.assertEqual(len(parsed.extra_edges),len(network.extra_edges))

    def test_counts(self):
        counts=[]
        for i in range(180):
            result = evolve(control_tree(2),CharacterConfig(n_characters=5,rate=0),seed=i,
                            horizontal=HorizontalConfig(borrowing_rate=2))
            counts.append(len([e for e in result.events if e.donor]))
        self.assertAlmostEqual(statistics.mean(counts),20,delta=1)

    def test_geography_and_observation(self):
        result = simulate(seed=4,config={'tree':{'n_taxa':8,'shape':'balanced'},
            'horizontal':{'borrowing_rate':1,'distance':'geography'},'observation':{'missing':.2}})
        self.assertTrue(any(e.donor for e in result.events))
