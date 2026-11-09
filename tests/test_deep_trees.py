"""Deep histories must not hit Python's recursion limit (review item 8)."""

import copy
import sys

from ngesh import Tree, simulate
from ngesh.models import TreeConfig, generate_tree
from ngesh.tree import copy_subtree


def caterpillar_newick(n: int) -> str:
    text = "T0:1"
    for i in range(1, n):
        text = f"({text},T{i}:1):1"
    return text + ";"


def test_newick_parser_handles_deep_nesting():
    depth = 5 * sys.getrecursionlimit()
    tree = Tree.from_newick(caterpillar_newick(depth))
    assert len(tree.tips()) == depth
    again = Tree.from_newick(tree.newick())
    assert again.newick() == tree.newick()


def test_newick_parser_keeps_error_messages():
    for text, message in [("((A,B);", "unclosed"), ("(A,B));", "semicolon"), ("(A:x);", "length")]:
        try:
            Tree.from_newick(text)
        except ValueError as error:
            assert message in str(error)
        else:
            raise AssertionError(text)


def test_copy_subtree_matches_deepcopy():
    tree = generate_tree(TreeConfig(n_taxa=30, death=0.5), seed="copy")
    reference = copy.deepcopy(tree.root)
    copied, mapping = copy_subtree(tree.root)
    pairs = [(reference, copied)]
    while pairs:
        a, b = pairs.pop()
        assert (a.name, a.length, a.time, a.attrs) == (b.name, b.length, b.time, b.attrs)
        assert len(a.children) == len(b.children)
        pairs.extend(zip(a.children, b.children))
    assert all(mapping[n].parent is (mapping[n.parent] if n.parent else None) for n in mapping)
    assert copied.attrs is not tree.root.attrs


def test_high_extinction_trees_generate():
    tree = generate_tree(TreeConfig(n_taxa=500, birth=1.0, death=0.9), seed="deep")
    assert len(tree.tips(extant_only=True)) == 500


def test_large_stemma_simulation():
    # Seeds that raised RecursionError in 2.0.0a4 (GSA snapshot via deepcopy).
    for seed in ["heldout-v1:stemma:100:1", "heldout-v1:stemma:100:11"]:
        evolution = simulate(
            preset="stemma",
            seed=seed,
            config={"tree": {"n_taxa": 100}, "characters": {"n_characters": 300}},
        )
        assert evolution is not None
