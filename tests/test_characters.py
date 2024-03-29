"""
test_characters
===============

Tests for character generation in the `ngesh` package.
"""

# Import Python standard libraries
import hashlib

# Import the library being tested
import ngesh


def test_add_characters():
    """
    Test the addition of random characters to trees.
    """

    # gamma parameters
    concepts = 10
    k = 4.0  # shape
    th = 1.0  # scale
    z = 1.045  # "zipf" correction

    # Create trees, also for coverage/profile
    _trees = [
        ngesh.gen_tree(
            1.0, 0.5, max_time=3.0, labels="bio", method="standard", seed=123
        ),
        ngesh.gen_tree(1.0, 0.5, min_leaves=25, lam=2.5, method="fast", seed="myseed"),
        ngesh.gen_tree(1.1, 0.0, min_leaves=7, lam=1.5, method="fast", seed="uppsala"),
    ]

    # Add characters to all trees, for coverage
    trees = [
        ngesh.add_characters(tree, concepts, k=k, th=th, mut_exp=z, seed="myseed")
        for tree in _trees
    ]

    # Assert the trees by digest
    digest0 = hashlib.md5(
        str(ngesh.tree2wordlist(trees[0])).encode("utf-8")
    ).hexdigest()
    digest1 = hashlib.md5(
        str(ngesh.tree2wordlist(trees[1])).encode("utf-8")
    ).hexdigest()
    digest2 = hashlib.md5(
        str(ngesh.tree2wordlist(trees[2])).encode("utf-8")
    ).hexdigest()

    assert digest0 == "e3cdead9e32700032cc0d5adbf0b45b0"
    assert digest1 == "63138f5a7080c2d9accd30e9d151053b"
    assert digest2 == "e451e777e1fc501f5ccc51818e4475a6"


def test_add_characters_with_hgt():
    """
    Test the addition of random characters to trees with HGT.
    """

    # Create tree
    tree_hgt = ngesh.gen_tree(
        1.0, 0.5, min_leaves=15, lam=2.5, method="fast", seed="myseed"
    )

    ngesh.add_characters(
        tree_hgt,
        10,
        k=4.0,
        th=1.0,
        mut_exp=1.045,
        k_hgt=2.0,
        th_hgt=1.1,
        seed="myseed",
    )
    digest = hashlib.md5(str(ngesh.tree2wordlist(tree_hgt)).encode("utf-8")).hexdigest()

    assert digest == "b7333cdb63e0c614f01a5384b7891a4e"
