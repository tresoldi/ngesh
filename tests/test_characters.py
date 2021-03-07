import hashlib
from ete3 import Tree
import ngesh


def test_add_characters():
    # gamma parameters
    NUM_CONCEPTS = 10
    k = 4.0  # shape
    th = 1.0  # scale
    z = 1.045  # "zipf" correction

    # Create trees, also for coverage/profile
    _trees = [
        ngesh.gen_tree(
            1.0, 0.5, max_time=3.0, labels="bio", method="standard", seed=1234
        ),
        ngesh.gen_tree(1.0, 0.5, min_leaves=25, lam=2.5, method="fast", seed="myseed"),
        ngesh.gen_tree(1.1, 0.0, min_leaves=7, lam=1.5, method="fast", seed="uppsala"),
    ]

    # Add characters to all trees, for coverage
    trees = [
        ngesh.add_characters(tree, NUM_CONCEPTS, k=k, th=th, mut_exp=z, seed="myseed")
        for tree in _trees
    ]

    # Assert the trees by digest
    digest0 = hashlib.sha256(
        str(ngesh.tree2wordlist(trees[0])).encode("utf-8")
    ).digest()
    digest1 = hashlib.sha256(
        str(ngesh.tree2wordlist(trees[1])).encode("utf-8")
    ).digest()
    digest2 = hashlib.sha256(
        str(ngesh.tree2wordlist(trees[2])).encode("utf-8")
    ).digest()

    assert (
        digest0
        == b"9\xbf\x8e\x1a,n\x19\x9c\xb6Z\xa3\x9b\x98\x99\x08\xe6\xe6\x94\xb0]\xa00\xd6s\xe3\x03\xc5\xad\x8a\xfc_\x98"
    )

    assert (
        digest1
        == b"\xd7\xa6'9+dni\xbc\xa6s\xa0\x81\xdf\xcb\x00K$\xf2(\xdc\xf0K\\\x1e/\xaa\x1d\xd7gQ\x18"
    )

    assert (
        digest2
        == b"\xda\xbe;B\x04\xdf\xde\xdd\x88\xd3\x95\xbcL\xca\xe8\xc7'\xde\xd6G@\x03\x16\t\xc1\xe4\x07\x90C\x00M{"
    )


def test_add_characters_with_hgt():

    # Create tree
    tree_hgt = ngesh.gen_tree(
        1.0, 0.5, min_leaves=15, lam=2.5, method="fast", seed="myseed"
    )

    ngesh.add_characters(
        tree_hgt,
        10,
        k=4.0,
        th=1.0,
        mut_exp=1.0,  # z=1.045, TODO: fix test later
        k_hgt=2.0,
        th_hgt=1.1,
        seed="myseed",
    )
    digest = hashlib.sha256(str(ngesh.tree2wordlist(tree_hgt)).encode("utf-8")).digest()

    assert (
        digest
        == b"@\x1c\xf7\xc7\xe1]]\x1e\x10\xd8\xd8<\x8aH\xaa\x1b\xdb\xd6\xfb:\xbbD\xea\x053\x10\x96\x85(@y\xeb"
    )
