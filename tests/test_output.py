"""
test_output
===========

Tests for various textual outputs of the `ngesh` package.
"""

# Import Python standard libraries
import hashlib
import pytest

# Import 3rd-party libraries
from ete3 import Tree

# Import the library being tested
import ngesh


@pytest.mark.parametrize(
    "newick,ref_nx_nochar,ref_nx,ref_wl",
    [
        [
            "(((Nucroto zolos:0.339415,Coddopus zoggaus:0.339415)1:2.29301,Aporos "
            "oiasis:2.63243)1:0.511706,Spetitis mubvoppis:3.14413);",
            "e253eb3a1099b49227478f4de3fcdae0",
            "3c5c0fa6e8d63574a699f95ba91aa2c2",
            "aeb95294eff4e2d9fab041b873c5b46e",
        ],
        [
            "((((Ataba eolus:0.274414,Dasoros audus:0.274414)1:1.59309,(Uvuros "
            "spalus:1.63679,(Zilavis sicagas:0.158265,Uzazopus aolo:0.158265)"
            "1:1.47853)1:0.230708)1:0.976915,Spempo gipus:2.84442)1:1.08647,"
            "(((Cobbas linis:0.242355,Ciggus sopebbas:0.242355)1:1.33741,"
            "(Vaoras ovamla:0.235349,Nirceo spemgazzo:0.235349)1:1.34442)"
            "1:0.796904,(Wiopepus spiparzas:1.86067,Eavoros airos:0.906654)"
            "1:0.515998)1:1.55422);",
            "1729eea617c0f4109c84a21b63820c1b",
            "ccf3fa67b5196f4c414d3cf2bf7aab3f",
            "741715e94dd5a03d8e6cbb760d6ee4a0",
        ],
    ],
)
def test_tree_output(newick: str, ref_nx_nochar: str, ref_nx: str, ref_wl: str):
    """
    Test various output methods.
    """

    # Test output without characters
    # Test output for a tree without characters
    tree_nochar = Tree(newick)
    digest_nx_nochar = hashlib.md5(
        str(ngesh.tree2nexus(tree_nochar)).encode("utf-8")
    ).hexdigest()
    assert digest_nx_nochar == ref_nx_nochar

    # Add characters to the tree
    tree = ngesh.add_characters(
        Tree(newick), 100, k=4.0, th=1.0, mut_exp=1.05, seed="myseed"
    )

    # Build nexus and wordlist output for assessment
    digest_nx = hashlib.md5(str(ngesh.tree2nexus(tree)).encode("utf-8")).hexdigest()
    digest_wl = hashlib.md5(str(ngesh.tree2wordlist(tree)).encode("utf-8")).hexdigest()

    assert digest_nx == ref_nx
    assert digest_wl == ref_wl
