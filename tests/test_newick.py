"""
test_newick
===========

Tests for the `newick` module.
"""

# Import the library being tested
import ngesh


def test_sorted_newick():
    """
    Tests if sorted_newick() is returning the expected string.
    """

    original_nw = "(((Pulet:0.14,Srufo:0.14):0.24,Mepale:0.39):0.58,Ei:0.98);"
    sorted_nw = "(Ei:0.98,(Mepale:0.39,(Pulet:0.14,Srufo:0.14):0.24):0.58);"

    assert ngesh.sorted_newick(original_nw) == sorted_nw
