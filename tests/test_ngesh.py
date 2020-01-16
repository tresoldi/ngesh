#!/usr/bin/env python3
# encoding: utf-8

"""
test_ngesh
==========

Tests for the `ngesh` package.
"""

# Import third-party libraries
from ete3 import Tree
import hashlib
import logging
import unittest
import sys

# Import the library being tested
import ngesh

# Setup the logger
logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
LOGGER = logging.getLogger("TestLog")

# Some pre-generated newick trees for testing
_TREES = [
    "(((Nucroto zolos:0.339415,Coddopus zoggaus:0.339415)1:2.29301,Aporos "
    "oiasis:2.63243)1:0.511706,Spetitis mubvoppis:3.14413);",
    "((((Ataba eolus:0.274414,Dasoros audus:0.274414)1:1.59309,(Uvuros "
    "spalus:1.63679,(Zilavis sicagas:0.158265,Uzazopus aolo:0.158265)"
    "1:1.47853)1:0.230708)1:0.976915,Spempo gipus:2.84442)1:1.08647,"
    "(((Cobbas linis:0.242355,Ciggus sopebbas:0.242355)1:1.33741,"
    "(Vaoras ovamla:0.235349,Nirceo spemgazzo:0.235349)1:1.34442)"
    "1:0.796904,(Wiopepus spiparzas:1.86067,Eavoros airos:0.906654)"
    "1:0.515998)1:1.55422);",
    "((((Egasis reggasas:0.0747242,(Niponis ecales:0.381168,Vossignis "
    "spepupapes:0.381168)1:0.0385491)1:0.697337,Tuapebbos "
    "eppozas:0.457931)1:2.02569,(Vucotes sparas:3.09293,((Ruttinges "
    "daerus:0.300024,(Ovugas setlanes:0.297199,Rezas "
    "emdis:0.260688)1:1.20787)1:1.35227,(Gemdas naoris:1.77171,(Tirreres "
    "sbissas:0.161549,Speppizipus spebbes:0.33078)1:1.44093)1:1.08563)"
    "1:0.235588)1:0.0498126)1:0.593566,Seti spines:0.874473);",
    "(((Uorus sbepis:1.00768,Winozbus veses:1.00768)1:2.02024,(((((Ragis "
    "zubarivos:0.239467,(Zoprotollas naessadda:0.032829,Sauvo "
    "tibbunnis:0.032829)1:0.206638)1:0.298691,Mirba "
    "ummes:0.538158)1:0.751761,Gitopes guaras:1.28992)1:0.219143,(Uvoros "
    "sbego:0.429929,Rouslis novuddenus:0.429929)1:1.07913)1:0.0145254,"
    "(((Spaupis gieves:0.155583,Ruasus coeddas:0.155583)1:0.222901,Tipapo "
    "vati:0.378484)1:0.044333,Zuvcus luruzuco:0.422817)1:1.10077)"
    "1:1.50433)1:0.97142,((Nunpis ozabmis:0.491355,Sbigirrecas sbilpo:"
    "0.491355)1:0.33275,Agolus galgis:0.824105)1:3.17523);",
]


class TestTree(unittest.TestCase):
    """
    Class for `ngesh` tests related to tree generation.
    """

    # We first test generations, just to see if no exception is thrown, etc.
    # For all tests, we use the same parameters of birth rate as 1.0
    # and death rate as 0.5.
    def test_generation_min_leaves(self):
        """
        Tests tree generation with minimum leaf number stop criterion.
        """

        tree = ngesh.gen_tree(1.0, 0.5, min_leaves=5, seed="myseed")
        assert (
            tree.write()
            == "(L1:2.19578,(((L2:1.62092,(L3:0.15995,L4:0.15995)1:1.46097)1:0.327325,L5:1.94825)1:2.23317,L6:4.18142)1:2.49625);"
        )

    def test_generation_max_time(self):
        """
        Tests tree generation with maximum_time stop criterion.
        """

        tree = ngesh.gen_tree(1.0, 0.5, max_time=2.5, seed="myseed")
        assert (
            tree.write()
            == "(L1:0.473217,(L2:1.81353,(L3:0.57295,L4:0.57295)1:1.24058)1:0.50789);"
        )

    def test_generation_yule_model(self):
        """
        Tests tree generation in a birth-only model.
        """

        tree = ngesh.gen_tree(1.0, 0.0, max_time=2.5, seed="myseed")
        assert tree.write() == "((L1:1.05875,L2:1.05875)1:0.881094,L3:1.93984);"

    def test_generation_labelling(self):
        """
        Tests tree generation with all the label models.
        """

        e_tree = ngesh.gen_tree(
            1.0, 0.5, max_time=0.5, labels="enum", seed="myseed"
        )
        h_tree = ngesh.gen_tree(
            1.0, 0.5, max_time=0.5, labels="human", seed="myseed"
        )
        b_tree = ngesh.gen_tree(
            1.0, 0.5, max_time=0.5, labels="bio", seed="myseed"
        )

        assert (
            e_tree.write()
            == "((L1:0.336704,L2:0.336704)1:0.161799,L3:0.498503);"
        )
        assert (
            h_tree.write()
            == "((Ugeg:0.336704,Hofgid:0.336704)1:0.161799,Vavi:0.498503);"
        )
        assert (
            b_tree.write()
            == "((Opgidas taninnos:0.336704,Uggegus atupuvus:0.336704)1:0.161799,Vavi ragolus:0.498503);"
        )

    def test_generation_pruning(self):
        """
        Tests tree generation with pruning in a birth-death model.
        """

        tree = ngesh.gen_tree(1.0, 0.5, max_time=5.0, prune=True, seed="myseed")
        assert (
            tree.write()
            == "((((((L01:0.424746,L02:0.424746)1:1.18203,L03:1.07416)1:0.0460857,L04:0.887676)1:0.516424,L05:1.6383)1:0.276193,(((L06:0.0121126,L07:0.0121126)1:1.17026,L08:1.18237)1:0.676307,L09:1.61378)1:0.118574)1:0.207203,(L10:0.466091,(L11:0.118302,L12:0.118302)1:0.347789)1:1.08924);"
        )

    def test_generation_polytomy(self):
        """
        Tests tree generation with hard politomy.
        """

        tree = ngesh.gen_tree(1.0, 0.5, min_leaves=25, lam=2.5, seed="myseed")
        assert (
            tree.write()
            == "(L01:2.19578,((L02:0.0682093,L03:0.0682093,L04:0.0682093,L05:0.0682093,L06:0.0682093)1:0.722613,L07:0.790822,(L08:0.691663,(L09:0.293216,L10:0.293216)1:0.398447,L11:0.691663)1:0.099159,L12:0.790822,L13:0.790822,L14:0.790822)1:3.32517,L15:4.116,(((L16:0.0125088,L17:0.0125088,L18:0.0125088,L19:0.0125088,L20:0.0125088,L21:0.0125088)1:0.0482921,L22:0.0608008,L23:0.0608008,L24:0.0608008,L25:0.0608008)1:1.70918,L26:1.76998,L27:1.52598,L28:1.6654,L29:0.74439)1:2.34601);"
        )

    def test_generation_seed_no_label(self):
        """
        Test equality of trees generated with the same seed, no label.
        """

        # No label
        t1 = ngesh.gen_tree(1.0, 0.5, max_time=3.0, labels=None, seed=1234)
        t2 = ngesh.gen_tree(1.0, 0.5, max_time=3.0, labels=None, seed=1234)
        assert t1.write() == t2.write()

    def test_generation_seed_enum_label(self):
        """
        Test equality of trees generated with the same seed, enum label.
        """

        # Enumerating label
        t1 = ngesh.gen_tree(1.0, 0.5, max_time=3.0, labels="enum", seed=1234)
        t2 = ngesh.gen_tree(1.0, 0.5, max_time=3.0, labels="enum", seed=1234)
        assert t1.write() == t2.write()

    def test_generation_seed_human_label(self):
        """
        Test equality of trees generated with the same seed, human label.
        """

        # Enumerating label
        t1 = ngesh.gen_tree(1.0, 0.5, max_time=3.0, labels="human", seed=1234)
        t2 = ngesh.gen_tree(1.0, 0.5, max_time=3.0, labels="human", seed=1234)
        assert t1.write() == t2.write()

    def test_generation_seed_bio_label(self):
        """
        Test equality of trees generated with the same seed, bio label.
        """

        # Enumerating label
        t1 = ngesh.gen_tree(1.0, 0.5, max_time=3.0, labels="bio", seed=1234)
        t2 = ngesh.gen_tree(1.0, 0.5, max_time=3.0, labels="bio", seed=1234)
        assert t1.write() == t2.write()


class TestCharacters(unittest.TestCase):
    """
    Class for `ngesh` tests related to character generation.
    """

    def test_add_characters(self):
        # gamma parameters
        NUM_CONCEPTS = 10
        k = 4.0  # shape
        th = 1.0  # scale
        z = 1.045  # "zipf" correction

        # Add characters to all trees, for coverage
        trees = [
            ngesh.add_characters(
                Tree(newick), NUM_CONCEPTS, k, th, z, seed="myseed"
            )
            for newick in _TREES
        ]

        # Assert the first one
        digest = hashlib.sha256(
            str(ngesh.tree2wordlist(trees[0])).encode("utf-8")
        ).digest()
        assert (
            digest
            == b"\x8a&\x18'\xe9{\x88'\x90U\x1b`\"[\nd\xb0\xbb\x82\x08\xc6$\xb9v\x98\x92\xcc|\xc0\xb7k\xed"
        )


class TestOutput(unittest.TestCase):
    """
    Class for `ngesh` tests related to output generation.
    """

    def test_tree_output(self):
        # Add characters to all test trees
        trees = [
            ngesh.add_characters(
                Tree(newick), 100, 4.0, 1.0, 1.05, seed="myseed"
            )
            for newick in _TREES
        ]

        # Assert the first one
        digest_nx = hashlib.sha256(
            str(ngesh.tree2nexus(trees[0])).encode("utf-8")
        ).digest()
        digest_wl = hashlib.sha256(
            str(ngesh.tree2wordlist(trees[0])).encode("utf-8")
        ).digest()

        assert (
            digest_nx
            == b"\x9d\x81\xd6\xa2/(\x0b\xcb[*7r=9eH\x83\xc7>\x80e\x7f\x93\x93'r\xc5{\ts\x9f\xf3"
        )
        assert (
            digest_wl
            == b".N\x04\xf2a@d~\x90 \x9e2\xc7\x07\xd6\xfb\xc6\xacJ\xad\x1a\xaal%1C\x14\xac\x11\xf6\x86\x99"
        )


if __name__ == "__main__":
    import sys

    sys.exit(unittest.main())
