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
            == "(L1:0.532642,((L2:0.0290176,L3:0.0290176)1:0.411228,((L4:0.0111507,L5:0.0111507)1:0.00297779,L6:0.0141285)1:0.426117)1:0.131906);"
        )

    def test_generation_max_time(self):
        """
        Tests tree generation with maximum_time stop criterion.
        """

        tree = ngesh.gen_tree(1.0, 0.5, max_time=1.0, seed="myseed")
        assert (
            tree.write()
            == "(L1:0.532642,(((L2:0.183099,L3:0.183099)1:0.0312144,L4:0.0414806)1:0.411228,(((L5:0.0772764,L6:0.185296)1:0.0111507,L7:0.111863)1:0.00297779,L8:0.199425)1:0.426117)1:0.131906);"
        )

    def test_generation_yule_model(self):
        """
        Tests tree generation in a birth-only model.
        """

        tree = ngesh.gen_tree(1.0, 0.0, max_time=1.0, seed="myseed")
        assert (
            tree.write()
            == "((((L01:0.0210924,L02:0.0210924)1:0.114516,((L03:0.108258,L04:0.108258)1:0.0135044,(L05:0.110611,L06:0.110611)1:0.0111507)1:0.0138468)1:0.0157386,(L07:0.124442,L08:0.124442)1:0.0269054)1:0.798963,((L09:0.096708,L10:0.096708)1:0.655743,L11:0.752451)1:0.197859);"
        )

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
            == "(L1:0.449746,((L2:0.0968117,L3:0.0968117)1:0.330426,L4:0.26938)1:0.0225083);"
        )
        assert (
            h_tree.write()
            == "(Hifvepo:0.449746,((Bibeu:0.0968117,Pelbe:0.0968117)1:0.330426,Fuzegpu:0.26938)1:0.0225083);"
        )
        assert (
            b_tree.write()
            == "(Sbibeus neartas:0.449746,((Spelbes rempucis:0.0968117,Spuzegpus spicus:0.0968117)1:0.330426,Wipepo uales:0.26938)1:0.0225083);"
        )

    def test_generation_pruning(self):
        """
        Tests tree generation with pruning in a birth-death model.
        """

        tree = ngesh.gen_tree(1.0, 0.5, max_time=5.0, prune=True, seed="myseed")
        assert (
            tree.write()
            == "(L01:2.08229,((((L02:0.836492,(L03:0.481244,L04:0.481244)1:0.266392)1:0.256375,L05:1.09287)1:0.325198,(((L06:0.0445452,L07:0.0445452)1:0.264622,L08:0.309167)1:0.323491,L09:0.632659)1:0.297593)1:1.46117,(((L10:0.925725,(L11:0.743874,L12:0.18682)1:0.0634664)1:0.304216,L13:0.490981)1:0.0129382,L14:0.740564)1:0.385008)1:0.692202);"
        )

    def test_generation_polytomy(self):
        """
        Tests tree generation with hard politomy.
        """

        tree = ngesh.gen_tree(1.0, 0.5, min_leaves=25, lam=2.5, seed="myseed")
        assert (
            tree.write()
            == "((L01:0.278359,L02:0.267157,((L03:0.329973,L04:0.329973,L05:0.329973)1:0.00737038,L06:0.337343,L07:0.00637183,(L08:0.154122,(L09:0.286005,L10:0.286005,L11:0.286005,L12:0.108018,L13:0.010044,L14:0.286005,L15:0.03119)1:0.0378078,L16:0.323813)1:0.0135302,L17:0.337343)1:0.279476,((L18:0.127314,L19:0.127314,L20:0.127314,((L21:0.0391456,L22:0.0391456,L23:0.0391456,L24:0.0391456)1:0.026045,L25:0.0651906,L26:0.0651906)1:0.0621232,L27:0.127314)1:0.217102,L28:0.344416,L29:0.344416,L30:0.344416)1:0.272403)1:0.131906,L31:0.748725,L32:0.748725);"
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
                Tree(newick), NUM_CONCEPTS, k=k, th=th, z=z, seed="myseed"
            )
            for newick in _TREES
        ]

        # Assert the first one
        digest = hashlib.sha256(
            str(ngesh.tree2wordlist(trees[0])).encode("utf-8")
        ).digest()
        assert (
            digest
            == b"\xdc<\x1f\x10N\xbf\xcc\xc4|l26\x10\xfc\xbaN\xb7\\c\x8bB\xca\x95.\xcbH\x82T\xa3\xbd\xff\x15"
        )


class TestOutput(unittest.TestCase):
    """
    Class for `ngesh` tests related to output generation.
    """

    def test_tree_output(self):
        # Add characters to all test trees
        trees = [
            ngesh.add_characters(
                Tree(newick), 100, k=4.0, th=1.0, e=1.05, seed="myseed"
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
            == b"E\xf8\x97\xb6*\x7f\xf4_j\x89\x02dn\x1d\xbe\xb0\xb6\xcd\xd9.\xca:\x9ft\xe2m\xc5y\xa5\xaa\x0fa"
        )
        assert (
            digest_wl
            == b"\xb6\xc8i\xf4!\xf4l\x91\xb9\x8d\xb5Kae\x1aF\x9c\xfd\n \x06\xf2D\x1e<\xdd(U]6(\xbf"
        )


if __name__ == "__main__":
    import sys

    sys.exit(unittest.main())
