"""
test_ngesh
==========

Tests for the `ngesh` package.
"""

# TODO: add test for `fast` method
# TODO: parametrize with pytest

# Import third-party libraries
from ete3 import Tree
import hashlib
import logging
import pytest
import sys

# Import the library being tested
import ngesh

# Setup the logger
logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
LOGGER = logging.getLogger("TestLog")


def test_generation_no_stop():
    """
    Tests an assert is raised if no stopping criteria are given.
    """

    with pytest.raises(ValueError):
        ngesh.gen_tree(1.0, 0.5)


@pytest.mark.parametrize(
    "birth,death,min_leaves,method,seed,expected",
    [
        [
            1.0,
            0.5,
            5,
            "standard",
            "myseed",
            "(L1:0.532642,((L2:0.0290176,L3:0.0290176)1:0.411228,((L4:0.0111507,L5:0.0111507)1:0.00297779,L6:0.0141285)1:0.426117)1:0.131906);",
        ],
        [
            1.0,
            0.5,
            5,
            "fast",
            "myseed",
            "(L1:0.411228,((L2:0.0141285,((L3:0,L4:0)1:0.0111507,L5:0.0111507)1:0.00297779)1:0.0253815,L6:0.03951)1:0.400736);",
        ],
        [
            1.2,
            0.3,
            7,
            "standard",
            12345,
            "(((((L1:1.54313,L2:1.54313)1:0.348404,((L3:0.0243007,L4:0.0243007)1:0.580711,L5:0.605011)1:1.28653)1:0.991873,L6:0.622061)1:0.579099,(L7:1.39093,L8:1.39093)1:2.07157)1:0.822327,L9:0.395632);",
        ],
        [
            1.2,
            0.3,
            7,
            "fast",
            12345,
            "(((((L1:1.39093,(L2:0.605011,L3:0.605011)1:0.785923)1:0.152198,((L4:0,L5:0)1:0.0243007,L6:0.0243007)1:1.51883)1:0.718217,L7:2.26135)1:0.622061,L8:0.991873)1:1.00579,L9:0.426695);",
        ],
    ],
)
def test_generation_min_leaves(birth, death, min_leaves, method, seed, expected):
    """
    Tests tree generation with minimum leaf number stop criterion.
    """

    tree = ngesh.gen_tree(birth, death, min_leaves=min_leaves, method=method, seed=seed)
    assert tree.write() == expected


@pytest.mark.parametrize(
    "birth,death,max_time,method,seed,expected",
    [
        [
            1.0,
            0.5,
            0.1,
            "standard",
            "myseed",
            "((L1:0.0437074,L2:0.0437074)1:0.0503426,L3:0.09405);",
        ],
        [
            1.0,
            0.5,
            0.1,
            "fast",
            "myseed",
            "((L1:0.00594995,L2:0.00594995)1:0.0437074,L3:0.0496574);",
        ],
        [
            1.2,
            0.3,
            0.2,
            "standard",
            12345,
            "(L1:0.0888205,(L2:0.0679173,L3:0.0679173)1:0.0209032);",
        ],
        [
            1.2,
            0.3,
            0.2,
            "fast",
            12345,
            "(L1:0.179097,(L2:0.111179,L3:0.111179)1:0.0679173);",
        ],
    ],
)
def test_generation_max_time(birth, death, max_time, method, seed, expected):
    """
    Tests tree generation with maximum_time stop criterion.
    """

    tree = ngesh.gen_tree(birth, death, max_time=max_time, method=method, seed=seed)
    assert tree.write() == expected


@pytest.mark.parametrize(
    "birth,method,seed,expected",
    [
        [
            1.0,
            "standard",
            "myseed",
            "(L1:0.226171,(L2:0.123796,L3:0.123796)1:0.102375);",
        ],
        [
            1.0,
            "fast",
            "myseed",
            "(L1:0.197625,(L2:0.073829,L3:0.073829)1:0.123796);",
        ],
    ],
)
def test_generation_yule_model(birth, method, seed, expected):
    """
    Tests tree generation in a birth-only model.
    """

    tree = ngesh.gen_tree(birth, 0.0, max_time=0.3, method=method, seed=seed)
    assert tree.write() == expected


def test_generation_labelling():
    """
    Tests tree generation with all the label models.
    """

    e_tree = ngesh.gen_tree(
        1.0, 0.5, max_time=0.5, labels="enum", method="fast", seed="myseed"
    )
    h_tree = ngesh.gen_tree(
        1.0, 0.5, max_time=0.5, labels="human", method="fast", seed="myseed"
    )
    b_tree = ngesh.gen_tree(
        1.0, 0.5, max_time=0.5, labels="bio", method="fast", seed="myseed"
    )

    assert (
        e_tree.write()
        == "(L1:0.477492,((L2:0.0502537,L3:0.0502537)1:0.157858,L4:0.0610461)1:0.26938);"
    )
    assert (
        h_tree.write()
        == "(Hifvepo:0.477492,((Bibeu:0.0502537,Pelbe:0.0502537)1:0.157858,Fuzegpu:0.0610461)1:0.26938);"
    )
    assert (
        b_tree.write()
        == "(Sbibeus neartas:0.477492,((Spelbes rempucis:0.0502537,Spuzegpus spicus:0.0502537)1:0.157858,Wipepo uales:0.0610461)1:0.26938);"
    )

    # Assert error
    with pytest.raises(ValueError):
        ngesh.gen_tree(1.0, 0.5, max_time=0.5, labels="XXX")


def test_generation_pruning():
    """
    Tests tree generation with pruning in a birth-death model.
    """

    tree = ngesh.gen_tree(
        1.0, 0.5, max_time=5.0, prune=True, method="fast", seed="myseed"
    )
    assert (
        tree.write()
        == "((L01:1.4913,(L02:2.46458,((L03:0.696441,L04:0.696441)1:0.502071,((L05:0.0794139,L06:0.0794139)1:0.286782,L07:0.366196)1:0.543533)1:0.0696582)1:0.1741)1:0.329019,(((L08:0.81711,L09:0.533336)1:0.499006,L10:1.31612)1:0.378263,((L11:0.554481,(L12:0.188303,L13:0.188303)1:0.366178)1:0.25932,L14:0.794803)1:0.0667755)1:0.00717659);"
    )


def test_generation_polytomy():
    """
    Tests tree generation with hard politomy.
    """

    tree = ngesh.gen_tree(
        1.0, 0.5, min_leaves=25, lam=2.5, method="fast", seed="myseed"
    )
    assert (
        tree.write()
        == "((L01:0.0123185,(L02:0.33846,L03:0.33846,(L04:0.158692,L05:0.286005,L06:0.286005)1:0.0524546,(L07:0.275961,(L08:0.0651906,L09:0.0651906,L10:0.0651906,L11:0.0651906,L12:0.0651906)1:0.210771,L13:0.275961,L14:0.275961,((L15:0,L16:0,L17:0,L18:0)1:0.0391456,L19:0.0391456,L20:0.0391456)1:0.236816,L21:0.275961,L22:0.275961)1:0.0624986)1:0.0112018,L23:0.0196889,L24:0.349662)1:0.267157,L25:0.272403,((L26:0.323813,L27:0.323813,L28:0.145826)1:0.00715832,L29:0.330971,L30:0.330971,L31:0.161281,L32:0.0761561)1:0.285848);"
    )


def test_generation_seed_no_label():
    """
    Test equality of trees generated with the same seed, no label.
    """

    # No label
    t1 = ngesh.gen_tree(1.0, 0.5, max_time=3.0, labels=None, method="fast", seed=1234)
    t2 = ngesh.gen_tree(1.0, 0.5, max_time=3.0, labels=None, method="fast", seed=1234)
    assert t1.write() == t2.write()


def test_generation_seed_enum_label():
    """
    Test equality of trees generated with the same seed, enum label.
    """

    # Enumerating label
    t1 = ngesh.gen_tree(1.0, 0.5, max_time=3.0, labels="enum", method="fast", seed=1234)
    t2 = ngesh.gen_tree(1.0, 0.5, max_time=3.0, labels="enum", method="fast", seed=1234)
    assert t1.write() == t2.write()


def test_generation_seed_human_label():
    """
    Test equality of trees generated with the same seed, human label.
    """

    # Enumerating label
    t1 = ngesh.gen_tree(
        1.0, 0.5, max_time=3.0, labels="human", method="fast", seed=1234
    )
    t2 = ngesh.gen_tree(
        1.0, 0.5, max_time=3.0, labels="human", method="fast", seed=1234
    )
    assert t1.write() == t2.write()


def test_generation_seed_bio_label():
    """
    Test equality of trees generated with the same seed, bio label.
    """

    # Enumerating label
    t1 = ngesh.gen_tree(1.0, 0.5, max_time=3.0, labels="bio", method="fast", seed=1234)
    t2 = ngesh.gen_tree(1.0, 0.5, max_time=3.0, labels="bio", method="fast", seed=1234)
    assert t1.write() == t2.write()


def test_bad_sampling():
    """
    Test bad sampling simulation on an existing tree.
    """

    tree = ngesh.gen_tree(1.0, 0.666, max_time=10, method="fast", seed="uppsala")
    ngesh.add_characters(
        tree,
        10,
        k=4.0,
        th=1.0,
        mut_exp=1.045,
        k_hgt=2.0,
        th_hgt=1.1,
        seed="uppsala",
    )

    previous = tree.write()
    ngesh.simulate_bad_sampling(tree, 0.5, seed="uppsala")

    digest = hashlib.sha256(str(ngesh.tree2wordlist(tree)).encode("utf-8")).digest()

    assert tree.write() != previous
    assert (
        digest
        == b"eA\x8e\xcf\xa2\xf9\xc1\x0b\xf5xF\x04\xf1\xe7q\xde\x8d\x01\xdd7\xec\x99\xcd\xb0F\xf8\xa5\xfd\xb2\x9f\xc0\xd8"
    )
