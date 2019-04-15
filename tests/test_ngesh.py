#!/usr/bin/env python3
# encoding: utf-8

"""
test_ngesh
==========

Tests for the `ngesh` package.
"""

import unittest
import ngesh

class TestTextGen(unittest.TestCase):
    """
    Class for `ngesh` tests related to textual generation.
    """

    # We run each generation twice with a seed, checking if it is
    # reproducible
    def test_label_gen(self):
        seq1 = ngesh.textgen.random_labels(size=5, seed=42)
        seq2 = ngesh.textgen.random_labels(size=5, seed=42)
        
        assert tuple(seq1) == tuple(seq2)


    def test_species_gen(self):
        seq1 = ngesh.textgen.random_species(size=5, seed=42)
        seq2 = ngesh.textgen.random_species(size=5, seed=42)
        
        assert tuple(seq1) == tuple(seq2)


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

        ngesh.gen_tree(1.0, 0.5, min_leaves=5)


    def test_generation_max_time(self):
        """
        Tests tree generation with maximum_time stop criterion.
        """
        
        ngesh.gen_tree(1.0, 0.5, max_time=2.5)


    def test_generation_yule_model(self):
        """
        Tests tree generation in a birth-only model.
        """

        ngesh.gen_tree(1.0, 0.0, max_time=2.5)


    def test_generation_labelling(self):
        """
        Tests tree generation with all the label models.
        """
        
        ngesh.gen_tree(1.0, 0.5, max_time=2.5, labels="enum")
        ngesh.gen_tree(1.0, 0.5, max_time=2.5, labels="human")
        ngesh.gen_tree(1.0, 0.5, max_time=2.5, labels="bio")
        ngesh.gen_tree(1.0, 0.5, max_time=2.5, labels=None)


    def test_generation_pruning(self):
        """
        Tests tree generation with pruning in a birth-death model.
        """
        
        ngesh.gen_tree(1.0, 0.5, max_time=5.0, prune=True)


    def test_generation_polytemy(self):
        """
        Tests tree generation with hard politemy.
        """
        
        ngesh.gen_tree(1.0, 0.5, min_leaves=25, lam=2.5)


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


if __name__ == "__main__":
    import sys
    sys.exit(unittest.main())
