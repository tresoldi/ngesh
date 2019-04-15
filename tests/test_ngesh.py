#!/usr/bin/env python3
# encoding: utf-8

"""
test_ngesh
==========

Tests for the `ngesh` package.
"""

# Import third-party libraries
from ete3 import Tree
import unittest

# Import the library being tested
import ngesh

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


class TestCharacters(unittest.TestCase):
    """
    Class for `ngesh` tests related to character generation.
    """
    
    def test_add_characters(self):
        # gamma parameters
        NUM_CONCEPTS = 10
        k = 4.0 # shape
        th = 1.0 # scale
        z = 1.045 # "zipf" correction

        for newick in _TREES:
            tree = Tree(newick)
            tree = ngesh.add_characters(tree, NUM_CONCEPTS, k, th, z)
            
            #buf = tree2nexus(tree)
            #print(buf)    

if __name__ == "__main__":
    import sys
    sys.exit(unittest.main())
