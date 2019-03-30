#!/usr/bin/env python3
# encoding: utf-8

"""
test_ngesh
==========

Tests for the `ngesh` package.
"""

import unittest
import ngesh

class TestNgesh(unittest.TestCase):
    """
    Class for `ngesh` tests.
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


if __name__ == "__main__":
    import sys
    sys.exit(unittest.main())
    
#    for i in range(2):
#        yuletree = birth_only(1.0, min_leaves=10, human=True)
#        yuletree = birth_only(1.0, max_time=3)

#    for i in range(2):
#        print("leaves")
#        bdtree = mine(1.0, 0.5, min_leaves=10, human=True)
#        if bdtree:
#            bdtree.show()
            
#        bdtree = gen_tree(1.0, 0.5, max_time=3, human_labels=True)
#        if bdtree:
#            bdtree.show()

    #print(dir(yuletree))
    #print(dir(bdtree))
    #bdtree.show()
