# __init__.py

"""
ngesh __init__.py
"""

# Version of the ngesh package
__version__ = "0.4.1"
__author__ = "Tiago Tresoldi"
__email__ = "tresoldi@shh.mpg.de"

# Build the namespace
from ngesh.random_tree import gen_tree, add_characters, simulate_bad_sampling
from ngesh.output import tree2nexus, tree2wordlist
from ngesh.utils import set_seeds, random_labels, random_species
