# __init__.py

# Version of the ngesh package
__version__ = "0.3.1"
__author__ = "Tiago Tresoldi"
__email__ = "tresoldi@shh.mpg.de"

import random

# Build the namespace
from ngesh.random_tree import gen_tree, add_characters
from ngesh.output import tree2nexus, tree2wordlist

# Generates and displays a random tree
def display_random_tree():
    # We make sure the death rate is less than the birth one.
    l = random.random()
    mu = random.random() * l * 0.9

    tree = gen_tree(l, mu, max_time=4.0, labels="human")
    tree.show()
