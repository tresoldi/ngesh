# __init__.py

# Version of the ngesh package
__version__ = "0.4"
__author__ = "Tiago Tresoldi"
__email__ = "tresoldi@shh.mpg.de"

# Build the namespace
from ngesh.random_tree import gen_tree, add_characters
from ngesh.output import tree2nexus, tree2wordlist
from ngesh.utils import set_seeds, random_labels, random_species

# Generates and displays a random tree
def display_random_tree():
    import numpy as np

    # We make sure the death rate is less than the birth one.
    l = np.random.random()
    mu = np.random.random() * l * 0.9

    tree = gen_tree(l, mu, max_time=4.0, labels="human")
    tree.show()
