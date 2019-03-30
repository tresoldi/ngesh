# __init__.py

# Version of the ngesh package
__version__ = "0.1"
__author__ = "Tiago Tresoldi"
__email__ = "tresoldi@shh.mpg.de"

# Build the namespace
from ngesh.random_tree import gen_tree, random_labels

# Generates and displays a random tree
def display_random_tree():
    tree = gen_tree(1.0, 0.5, max_time=3, labels="bio")
    tree.show()

