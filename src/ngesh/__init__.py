"""
ngesh __init__.py
"""

# Version of the ngesh package
__version__ = "0.4.1"
__author__ = "Tiago Tresoldi"
__email__ = "tiago.tresoldi@lingfil.uu.se"

# Import from local modules
from .random_tree import gen_tree, add_characters, simulate_bad_sampling
from .output import tree2nexus, tree2wordlist

# Build the namespace
__all__ = [
    "gen_tree",
    "add_characters",
    "simulate_bad_sampling",
    "tree2nexus",
    "tree2wordlist",
]
