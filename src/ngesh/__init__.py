"""
ngesh __init__.py
"""

# Version of the ngesh package
__version__ = "1.0"
__author__ = "Tiago Tresoldi"
__email__ = "tiago.tresoldi@lingfil.uu.se"

# Import from local modules
from .common import set_seeds
from .random_tree import gen_tree, add_characters, simulate_bad_sampling
from .output import tree2nexus, tree2wordlist

# TODO: move to some "extra" module?
def show_random_tree():
    """
    Shows a random tree using ETE3.

    This function is intended for a quick demonstration of the library.
    """

    # Import the libraries required by the demo, including PyQt5 whose
    # presence is checked here
    import random

    try:
        import PyQt5
    except ImportError:
        raise ImportError("Unable to import `PyQt5` (is it installed?).")

    # We make sure the death rate is less than the birth one.
    l = random.random()
    mu = random.random() * l * 0.9

    tree = gen_tree(l, mu, max_time=4.0, labels="human")
    tree.show()


# Build the namespace
__all__ = [
    "gen_tree",
    "add_characters",
    "simulate_bad_sampling",
    "tree2nexus",
    "tree2wordlist",
    "show_random_tree",
]
