ngesh library
=============

Lorem ipsum


# Generates and displays a random tree
def display_random_tree():
    import numpy as np

    # We make sure the death rate is less than the birth one.
    l = np.random.random()
    mu = np.random.random() * l * 0.9

    tree = gen_tree(l, mu, max_time=4.0, labels="human")
    tree.show()
