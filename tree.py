#!/usr/bin/env python
# encoding: utf-8

"""
Random Phylogenetic Tree Generator.

This script provides function to generate random phylogenetic trees in
a Yule (birth only) or Birth-Death model, setting different generation
parameters and limiting the tree in terms of number of leaves and/or
evolution time.

The implementation is based on the code and instructions originally
provided by `mrnouthahi` at
https://mrnoutahi.com/2017/12/05/How-to-simulate-a-tree/.
"""

# Import Python standard libraries
import math
import random

# Import 3rd party libraries
import numpy as np
from ete3 import Tree

def _extant(tree):
    """
    Internal function returning a list of non-extinct leaves in a tree.
    """

    # Return a filtered list compiled with a list comprehension; the
    # 'extinct' field is not part of ETE3 defaults, but we use here in
    # order to easily differentiate between alive and extinct leaves in
    # Birth-Death models.
    return [
        leave for leave in tree.get_leaves()
        if leave.extinct is False
    ]


# TODO: there is an infinitesimal possibility that this might get stuck in
# a loop, fix it
# TODO: remove repeated, maybe also with semi-vowels
# TODO: /h/ only at the beginning of after vowels (i.e not after consoannts)
def random_labels(size=1, seed=None):
    """
    Returns a list of unique random pronounceable labels.

    Parameters
    ----------
    size : int
        The number of labels in the returned set. Defaults to one.
    seed : value
        An optional seed for the random number generator. Defaults to None.

    Returns
    -------
    labels : list of strings
        The list of unique labels.
    """

    # Set up the basic data
    sounds = {
        'C' : [c for c in 'bpdtfvszrlgkmnh'],
        'V' : [v for v in 'aeiou'],
    }

    patterns = ['V', 'CV', 'CV', 'CVC']

    # Initialize the RNG
    random.seed(seed)

    # Collect data until enough labels have been found
    ret = []
    while True:
        # draw a random number of syllables
        syllables = []

        # Map each syllable to random sounds
        for _ in range(random.randint(2, 4)):
            syllable = ''.join([
                random.choice(sounds[sound_class])
                for sound_class in random.choice(patterns)
            ])

            syllables.append(syllable)

        # Build the string and add it if not used yet
        label = ''.join(syllables).capitalize()
        if label not in ret:
            ret.append(label)

        # Check if we have enough data
        if len(ret) == size:
            break

    # Return the list of labels
    return ret

# TODO: add species model?
def label_tree(tree, prefix='L', human=False, seed=None):
    """
    Labels the nodes of a tree in an enumerating or linguistic way.

    Linguistic labels are unique names generated in a way intended to be
    readable. Enumerating labels are, as expected, a plain enumeration.
    The `tree` object is changed in place (no return).

    Parameters
    ----------

    prefix: str
        The prefix string in case of enumerating labels. Defaults to "L".
    human : bool
        Whether to use or not random linguistic labels. Defaults to False.
    seed : value
        An optional seed for the random number generator, only used in case
        of linguistic labels. Defaults to None.
    """

    # Cache the leaves, so we can also obtain their number
    leaves = tree.get_leaves()

    # Add enumerating or linguistic labels
    if human:
        for leaf_node, name in zip(leaves, random_labels(len(leaves), seed)):
            leaf_node.name = name
    else:
        # Build the pattern for the label, including the computation of the
        # number of padding zeros needed
        pattern = '%s%%0%ii' % (
            prefix,
            1 + math.floor(math.log10(len(leaves))))

        # Label all leaves first
        for leaf_idx, leaf_node in enumerate(leaves):
            leaf_node.name = pattern % (leaf_idx+1)


# TODO: add prune param
def gen_tree(birth, death, min_leaves=None, max_time=None,
         human_labels=False, seed=None, lam=0, prune=False):
    """
    Returns a random Birth-Death tree.    

    At least one stopping criterion must be informed, with the tree being
    returned when the either is met.

    Parameters
    ----------

    birth : float
        The birth rate for the generated tree.
    death : float
        The death rate for the generated tree. Must be explicitly set to zero
        for Yule model (i.e., birth only).
    num_leaves : int
        A stopping criterion with the desired number of leaves. Defaults to
        None.
    max_time : float
        A stopping criterion with the maximum allowed time for evolution.
        Defaults to None.
    human_labels : bool
        Whether to use or not random linguistic labels. Defaults to False.
    seed : value
        An optional seed for the random number generator. Defaults to None.
    lam : float
        The expectation of interval for sampling a Poisson distribution
        during speciation, with a minimum of two descendants. Should be used
        if more than two descendants are to be allowed. Defaults to zero,
        meaning that all speciation events will have two and only two
        descendents.

    Returns
    -------

    tree : ete3 tree
        The tree randomly generated according to the parameters.
    """

    # Confirm that at least one stopping condition was provided
    if not (min_leaves or max_time):
        raise ValueError('At least one stopping criterion is required.')

    # Compute the overall event rate (birth plus death), from which the
    # random expovariate will be drawn. `birth` is here normalized in range
    # [0..1] so that we can directly compare with the results of
    # random.random() and decide if the event is a birth or a death.
    # `death` does not need to be normalized, as it is not used anymore. 
    event_rate = birth + death
    birth = birth / event_rate

    # Initialize the RNG
    random.seed(seed)

    # Create the tree root as a node and sets its branch length to 0.0.
    # The root is at first set to non-extinct with our custom "extinct"
    # feature.
    tree = Tree()
    tree.dist = 0.0
    tree.extinct = False

    # Repeat until an acceptable tree is generated; `total_time`, of which
    # we keep track in case `max_time` is provided, is the total time of
    # evolution that, in this case of a Yule tree, is the
    # uniform distance of all leaves at any given moment.
    total_time = 0.0
    while True:
        # Get the list of extant species/languages
        leaf_nodes = _extant(tree)

        # Compute a random expovariate event time before the next speciation
        # from the number of leaves and the provided birth rate
        event_time = random.expovariate(len(leaf_nodes) / birth)

        # Update the total evolution time, rescaling `event_time` if needed
        if max_time and (total_time + event_time) > max_time:
            event_time = max_time - total_time

        total_time += event_time

        # Draw a random node among the extant ones, and set it as extinct,
        # optionally (depending of the value of random.random(), which is
        # compared to `birth` as already normalized in relation to
        # `event_rate`.
        node = random.choice(leaf_nodes)
        node.extinct = True
        if random.random() <= birth:
            # The event will be a birth (i.e., speciation one), with at least
            # two children (the number is increased by a random sample from a
            # Poisson distribution using the `lam` parameter). The distance
            # of the children is here first set to zero, and will be
            # increased by `event_time` in the loop below, along with all
            # other extant nodes.
            for _ in range(2 + np.random.poisson(lam)):
                child_node = Tree()
                child_node.dist = 0
                child_node.extinct = False
                
                node.add_child(child_node)

        # Extract the list of extant nodes now that we might have new children
        # and that the randomly selected node went extinct (easier than
        # manipulating the list). From the updated list, we will update
        # ("extend") the branch length of all extant leaves by the
        # `event_time` computed above.
        leaf_nodes = _extant(tree) 
        for leaf in leaf_nodes:
            if leaf != node:
                new_leaf_dist = leaf.dist + event_time
                leaf.dist = min(new_leaf_dist, (max_time or new_leaf_dist))

        # If the event above was a death event, we might be in the undesirable
        # situation where all lineages went extinct before we
        # could finish the random generation according to the
        # user-requested parameters. A solution could
        # be to recursively call the function, with the same
        # parameters, until a valid tree is found, but this is not
        # optimal (nor elegant) and might get us stuck in a
        # loop if we don't keep track of the number of iterations
        # (especially if we got to this point by using a
        # user-provided random seed). In face of that,
        # it is preferable to be explicit about the problem by
        # returning a None value, with the user (or a wrapper
        # function) in charge of taking care that the desired
        # number of random trees is collected (even if it is only one).
        if len(leaf_nodes) == 0:
            tree = None
            break

        # Check whether one or both the stopping criteria were reached
        if min_leaves and len(leaf_nodes) >= min_leaves:
            break

        if max_time and total_time >= max_time:
            break

    # Prune the tree, removing extinct leaves, if requested and if a
    # tree was found
    if prune and tree:
        tree.prune(_extant(tree))

    # Label the tree before returning it, if it was provided
    if tree:
        label_tree(tree, human=human_labels, seed=seed)

    return tree


def gen_tree_safe(birth, death, min_leaves=None, max_time=None,
         human_labels=False, seed=None, lam=0, prune=False):
    max_attempts = 10000

    cur_attempt = 0
    tree = None
    while True:
        # Update the current attempt and drop out if maximum was reached
        # In the extremely unlikely event we were not able to generate a tree, let's
        # fail explicitly
        cur_attempt += 1
        if cur_attempt == max_attempts:
            raise RuntimeError("Unable to generate a valid birth-death tree.")

        # Generate a new seed -- if a seed was provided and it fails (no
        # tree generation), there is no point in trying with the same number;
        # this way, we feed a different but predictable seed each time
        random.seed(seed)
        seed = random.random()

        tree = gen_tree(birth, death, min_leaves, max_time,
                    human_labels, seed, lam, prune)

        if tree:
            break

    return tree



if __name__ == '__main__':
#    for i in range(2):
#        yuletree = birth_only(1.0, min_leaves=10, human=True)
#        yuletree = birth_only(1.0, max_time=3)

    for i in range(2):
#        print("leaves")
#        bdtree = mine(1.0, 0.5, min_leaves=10, human=True)
#        if bdtree:
#            bdtree.show()
            
        bdtree = gen_tree_safe(1.0, 0.5, max_time=3, human_labels=True)
        if bdtree:
            bdtree.show()

    #print(dir(yuletree))
    #print(dir(bdtree))
    #bdtree.show()
