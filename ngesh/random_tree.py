#!/usr/bin/env python
# encoding: utf-8

"""
Random Phylogenetic Tree Generator.

This script provides function to generate random phylogenetic trees in
a Yule (birth only) or Birth-Death model, setting different generation
parameters and limiting the tree in terms of number of leaves and/or
evolution time.
"""

# Import Python standard libraries
import math
import random

# Import 3rd party libraries
import numpy as np
from ete3 import Tree

def __extant(tree):
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
# TODO: turn into a generator?
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

# TODO: add species model
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

# TODO: note on hard polytomies
def __gen_tree(birth, death, min_leaves, max_time, labels, lam, prune, seed):
    """
    Internal function for tree generation.
    
    This is an internal function for the tree generation, whose main
    difference to `gen_tree()`, the one offered to the user, is that it
    does not guarantee that a tree will be generated, as the parameters and
    the random sampling might lead to dead-ends where all the leaves in
    a tree are extinct before any or all the stopping criteria are met.
    
    As an internal function, it does not set default values to the arguments
    and does not perform any checking on the values. Information on the
    arguments, which have the same names and properties, are given in the
    documentation for `gen_tree()`.
    """

    # Compute the overall event rate (birth plus death), from which the
    # random expovariate will be drawn. `birth` is here normalized in range
    # [0..1] so that we can directly compare with the results of
    # random.random() and decide if the event is a birth or a death.
    # `death` does not need to be normalized, as it is not used anymore (the
    # only check, below, is `random.random() <= birth`).
    event_rate = birth + death
    birth = birth / event_rate

    # Initialize the RNG
    random.seed(seed)

    # Create the tree root as a node. As the root is at first set as
    # non-extinct and with a branch length of 0.0, it will be immediately
    # subject to either a speciation or extinction event. 
    tree = Tree()
    tree.dist = 0.0
    tree.extinct = False

    # Iterate until an acceptable tree is generated (breaking the loop with
    # a tree) or all leaves go extinct (breaking the loop with `tree` as None).
    # `total_time`, of which we keep track in case `max_time` is provided,
    # is the total evolution time (sum of branch lengths) from the root to the
    # extant nodes.
    total_time = 0.0
    while True:
        # Get the list of extant leaves
        leaf_nodes = __extant(tree)

        # Compute the event time before the next birth/death event from a
        # random exporaviate reflecting the number of extant leaves and the
        # combined event probability.
        event_time = random.expovariate(len(leaf_nodes) * event_rate)

        # Update the total evolution time, rescaling `event_time` if needed
        # (cases when we would overshoot the maximum alloted time)
        if max_time and (total_time + event_time) > max_time:
            event_time = max_time - total_time

        total_time += event_time

        # Select a random node among the extant ones and set it as extinct
        # before simulating either a birth or death event; the type of
        # event is decided based on the comparison of the result of a
        # `random.random()` call with `birth` (here already normalized in
        # relation to `event_rate`)
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
        # manipulating the list). From the updated list, we will extend
        # the branch length of all extant leaves (including any new
        # children) by the `event_time` computed above.
        leaf_nodes = __extant(tree) 
        for leaf in leaf_nodes:
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
    # tree was found. Remember that the ete3 `prune()` method takes a list
    # of the nodes that will be kept, removing the other ones.
    if prune and tree:
        tree.prune(__extant(tree))

    # Label the tree before returning it, if it was provided
    # TODO: fix
    if tree:
        label_tree(tree, human=True, seed=seed)

    return tree

# TODO: note on unpossible generation
# TODO: check if `labels` is good.
# TODO: move to lambda (birth) and mu (death); specify they are constant
# TODO: should move from expovariate to poisson?
# TODO: note that there is no guarantee all trees will be different, of course
def gen_tree(birth, death, min_leaves=None, max_time=None,
         labels=False, lam=0.0, prune=False, seed=None):
    """
    Returns a random phylogenetic tree.    

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
    labels : str or None
        The model to be used for generating random labels, either
        "enum" (for enumerated labels), "human" (for random single names),
        "bio" (for random biological names" or None. Defaults to "enum".
    seed : value
        An optional seed for the random number generator. Defaults to None.
    lam : float
        The expectation of interval for sampling a Poisson distribution
        during speciation, with a minimum of two descendants. Should be used
        if more than two descendants are to be allowed. Defaults to zero,
        meaning that all speciation events will have two and only two
        descendents.
    prune : bool
        A flag indicating whether any non-extant leaves should be pruned from
        the tree before it is returned.

    Returns
    -------

    tree : ete3 tree
        The tree randomly generated according to the parameters.
    """
         
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

        tree = __gen_tree(birth, death, min_leaves, max_time,
                    labels, lam, prune, seed)

        # TODO: if notree
        if tree:
            break

    return tree
