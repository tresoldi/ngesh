# encoding: utf-8

"""
Random Phylogenetic Tree Generator.

This script provides function to generate random phylogenetic trees in
a Yule (birth only) or Birth-Death model, setting different generation
parameters and limiting the tree in terms of number of leaves and/or
evolution time.
"""

# Import Python standard libraries
import itertools
import math
import random

# Import 3rd party libraries
import numpy as np
from ete3 import Tree

# Import other modules from this library
from ngesh.textgen import *

# Define the maximum number of tries for generation
__MAX_ATTEMPTS = 3000


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


def label_tree(tree, model, seed=None):
    """
    Labels the nodes of a tree according to a model.

    Linguistic labels are unique names generated in a way intended to be
    readable.

    Please note that the `tree` object is changed in place (no return).

    Parameters
    ----------

    tree: ete3 tree object
        The tree whose nodes will be labeled.
    model : str
        A string indicating which model for label generation should be
        used. Possible values are "enum" (for enumerated labels), "human"
        (for random single names), and "bio" (for random biological names).
    seed : value
        An optional seed for the random number generator, only used in case
        of linguistic labels. Defaults to None.
    """

    # Cache the leaves, so we can also obtain their number
    leaves = tree.get_leaves()

    if model == 'bio':
        # As we are using a simple model with replacements, even if
        # extremely unlikely, we might have repeated items in the labels.
        # The execution would not fail as we are using `zip()`, only items
        # would be unnamed, but we are manually adding missing labels as
        # enumerations to make sure there are no anynomous nodes.
        species = list(set(random_species(len(leaves), seed)))
        species += ['L%i' % i for i in range(len(leaves) - len(species))]

        for leaf_node, name in zip(leaves, species):
            leaf_node.name = name

    elif model == 'human':
        for leaf_node, name in zip(leaves, random_labels(len(leaves), seed)):
            leaf_node.name = name

    else:
        # Build the pattern for the label, including the computation of the
        # number of padding zeros needed
        pattern = 'L%%0%ii' % (1 + math.floor(math.log10(len(leaves))))

        # Label all leaves first
        for leaf_idx, leaf_node in enumerate(leaves):
            leaf_node.name = pattern % (leaf_idx+1)


def __gen_tree(birth, death, min_leaves, max_time, labels, lam, prune, seed):
    """
    Internal function for tree generation.

    This is an internal function for the tree generation, whose main
    difference to `gen_tree()`, the one exposed to the user, is that it
    does not guarantee that a tree will be generated, as the parameters and
    the random sampling might lead to dead-ends where all the leaves in
    a tree are extinct before any or all the stopping criteria are met.

    As an internal function, it does not set default values to the arguments
    and does not perform any checking on the values. Information on the
    arguments, which have the same variable names and properties, are given
    in the documentation for `gen_tree()`.
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

    # Create the tree root as a node. Given that the root is at first set as
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

        # Update the total evolution time. If a maximum alloted time
        # `max_time` is provided and we overshoot it, break the loop
        # without implementing the event (as, by the random event time, it
        # would take place *after* our maximum time, in the future).
        total_time += event_time
        if max_time and total_time > max_time:
            break

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
            # Poisson distribution using the `lam` parameter, so that
            # hard politomies are possible). The distance
            # of the children is here initially set to zero, and will be
            # increased by `event_time` in the loop below, along with all
            # other extant nodes.
            for _ in range(2 + np.random.poisson(lam)):
                child_node = Tree()
                child_node.dist = 0
                child_node.extinct = False

                node.add_child(child_node)

        # (Re)Extract the list of extant nodes, now that we might have new
        # children and that the randomly selected node went extinct
        # (easier than directly manipulating the Python list). From the
        # updated list, we will extend the branch length of all extant leaves
        # (thus including any new children) by the `event_time` computed
        # above.
        leaf_nodes = __extant(tree)
        for leaf in leaf_nodes:
            new_leaf_dist = leaf.dist + event_time
            leaf.dist = min(new_leaf_dist, (max_time or new_leaf_dist))

        # If the event above was a death event, we might be in the undesirable
        # situation where all lineages went extinct before we
        # could finish the random generation according to the
        # user-requested parameters, so that one or both stopping criteria
        # cannot be satisfied. A solution could
        # be to recursively call this function, with the same
        # parameters, until a valid tree is found, but this is not
        # optimal (nor elegant) and might get us stuck in a
        # loop if we don't keep track of the number of iterations
        # (especially if we got to this point by using a
        # user-provided random seed). In face of that,
        # it is preferable to be explicit about the problem by
        # returning a None value, with the user (or a wrapper
        # function) being in charge of asserting that the desired
        # number of random trees is collected (even if it is only one).
        if not leaf_nodes:
            tree = None
            break

        # Check whether one or both the stopping criteria were reached
        if min_leaves and len(leaf_nodes) >= min_leaves:
            break

        if max_time and total_time >= max_time:
            break

    # In some cases we might end up with technically valid trees composed
    # only of the root. We make sure at least one speciation event took
    # place, returning `None` as failure in other cases.
    if tree and len(__extant(tree)) <= 2:
        tree = None

    # Prune the tree, removing extinct leaves, if requested and if a
    # tree was found. Remember that the ete3 `prune()` method takes a list
    # of the nodes that will be kept, removing the other ones.
    if prune and tree:
        tree.prune(__extant(tree))

    # Label the tree before returning it, if it was provided
    if labels and tree:
        label_tree(tree, labels, seed=seed)

    return tree


def gen_tree(birth, death, min_leaves=None, max_time=None,
             labels="enum", lam=0.0, prune=False, seed=None):
    """
    Returns a random phylogenetic tree.

    At least one stopping criterion must be informed, with the tree being
    returned when the either is met.

    This function wraps the internal `__gen_tree()` function which cannot
    guarantee that a valid tree will be generated given the user
    parameters and the random sampling. It will try as many times as
    necessary to provide a valid (and reproducible, given a `seed`) tree,
    within the limits of an internal parameter for maximum number of
    attempts.

    Parameters
    ----------

    birth : float
        The birth rate (lambda) for the generated tree.
    death : float
        The death rate (mu) for the generated tree. Must be explicitly set
        to zero for Yule model (i.e., birth only).
    num_leaves : int
        A stopping criterion with the minimum number of extant leaves.
        The generated tree will have at least the number of requested
        extant leaves (possibly more, as the last speciation event might
        produce more leaves than the minimum specified.
        Defaults to None.
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

    # Confirm that at least one stopping condition was provided
    if not (min_leaves or max_time):
        raise ValueError('At least one stopping criterion is required.')

    # Confirm that a valid `labels` was passed
    if labels not in ["enum", "human", "bio", None]:
        raise ValueError("Invalid label model provided ('%s')" % labels)

    # Generate the random tree
    cur_attempt = 0
    while True:
        # Update the current attempt and drop out if maximum was reached.
        # In the extremely unlikely event we were not able to generate a tree,
        # after all these iterations (probably due to the user parameters),
        # we fail explicitly.
        cur_attempt += 1
        if cur_attempt == __MAX_ATTEMPTS:
            raise RuntimeError("Unable to generate a valid tree.")

        # Seed the RNG with the current seed (which, in the first run, will
        # either be the one provided by the user or `None`) and extract a
        # new seed for future iterations if needed (if the provided seed fails,
        # in case of no tree generation, there is no point in trying the
        # feed again).
        random.seed(seed)
        seed = random.random()

        # Ask for a new tree
        tree = __gen_tree(birth, death, min_leaves, max_time, labels,
                          lam, prune, seed)

        # Break out of the loop if a valid tree was found, as in most of the
        # cases; if no tree could be generated, `__gen_tree()` will return
        # `None`.
        if tree:
            break

    return tree
