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
from typing import Hashable, List, Optional

# Import 3rd-party libraries
from ete3 import Tree, TreeNode
import numpy as np

# Import other modules from this library
from . import common

# Define the maximum number of tries for generation
__MAX_ATTEMPTS = 3000


def __extant(tree: Tree) -> List[TreeNode]:
    """
    Internal function returning a list of non-extinct leaves in a tree.

    :param tree: The tree whose nodes will be checked.
    :return: List of extant leaves.
    """

    # Return a filtered list compiled with a list comprehension; the
    # 'extinct' field is not part of ETE3 defaults, but we use here in
    # order to easily differentiate between alive and extinct leaves in
    # Birth-Death models.
    return [leaf for leaf in tree.get_leaves() if leaf.extinct is False]


def label_tree(tree: Tree, model: str = "enum", seed: Optional[Hashable] = None):
    """
    Labels the nodes of a tree according to a model.

    Linguistic labels are unique names generated in a way intended to be
    readable.

    Please note that the `tree` object is changed in place (no return).

    :param tree: The tree whose nodes will be labeled in place.
    :param model: A string indicating which model for label generation should be
        used. Possible values are "enum" (for enumerated labels), "human"
        (for random single names), and "bio" (for random biological names).
    :param seed: An optional seed for the random number generator, only used in case
        of linguistic and biological labels. Defaults to `None`.
    """

    # Cache the leaves, so we can also obtain their number
    leaves = tree.get_leaves()

    if model == "bio":
        # As we are using a simple model with replacements, even if
        # extremely unlikely, we might have repeated items in the labels.
        # The execution would not fail as we are using `zip()`, only items
        # would be unnamed, but we are manually adding missing labels as
        # enumerations to make sure there are no anonymous nodes.
        species = sorted(set(common.random_species(len(leaves), seed)))
        species += ["L%i" % i for i in range(len(leaves) - len(species))]

        for leaf_node, name in zip(leaves, species):
            leaf_node.name = name

    elif model == "human":
        for leaf, name in zip(leaves, common.random_labels(len(leaves), seed)):
            leaf.name = name

    else:
        # Build the pattern for the label, including the computation of the
        # number of padding zeros needed
        pattern = "L%%0%ii" % (1 + math.floor(math.log10(len(leaves))))

        # Label all leaves first
        for leaf_idx, leaf_node in enumerate(leaves):
            leaf_node.name = pattern % (leaf_idx + 1)


def _gen_tree_fast(
    birth: float,
    death: float,
    min_leaves: Optional[int] = None,
    max_time: Optional[float] = None,
    lam: float = 0.0,
    prune: bool = False,
    labels: str = "enum",
    seed: Optional[Hashable] = None,
):
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

    # Initialize the RNG
    common.set_seeds(seed)

    # Compute the overall event rate (birth plus death), from which the
    # random expovariate will be drawn. `birth` is here normalized in range
    # [0..1] so that we can directly compare with the results of
    # `.random()` and decide if the event is a birth or a death.
    # `death` does not need to be normalized, as it is not used anymore (the
    # only check, below, is `.random() <= birth`).
    event_rate = birth + death
    birth = birth / event_rate

    # Create the tree root as a node. Given that the root is at first set as
    # non-extinct and with a branch length of 0.0, it will be immediately
    # subject to either a speciation or extinction event.
    tree = Tree()
    tree.dist = 0.0
    tree.extinct = False

    # list of currently existing leaves.
    # this saves us recurring along the phylogeny each time a sample is collected.
    leaves = [tree]

    # Iterate until an acceptable tree is generated (breaking the loop with
    # a tree) or all leaves go extinct (breaking the loop with `tree` as None).
    # `total_time`, of which we keep track in case `max_time` is provided,
    # is the total evolution time (sum of branch lengths) from the root to the
    # extant nodes.
    total_time = 0.0
    while True:
        # Note that, in comparison to the default method, this implementation
        # does not check for extant leaves, as with the new `leaves` list we
        # can more efficiently keep track of the extant leaf nodes

        # Compute the event time before the next birth/death event from a
        # random expovariate reflecting the number of extant leaves and the
        # combined event probability.
        event_time = random.expovariate(len(leaves) * event_rate)

        # Update the total evolution time. If a maximum allotted time
        # `max_time` is provided and we overshoot it, break the loop
        # without implementing the event (as, by the random event time, it
        # would take place *after* our maximum time, in the future).
        total_time += event_time
        if max_time and total_time > max_time:
            # NICOLA: this is a trick to save time on the branch lengths updates.
            # node.dist is initialized with the time at which the nodes were created,
            # and when a node is terminated or the simulation is finished, node.dist
            # is updated to # represent the branch length.
            for node in leaves:
                node.dist = max_time - node.dist

            break

        # Select a random node among the extant ones and set it as extinct
        # before simulating either a birth or death event; the type of
        # event is decided based on the comparison of the result of a
        # `random.random()` call with `birth` (here already normalized in
        # relation to `event_rate`). In relation to the default implementation,
        # this one keeps track of the position of the parent node
        # within the list leaves.
        leaf_n = np.random.random_integers(len(leaves)) - 1
        node = leaves[leaf_n]
        node.extinct = True

        # As before, here we update node.dist for a node before putting it away
        node.dist = total_time - node.dist
        if np.random.random() <= birth:
            # The event will be a birth (i.e., speciation one), with at least
            # two children (the number is increased by a random sample from a
            # Poisson distribution using the `lam` parameter, so that
            # hard politomies are possible). The distance
            # of the children is here initially set to zero, and will be
            # increased by `event_time` in the loop below, along with all
            # other extant nodes.
            # Here, Nicola's implementation takes care of the first child
            # out of the loop just so it can replace its parent with it in
            # the leaves list - this was not necessary.
            child_node = Tree()

            # As mentioned before, in Nicola's implementation node.dist is
            # initialized as the current time - this is a speed-up trick but
            # the result is the same.
            child_node.dist = total_time
            child_node.extinct = False
            node.add_child(child_node)

            # Child replacing parent in the leaves list
            leaves[leaf_n] = child_node
            for _ in range(1 + np.random.poisson(lam)):
                child_node = Tree()
                child_node.dist = total_time
                child_node.extinct = False

                # Adding child to node and leaves
                node.add_child(child_node)
                leaves.append(child_node)

        else:
            # Remove node from leaves if it dies.
            del leaves[leaf_n]

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
        # user-provided random seed and/or set of unfortunate parameters).
        # In face of that, it is preferable to be explicit about the problem by
        # returning a `None` value, with the user (or a wrapper
        # function) being in charge of asserting that the desired
        # number of random trees is collected (even if it is a single one).
        if not leaves:
            tree = None
            break

        # Check whether the number of leaves stopping criterion was reached
        if min_leaves and len(leaves) >= min_leaves:
            for node in leaves:
                node.dist = total_time - node.dist
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
        tree.prune(leaves)

    # Label the tree before returning it, if it was provided
    if labels and tree:
        label_tree(tree, labels, seed=seed)

    return tree


def _gen_tree(
    birth: float,
    death: float,
    min_leaves: Optional[int] = None,
    max_time: Optional[float] = None,
    lam: float = 0.0,
    prune: bool = False,
    labels: Optional[str] = "enum",
    seed: Optional[Hashable] = None,
) -> Tree:
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

    :param birth: The birth rate (lambda) for the generated tree.
    :param death: The death rate (mu) for the generated tree. Must be explicitly set
        to zero for Yule model (i.e., birth only).
    :param min_leaves: A stopping criterion with the minimum number of extant leaves.
        The generated tree will have at least the number of requested
        extant leaves (possibly more, as the last speciation event might
        produce more leaves than the minimum specified.
        Defaults to `None`.
    :param max_time: A stopping criterion with the maximum allowed time for evolution.
        Defaults to `None`.
    :param lam: The expectation of interval for sampling a Poisson distribution
        during speciation, with a minimum of two descendants. Should be used
        if more than two descendants are to be allowed. Defaults to zero,
        meaning that all speciation events will have two and only two
        descendents.
    :param prune: A flag indicating whether any non-extant leaves should be
        pruned from the tree before it is returned.
    :param labels: The model to be used for generating random labels, either
        "enum" (for enumerated labels), "human" (for random single names),
        "bio" (for random biological names" or None. Defaults to "enum".
    :param seed: An optional seed for the random number generator. Defaults
        to `None`.
    :return: The tree randomly generated according to the parameters.
    """

    # Initialize the RNG
    common.set_seeds(seed)

    # Compute the overall event rate (birth plus death), from which the
    # random expovariate will be drawn. `birth` is here normalized in range
    # [0..1] so that we can directly compare with the results of
    # `.random()` and decide if the event is a birth or a death.
    # `death` does not need to be normalized, as it is not used anymore (the
    # only check, below, is `.random() <= birth`).
    event_rate = birth + death
    birth = birth / event_rate

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
        # random expovariate reflecting the number of extant leaves and the
        # combined event probability.
        event_time = random.expovariate(len(leaf_nodes) * event_rate)

        # Update the total evolution time. If a maximum allotted time
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
        node = np.random.choice(leaf_nodes)
        node.extinct = True
        if np.random.random() <= birth:
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
        # user-provided random seed and/or set of unfortunate parameters).
        # In face of that, it is preferable to be explicit about the problem by
        # returning a `None` value, with the user (or a wrapper
        # function) being in charge of asserting that the desired
        # number of random trees is collected (even if it is a single one).
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


def gen_tree(
    birth: float = 1.0,
    death: float = 0.5,
    method: str = "standard",
    min_leaves: Optional[int] = None,
    max_time: Optional[float] = None,
    lam: float = 0.0,
    prune: bool = False,
    labels: Optional[str] = "enum",
    seed: Optional[Hashable] = None,
) -> Tree:
    """
    Return a random phylogenetic tree.

    At least one stopping criterion must be informed, with the tree being
    returned when such criterion, or either criteria, is/are met.

    This function wraps the internal `__gen_tree()` function which cannot
    guarantee that a valid tree will be generated given the user
    parameters and the random sampling. It will try as many times as
    necessary to provide a valid (and reproducible, given a `seed`) tree,
    within the limits of an internal parameter for maximum number of
    attempts.

    :param birth: The birth rate (lambda) for the generated tree. Defaults
        to 1.0.
    :param death: The death rate (mu) for the generated tree. Must be explicitly set
        to zero for Yule model (i.e., birth only). Defaults to 0.5.
    :param method: The generation method to use. Available methods are "default" and
        "fast" (contributed by Nicola de Maio).
    :param min_leaves: A stopping criterion with the minimum number of extant leaves.
        The generated tree will have at least the number of requested
        extant leaves (possibly more, as the last speciation event might
        produce more leaves than the minimum specified.
        Defaults to `None`.
    :param max_time: A stopping criterion with the maximum allowed time for evolution.
        Defaults to `None`.
    :param lam: The expectation of interval for sampling a Poisson distribution
        during speciation, with a minimum of two descendants. Should be used
        if more than two descendants are to be allowed. Defaults to zero,
        meaning that all speciation events will have two and only two
        descendents.
    :param prune: A flag indicating whether any non-extant leaves should be
        pruned from the tree before it is returned.
    :param labels: The model to be used for generating random labels, either
        "enum" (for enumerated labels), "human" (for random single names),
        "bio" (for random biological names" or None. Defaults to "enum".
    :param seed: An optional seed for the random number generator. Defaults
        to `None`.
    :return: The tree randomly generated according to the parameters.
    """

    # Confirm that at least one stopping condition was provided
    if not (min_leaves or max_time):
        raise ValueError("At least one stopping criterion is required.")

    # Confirm that a valid `labels` was passed
    if labels not in ["enum", "human", "bio", None]:
        raise ValueError("Invalid label model provided ('%s')" % labels)

    # Decide on the method to use
    if method == "standard":
        gen_function = _gen_tree
    elif method == "fast":
        gen_function = _gen_tree_fast
    else:
        raise ValueError("Unknown generation method")

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
        common.set_seeds(seed)
        seed = np.random.random()

        # Ask for a new tree
        tree = gen_function(
            birth=birth,
            death=death,
            min_leaves=min_leaves,
            max_time=max_time,
            labels=labels,
            lam=lam,
            prune=prune,
            seed=seed,
        )

        # Break out of the loop if a valid tree was found, as in most of the
        # cases; if no tree could be generated, `__gen_tree()` will return
        # `None`.
        if tree:
            break

    return tree


def add_characters(
    tree: Tree,
    num_characters: int,
    k: float,
    th: float,
    mut_exp: float = 1.0,
    k_hgt: Optional[float] = None,
    th_hgt: Optional[float] = None,
    seed: Optional[Hashable] = None,
) -> Tree:
    """
    Add random characters to the nodes of a tree.

    Characters are added according to parameters of gamma distributions
    which are related to the length of each branch. The two possible
    events are mutation (assumed to be always to a new character, i.e., no
    parallel evolution) and horizontal gene transfer. No perturbation,
    such as the simulation of errors in sequencing/data collection, is
    performed by this function.

    :param tree: The ete3 tree to which characters will be added. Any previous
        characters will be overridden.
    :param num_characters: The number of characters for each taxa.
    :param k: The k parameter for the gamma distribution related to mutation
        events.
    :param th: The theta parameter for the gamma distribution related to mutation
        events.
    :param k_hgt: The k parameter for the gamma distribution related to horizontal
        gene transfer events. Defaults to None (in case HGT should be
        modelled but the user is unsure about an appropriate value for
        `k_hgt`, it suggested to set it to 1.5 times `k`).
    :param th_hgt: The theta parameter for the gamma distribution related to horizontal
        gene transfer events. Defaults to None (in case HGT should be
        modelled but the user is unsure about an appropriate value for
        `th_hgt`, it suggested to set it to the same value as `th`).
    :param mut_exp: The exponent for correction of mutation probability of each
        character. Defaults to 1.0 (no correction).
    :param seed: An optional seed for the random number generator. Defaults
        to `None`.
    :return: The provided tree, with random characters added.
    """

    # Set seeds
    common.set_seeds(seed)

    # cache a range with the number of characters, for some speed up
    char_range = list(range(num_characters))

    # build the k vector per character, with the correction
    k_vec = [k / (mut_exp ** idx) for idx in char_range]

    # build the k vector per character for horizontal gene transfer
    if k_hgt and th_hgt:
        k_hgt_vec = [k / (mut_exp ** idx) for idx in char_range]
    else:
        k_hgt_vec = None

    # When simulating the character evolution, we need to traverse the tree
    # in chronological order, so that when the characters for each taxon
    # are compiled/simulated, all the characters from chronologically
    # anterior taxa are already available (allowing to model processes such
    # as horizontal gene transfer, etc.). In order to do so, and later
    # facilitating the selection of nodes for some processes, we compile
    # a dictionary of distances from root for all nodes in the tree,
    # sorting them in ascending order for the tree traversal.
    # We also need to sort by distance, with shorter ones first, and then by
    # key, as some nodes (especially the final leaves) will have the exact
    # same distance and we could not otherwise guarantee the order.
    sorted_nodes = sorted(
        tree.traverse("preorder"), key=lambda x: (x.get_distance(tree), x.name)
    )

    # Traverse the tree from the root, adding characters to all
    # nodes (i.e., not only leaves); `states` is the number of states
    # currently used at each point, so we know which number to use in
    # case of a mutation event.
    states = num_characters
    for node in sorted_nodes:
        # Obtain the characters from the immediate ancestor; if there is no
        # immediate ancestor we are at the root, when we default to one
        # different state for each character (copying from the already
        # compiled `char_range` list -- a copy is not strictly needed,
        # but is a best practice here).
        # NOTE: We need the expensive operation of sorting the ancestors
        #       to guarantee reproducibility.
        ancestors = sorted(
            node.get_ancestors(),
            key=lambda x: (x.get_distance(tree), x.name),
            reverse=True,
        )
        if not ancestors:
            node.chars = char_range[:]
            continue

        # In case there are ancestors, we start by making a copy of the
        # states of the immediate ancestor to `chars`; the list can be
        # manipulated during the loop, and will only be attributed to the
        # node at the end.
        chars = ancestors[0].chars[:]

        # We build an array indicating whether each character should be
        # changed (a mutation), based on the comparison of a randomly
        # drawn number of a gamma distribution to the distance of the
        # current node to its immediate ancestor; the same is performed to
        # a borrowing (i.e., horizontal/lateral gene transfer). While
        # different logics could be taken here, this more explicit approach
        # makes debugging and coding of different methods easier. Also
        # note that, as we might have individual `k` parameters due to
        # the exponential correction, we cannot just ask for an array/list
        # of random numbers, but need to iterate one by one.
        mutation_event = [np.random.gamma(k_vec[i], th) < node.dist for i in char_range]

        # Mutate characters according to `mutation_event`.
        for idx, mutation in enumerate(mutation_event):
            if mutation:
                chars[idx] = states
                states += 1

        # Perform HGT simulation if requested
        if k_hgt_vec:
            # For HGT simulation, we first obtain a list of all donors, taxa
            # which can potentially be source for the borrowing
            # (those whose root distance is less to or equal to that of the
            # current taxon, less the current taxon itself), and compute the
            # distance form the current node to each donor (so we can favor
            # closer taxa in the borrowing).
            hgt_event = [
                np.random.gamma(k_hgt_vec[i], th_hgt) < node.dist for i in char_range
            ]

            # NOTE: sorting to guarantee reproducibility
            pot_source = sorted(
                [
                    donor
                    for donor in sorted_nodes
                    if donor.get_distance(tree) < node.get_distance(tree)
                    and donor != node
                ],
                key=lambda x: (x.get_distance(tree), x.name),
            )

            # Compute the probability of each donor, inversely proportional to
            # the current distance (thus the (min+max)-i). For the time being,
            # only a linear distribution of the probabilities is allowed.
            donor_prob = np.array([node.get_distance(donor) for donor in pot_source])
            donor_prob = (min(donor_prob) + max(donor_prob)) - donor_prob

            # Select a random donor for each character, which will be used
            # only if an hgt event is set for the character. While this logic
            # might seem inefficient (we could draw just those characters to
            # be replace), it leaves the code in place for future planned
            # developments.
            # Note that the probability for `np.random.choice()` (used here as
            # `random.choice()` from the standard library is not available in
            # Python 3.5) needs to normalized in range [0,1].
            donor_nodes = np.random.choice(
                pot_source, num_characters, p=(donor_prob / sum(donor_prob))
            )

            # Mutate characters by performing a horizontal gene transfer
            # according to `hgt_event`.
            for idx, hgt, donor_node in zip(char_range, hgt_event, donor_nodes):
                if hgt:
                    chars[idx] = donor_node.chars[idx]

        # Set the new characters.
        node.chars = chars

    return tree


def simulate_bad_sampling(
    tree: Tree, cutoff: Optional[float] = None, seed: Optional[Hashable] = None
):
    """
    Modify a tree in place simulating bad sampling.

    Bad sampling is currently simulated in an uniform distribution, i.e.,
    all existing leaves have the same probability of being removed. Note that
    if a full simulation of tree topology and characters is performed,
    this task must be carried out *after* the simulation of character
    evolution, as otherwise they would fit the sampled tree and not the
    original one.

    As the bad sampling simulation is also based on random numbers,
    while unlikely it is possible that no actual simulation takes place.

    :param tree: ETE3 Tree object for bad sampling simulation.
    :param cutoff: The approximate percentage of extant leaves to remove from the
        tree before returning, simulating uniform bad sampling. As this is
        performed randomly, there is no guarantee that any leaf will
        actually be removed. Default to `None` (no bad sampling simulation).
    :param seed: An optional seed for the random number generator. Defaults to `None`.
    """

    # Initialize the RNG
    common.set_seeds(seed)

    # Simulate bad sampling by building an array of random numbers, the same
    # length of the number of leaves, and dropping thus below a given
    # threshold.
    # NOTE: this is currently only operating on leaves (after pruning, if
    #       requested) and only removing (i.e., not detaching)
    if cutoff:
        rnd_remove_vec = [np.random.random() for _ in tree.get_leaves()]
        for leaf, rnd_remove in zip(tree.get_leaves(), rnd_remove_vec):
            if rnd_remove <= cutoff:
                leaf.delete()
