"""
Module with auxiliary function for output generation.
"""

# Import Python standard libraries
import itertools

# Import 3rd-party libraries
from ete3 import Tree


def tree2wordlist(tree: Tree) -> str:
    """
    Returns a string with the representation of a tree in wordlist format.

    :param tree: The ete3 tree whose CSV representation will be returned.
    :return: A string with the full representation of the tree in CSV format.
    """

    # The number of characters for each taxon; it will be set when the
    # first taxon data is read.
    num_chars = None

    # Buffer for the data to be returned
    buf = ["Language_ID,Feature_ID,Value"]

    # Iterate over all leaves, collecting data
    rows = []
    for leave in tree.get_leaves():
        # Get the number of characters if necessary
        if not num_chars:
            num_chars = len(leave.chars)

        # Iterate over all characters of the current leave
        rows += [
            [leave.name, "feature_%i" % idx, str(leave.chars[idx])]
            for idx in range(num_chars)
        ]

    # Sort the row info, so it better resembles real datasets, and
    # add it to the buffer
    rows = sorted(rows, key=lambda r: (r[1], r[0]))
    for row in rows:
        buf.append(",".join(row))

    # Join the buffer and return it
    return "\n".join(buf)


def tree2nexus(tree: Tree) -> str:
    """
    Returns a string with the representation of a tree in NEXUS format.

    :param tree: The ete3 tree whose NEXUS representation will be returned.
    :return: A string with the full representation of the tree in NEXUS format.
    """

    # Collect all taxa and their characters, provided the characters
    # exist
    try:
        data = {leaf.name: leaf.chars for leaf in tree.get_leaves()}
        missing_chars = False
    except AttributeError:
        data = {leaf.name: [] for leaf in tree.get_leaves()}
        missing_chars = True

    # Collect the number of states used per concept in the entire tree.
    concept_states = [
        set(concept)
        for concept in itertools.zip_longest(*[data[taxon] for taxon in sorted(data)])
    ]

    # Build the textual binary strings
    bin_strings = {}
    for taxon in sorted(data):
        # Build a sequence of booleans indicating whether the state is found
        seq = itertools.chain.from_iterable(
            [
                [concept_state == state for state in concept_states[concept_idx]]
                for concept_idx, concept_state in enumerate(data[taxon])
            ]
        )

        # Map the `seq`uence to a binary string
        bin_strings[taxon] = "".join(["01"[value] for value in seq])

    # Get the length of the longest taxon name for alignment
    # NOTE: This will result in `align_string` being something like "%-24s %s",
    #       guaranteeing left alignment of the taxa name and strings of
    #       characters starting at the same column.
    max_len = max([len(name) for name in data])
    align_string = "%%-%is %%s" % (max_len + 3)

    # Build the buffer string holding the entire NEXUS file.
    buf = ["#NEXUS", ""]

    if missing_chars:
        buf.append("[WARNING: characters missing from tree]\n")

    buf.append("begin data;")
    buf.append(
        "  dimensions ntax=%i nchar=%i;"
        % (len(bin_strings), len(list(bin_strings.values())[0]))
    )
    buf.append("  format datatype=standard missing=? gap=-;")
    buf.append("  matrix")

    for taxon in sorted(bin_strings):
        buf.append(align_string % (taxon.replace(" ", "_"), bin_strings[taxon]))

    buf.append("  ;")
    buf.append("end;")

    # Join the buffer and return it
    return "\n".join(buf)
