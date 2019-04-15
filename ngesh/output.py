# encoding: utf-8

"""
Module with auxiliary function for text generation.
"""

# Import Python standard libraries
import itertools

# TODO: only working with binary data, but this must be changed at
# character level which is the only supported so far
def tree2nexus(tree):
    """
    Returns a string with the representation of a tree in NEXUS format.
    
    Parameters
    ----------
    
    tree : ete3
        The ete3 tree whose NEXUS representation will be returned.
        
    Returns
    -------
    
    buf : string
        A string with the full representation of the tree in NEXUS format.
    """

    # Collect all taxa and their characters
    data = {leaf.name:leaf.chars for leaf in tree.get_leaves()}

    # Collect the number of states used per concept in the entire tree.
    concept_states = [
        set(concept) for concept
        in itertools.zip_longest(*data.values())
    ]    
    
    # Build the textual binary strings
    bin_strings = {}
    for taxon, char in data.items():
        # Build a sequence of booleans indicating whether the state is found
        seq = itertools.chain.from_iterable([
            [concept_state == state
            for state
            in concept_states[concept_idx]]
            for concept_idx, concept_state in enumerate(char)
        ])
        
        # Map the `seq`uence to a binary string
        bin_strings[taxon] = ''.join(['01'[value] for value in seq])

    # Get the length of the longest taxon name for alignment
    # NOTE: This will result in `align_string` being something like "%-24s %s",
    #       guaranteeing left alignment of the taxa name and strings of
    #       characters starting at the same column.
    max_len = max([len(name) for name in data])
    align_string = "%%-%is %%s" % (max_len+3)
    
    # Build the buffer string holding the entire NEXUS file.
    buf = ["#NEXUS", ""]
    buf.append("begin data;")
    buf.append("  dimensions ntax=%i nchar=%i;" % \
        (len(bin_strings), len(list(bin_strings.values())[0])))
    buf.append("  format datatype=standard missing=? gap=-;")
    buf.append("  matrix")
    
    for taxon, bs in bin_strings.items():
        buf.append(align_string % (taxon.replace(' ', '_'), bs))
    
    buf.append("  ;")
    buf.append("end;")
    
    # Join the buffer and return it
    return '\n'.join(buf)
