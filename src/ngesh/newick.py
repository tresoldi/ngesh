#!/usr/bin/env python3

"""
newick.py

This module implements a simple method for sorting Newick representation
of phylogenetic trees, allowing to normalize and visually compare the
results of different calls to methods that don't preserve the internal
order (as the branching patterns have no meaning, such as in the
common simile of a phylogenetic tree as a hanging mobile).

The module can be used as a stand-alone tool, called from the
command-line.

The module was written by Tiago Tresoldi (tiago.tresoldi@lingfil.uu.se),
as part of the forthcoming `phylovector` library. The library is developed
in the context of the Cultural Evolution of Texts project at the University
of Uppsala, with funding from the Riksbankens Jubileumsfond (grant
agreement ID: MXM19-1087:1), and released under the MIT license.

If you use this module, please cite the library it is part of:

> Tresoldi, T., (2021). Ngesh: a Python library for synthetic phylogenetic data.
  Journal of Open Source Software, 6(66), 3173, https://doi.org/10.21105/joss.03173
"""

# TODO: branches without length have a "None" string in the output (bug)

# Import Python standard libraries
import re

# Parsing function following the elegant approach by `tricot` at
# https://stackoverflow.com/a/51375562 (licensed under CC BY-SA 4.0)
def _parse(newick: str):
    """
    "Private" function for parsing a Newick string.

    This function is intedend for internal usage by the sorting
    method; if you need an actual parsing of Newick trees, please
    use one of the many stabler libraries.

    :param newick: The Newick tree to be parsed.
    :return: A nested dictionary representation of a Newick tree,
        as provided by the internal `traverse()` function.
    """

    # Find the tokens in the Newick string; we append a semicolon in
    # all cases, as having it makes the parser easier but users
    # might not provide it (or it might be missing in internal
    # calls)
    tokens = re.finditer(
        r"([^:;,()\s]*)(?:\s*:\s*([\d.]+)\s*)?([,);])|(\S)", newick + ";"
    )

    def traverse(next_id: int = 0, parent_id: int = -1):
        """
        Internal function for processing one node.
        """

        this_id = next_id
        children = []

        name, length, delim, char = next(tokens).groups(0)
        if char == "(":
            while char in "(,":
                node, char, next_id = traverse(next_id + 1, this_id)
                children.append(node)
            name, length, delim, char = next(tokens).groups(0)

        return (
            {
                "id": this_id,
                "name": name,
                "length": float(length) if length else None,
                "parent_id": parent_id,
                "children": children,
            },
            delim,
            next_id,
        )

    return traverse()[0]


def sorted_newick(newick):
    """
    Build a sorted representation of a Newick string.

    An internal function parses a Newick string by identifying tokens with
    a regular expression, which might fail for complex trees such as
    those carrying information other than branch length and node name.

    :param newick: The Newick tree to be sorted.

    :return: A corresponding but sorted Newick tree.
    """

    # Parse the tree in newick format
    tree = _parse(newick)

    # Sort the dictionary with the parsed tree
    def _node_sort(node):
        """
        Internal function for performing in-memory sorting.

        Note that, as intended, this will modify in-memory
        the original `node`.
        """
        if node["children"]:
            # First sort children recursively, if any
            for child in node["children"]:
                _node_sort(child)

            # Sort the current node
            node["children"] = sorted(
                node["children"], key=lambda n: (n["name"], -len(n["children"]))
            )

    _node_sort(tree)

    # Build the new textual representation
    def _tree2str(node):
        elms = []
        if not node["children"]:
            return f"{node['name']}:{node['length']}"
        else:
            elms += [_tree2str(n) for n in node["children"]]

        ret = "(%s)" % (",".join(elms))
        if node["length"]:
            ret += f":{node['length']}"

        return ret

    new_newick = _tree2str(tree) + ";"

    # Return
    return new_newick


# Deal with command-line, allowing to use this file as a script
if __name__ == "__main__":
    # Import Python standard libraries used for command-line execution
    import argparse
    import fileinput

    parser = argparse.ArgumentParser(description="Sort a Newick tree.")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        help="Path a file containing the Newick tree to be sorted. If not provided, will default to `stdin`.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Path a file where to write the sorted Newick tree. If not provided, will default to `stdout`.",
    )
    args = parser.parse_args()

    # Read source
    if args.input:
        with open(args.input, encoding="utf-8") as handler:
            newick = handler.read().strip()
    else:
        newick = ""
        for line in fileinput.input():
            newick += line.strip()

    # Convert
    s_newick = sorted_newick(newick)

    # Output
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handler:
            handler.write(s_newick)
    else:
        print(s_newick)
