"""
Common data structures and functions.
"""

import ete3


def sorted_newick(newick):
    """
    Build a sorted representation of a Newick string.

    An internal function parses a Newick string by identifying tokens with
    a regular expression, which might fail for complex trees such as
    those carrying information other than branch length and node name.

    :param newick: The Newick tree to be sorted.

    :return: A corresponding but sorted Newick tree.
    """

    tree = ete3.Tree(newick)
    tree.sort_descendants()
    tree.ladderize()

    return tree.write(format=1)


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
