#!/usr/bin/env python3
# encoding: utf-8

"""
__main__.py

Module for command-line execution and generation of random trees.
"""

# Import Python standard libraries
import argparse
import configparser

# Import our library
import ngesh


def new_tree(args):
    """
    Generates and returns a new tree.

    This function is a simple wrapper on the main `gen_tree()` and
    `add_characters()` functions. It is intented for command-line usage,
    but it can be used in code for quick prototyping.

    Parameters
    ----------

    args : Namespace
        A namespace with the parameters for tree generation.

    Returns
    -------
    tree: ete3 object
        The randomly generated tree.
    """

    # Generate the random tree
    tree = ngesh.gen_tree(
        args.birth,
        args.death,
        min_leaves=args.min_leaves,
        max_time=args.max_time,
        labels=args.labels,
        seed=args.seed,
    )

    # Add characters if requested
    if args.num_chars:
        tree = ngesh.add_characters(
            tree,
            args.num_chars,
            args.k_mut,
            args.th_mut,
            args.k_hgt,
            args.th_hgt,
            args.e_mut,
        )

    return tree


def parse_arguments():
    """
    Parses arguments and returns a namespace.

    This is a normal function for command-line argument passing relying on
    `argparse`. It allows to specify defaults in code, which are optionally
    overridden by a configuration file first (whose path is provided as
    an argument itself), and by command-line parameters second.

    Returns
    -------
    args : namespace
        A namespace with all the parameters.
    """

    # Specify the defaults for the options: these are overridden, in order,
    # by the configuration file and by the command line arguments
    options = {
        "labels": "human",
        "max_time": None,
        "min_leaves": 10,
        "seed": None,
        "birth": 1.0,
        "output": "newick",
        "num_chars": 0,
        "k_mut": 5.0,
        "th_mut": 1.0,
        "e_mut": 1.05,
    }

    # Parse any configuration speification first. Note that `add_help` is
    # set to False so that a call with `-h` is not parsed here; also
    # note that the `formatter_class` uses the raw description, so it
    # does not mess up with its own readability assumptions.
    conf_parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
    )
    conf_parser.add_argument(
        "-c", "--conf_file", help="Specify config file", metavar="FILE"
    )

    # Parse with `conf_parser`, looking for a configuration file; default
    # options will be those in `defaults` overridden by any value
    # specified in the configuration file
    args, remaining_argv = conf_parser.parse_known_args()
    if args.conf_file:
        config = configparser.SafeConfigParser()
        config.read([args.conf_file])
        options.update(dict(config.items("Config")))

    # Parse any remaining argument, this time without supressing the `--help`;
    # the options of `config_parser` are inherited through `parents`,
    # and defaults are set from `options` at the end of the new argument
    # creation.
    parser = argparse.ArgumentParser(parents=[conf_parser])

    parser.add_argument(
        "-b", "--birth", type=float, help="Set birth rate (l, default 1.0)"
    )
    parser.add_argument(
        "-d",
        "--death",
        type=float,
        help="Set death rate (mu, default 0.5 * `birth-rate`)",
    )
    parser.add_argument(
        "-t",
        "--max_time",
        type=float,
        help="Set maximum time stopping criterion",
    )
    parser.add_argument(
        "-l",
        "--min_leaves",
        type=float,
        help="Set minimum leaf number stopping criterion (defaults 10)",
    )
    parser.add_argument(
        "-x",
        "--labels",
        type=str,
        help='Set text generation model (defaults "human")',
    )
    parser.add_argument("-r", "--seed", type=str, help="Set RNG seed string")
    parser.add_argument(
        "-n",
        "--num_chars",
        type=int,
        help="Set random character number (default 0)",
    )
    parser.add_argument(
        "--k_mut",
        type=float,
        help="Set character mutation gamma `k` parameter (default 5.0)",
    )
    parser.add_argument(
        "--th_mut",
        type=float,
        help="Set character mutation gamma `theta` parameter (default 1.0)",
    )
    parser.add_argument(
        "--e_mut",
        type=float,
        help="Set character mutation `e` parameter (default 1.05)",
    )
    parser.add_argument(
        "--k_hgt",
        type=float,
        help="Set character HGT gamma `k` parameter (default 1.5 * `k_mut`)",
    )
    parser.add_argument(
        "--th_hgt",
        type=float,
        help="Set character HGT gamma `theta` parameter (default `th_mut`)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        choices=["newick", "ascii", "nexus", "wl"],
        help='Set output type (default "newick")',
    )

    parser.set_defaults(**options)

    args = parser.parse_args(remaining_argv)

    # Set death to half birth, if not provided
    if not args.death:
        args.death = args.birth / 2.0

    return args


def main():
    """
    Main function for tree generation from the command line.
    """

    # Parse command-line arguments
    args = parse_arguments()

    # Generate the tree
    tree = new_tree(args)

    # Output the tree according to the requested format
    if args.output == "newick":
        print(tree.write())
    elif args.output == "ascii":
        print(tree)
    elif args.output == "nexus":
        print(ngesh.tree2nexus(tree))
    elif args.output == "wl":
        print(ngesh.tree2wordlist(tree))


if __name__ == "__main__":
    main()
