#!/usr/bin/env python3
# encoding: utf-8
#
# __main__.py

# Import Python standard libraries
import argparse
import configparser
import random
import sys

# Import our library 
import ngesh

def generate_tree(args):
    # check arguments from the command line: at least one stopping
    # criterion, maybe death less than birth, etc.

    tree = ngesh.gen_tree(args.birth, args.death,
        min_leaves = args.min_leaves,
        max_time = args.max_time,
        labels = args.labels,
        seed = args.seed
        )

    return tree.write()


def parse_arguments():
    # Specify the defaults for the options: these are overridden, in order,
    # by the configuration file and by the command line arguments        
    options = {
        "labels" : "human",
        "max_time" : None,
        "num_leaves" : None,
        "seed" : None,
    }
        
    # Parse any configuration speification first. Note that `add_help` is
    # set to False so that a call with `-h` is not parsed here; also
    # note that the `formatter_class` uses the raw description, so it
    # does not mess up with its own readability assumptions.
    conf_parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class = argparse.RawDescriptionHelpFormatter,
        add_help = False
        )
    conf_parser.add_argument("-c", "--conf_file",
        help="Specify config file",
        metavar="FILE")

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
    parser = argparse.ArgumentParser(
        parents=[conf_parser]
        )
        
    parser.add_argument("-b", "--birth",
        type=float,
        help="Specify the birth rate (l)")
    parser.add_argument("-d", "--death",
        type=float,
        help="Specify the death rate (mu, default to half the birth-rate)")    
    parser.add_argument("-t", "--max_time",
        type=float,
        help="Specify the maximum time stopping criterion")  
    parser.add_argument("-l", "--min_leaves",
        type=float,
        help="Specify the minimum leaf number stopping criterion")
    parser.add_argument("-x", "--labels",
        type=str,
        help='Specify the model for text generation (defaults to "human")')
    parser.add_argument("-r", "--seed",
        type=str,
        help='Specify a string as the RNG seed')  
        
    parser.set_defaults(**options)
    
    args = parser.parse_args(remaining_argv)
    
    # Set death to half birth, if not provided
    if not args.death:
        args.death = args.birth / 2.0

    return args
    

def main():
    # Parse command-line arguments
    args = parse_arguments()

    # Generate the tree
    tree = generate_tree(args)


    print(tree)

if __name__ == "__main__":
    main()
    

