#!/usr/bin/env python3
# encoding: utf-8
#
# __main__.py

# Import Python standard libraries
import random

# Import our library 
import ngesh

def main():
    # We make sure the death rate is less than the birth one.
    l = random.random()
    mu = random.random() * l * 0.9

    tree = ngesh.gen_tree(l, mu, max_time=4.0, labels="human")
    print(tree.write())

    import sys
    print(sys.argv)

if __name__ == "__main__":
    main()
    

