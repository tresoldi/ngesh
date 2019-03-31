# __main__.py

import random

import ngesh

def main():
    # We make sure the death rate is less than the birth one.
    l = random.random()
    mu = random.random() * l * 0.9

    tree = ngesh.gen_tree(l, mu, max_time=4.0, labels="bio")
    print(tree.write())


if __name__ == "__main__":
    main()
