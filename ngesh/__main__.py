# __main__.py

import random_tree

def main():
    print("ngesh main")
    bdtree = random_tree.gen_tree_safe(1.0, 0.5, max_time=3, human_labels=True)
    print(bdtree)

if __name__ == "__main__":
    main()
