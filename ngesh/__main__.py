# __main__.py

import ngesh

def main():
    print("ngesh main")
    bdtree = ngesh.gen_tree_safe(1.0, 0.5, max_time=3, human_labels=True)
    print(bdtree)

if __name__ == "__main__":
    main()
