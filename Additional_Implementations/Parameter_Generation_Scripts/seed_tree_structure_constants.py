#!/usr/bin/python3
from sys import argv, exit
from math import log2, ceil

# /* level offsets contains the indexes of the first seed of a given level,
#  * nodes per level contains the number of nodes for each level. All these can
#  * be precomputed and declared as constants: a runtime computation is provided for
#  * reference purposes. The function returns the total number of tree nodes for
#  * VLA allocation */

def compute_nodes_per_level_and_offsets(T):
    nodes_per_level = []
    nodes_curr_level = T
    total_nodes = 0
    for level in range(ceil(log2(T)),-1,-1):
        nodes_per_level = [nodes_curr_level] + nodes_per_level
        total_nodes += nodes_curr_level
        nodes_curr_level = ceil(nodes_curr_level/2)

    cumulative_offset = 0
    total_saved_nodes = 0
    level_offsets = []
    # Level 0, the root has no previous levels, no missing nodes before it
    missing_nodes_before_level = []
    cumul_saved_nodes = []
    for level in range(ceil(log2(T))+1):
        level_offsets.append(cumulative_offset)
        cumulative_offset += nodes_per_level[level]
        if (level == 0):
            missing_nodes_before_level.append(0)
        else:
            missing_nodes_before_level.append(total_saved_nodes)
        total_saved_nodes += 2**level-nodes_per_level[level]

        cumul_saved_nodes.append(total_saved_nodes)

    return total_nodes,total_nodes-nodes_per_level[-1],nodes_per_level,level_offsets,missing_nodes_before_level

def worst_case_seed_tree_cost(t,w):
    num_not_to_reveal = t-w
    return ceil(num_not_to_reveal* log2(t/(num_not_to_reveal)))


def print_defines(nodes_to_store,npl,mnbl,num_nodes):
    print("\n/* determined via seed_tree_structure_constants.py script */")
    print(f"#define TREE_NODES_TO_STORE {nodes_to_store}")
    print(f"#define NUM_NODES_SEED_TREE {num_nodes}")
    c_array_npl = repr(npl).replace("[","{").replace("]","}")
    print(f"#define NODES_PER_LEVEL_ARRAY {c_array_npl}")
    c_array_mnbl = repr(mnbl).replace("[","{").replace("]","}")
    print(f"#define MISSING_NODES_BEFORE_LEVEL_ARRAY {c_array_mnbl}")

if __name__ == "__main__":
    if len(argv) != 3:
        print("Script to compute seed tree constants")
        print("Computes nodes per level and seed offsets")
        print(f"Usage: {argv[0]} <number_of_leaves_T> <num_set_bits_W>")
        exit()
    num_leaves_T = int(argv[1])
    num_leaves_W = int(argv[2])
    print(f"computing for {num_leaves_T} leaves")
    num_nodes,NUM_INNER_NODES_OF_SEED_TREE,npl,lo,mnbl = compute_nodes_per_level_and_offsets(num_leaves_T)
    nodes_to_store = worst_case_seed_tree_cost(num_leaves_T,num_leaves_W)
    print(f"ceil(log2(t)) {ceil(log2(num_leaves_T))}")
    print(f"nodes of seed tree {num_nodes}")
    print(f"inner nodes of seed tree {NUM_INNER_NODES_OF_SEED_TREE}")
    print(f"nodes_per_level {npl}")
    print(f"start of level {lo}")
    print(f"MNBL {mnbl}")
    print_defines(nodes_to_store,npl,mnbl,num_nodes)
