#!/usr/bin/python3
# This script computes all the derived parameters for CROSS which cannot be
# automated via C preprocessor macros
# Such parameters include the number of nodes for each level of the seed/Merkle
# trees and the amount of randomness required to be extracted from the CSPRNGs


# arbitrary precision floating point library
from mpmath import mp
from math import comb,log2,ceil,floor
import csv
from sys import argv,exit

# number of bits to represent a given number
def bits_to_represent(p):
    return max(1, p.bit_length())

# number of bits to represent an element of Z_p
def bits_for_element_mod(p):
    # it's the number of bits to represent p-1
    return bits_to_represent(p-1)

def mpmath_log2(a):
    return mp.log(a)/mp.log(mp.mpf(2))

def clog2(a):
    return max(int(ceil(log2(a))), 1)

def l_child(a):
    return 2*a + 1

def r_child(a):
    return 2*a + 2


# find number of bits to draw from the CSPRNG to obtain at least n
# values mod p via rejection sampling, with Pr at least 1-2**-lambda
def find_num_csprng_bits(n,p,lam):
    bitlen_elem = bits_for_element_mod(p)
    pr_succ = mp.mpf(p/2**bitlen_elem)
    pr_insucc = 1-pr_succ
    # we at least need n attempts in rejection sampling
    n_attempts = n
    # iterate until the cumulative probability of getting enough successes is
    # higher than 1-2**-lambda
    cumul_pr = mp.mpf(0)
    while(cumul_pr < 1-mp.mpf(2)**-lam):
        n_attempts += 1
        cumul_pr = mp.mpf(0)
        # add the probabilities of all the cases where at least n successes
        # take place, i.e., at least n good values are sampled
        for n_ok in range(n, n_attempts+1):
            cumul_pr += mp.mpf(comb(n_attempts,n_ok))* mp.power(pr_succ,n_ok) * mp.power(pr_insucc,(n_attempts-n_ok))
    # always shift by the amount of bits required for a sample
    total_bits_required = n_attempts*bitlen_elem
    return total_bits_required


# find number of bits to draw from the CSPRNG to obtain at least t-w distinct
# positions in a t bits long bitstring via rejection sampling,
# with Pr at least 1-2**-lambda
def find_num_csprng_bits_cw_string(t,p,lam):
    n_iter_FY = t-1
    successes_to_achieve = n_iter_FY
    cumul_pr = mp.mpf(0)
    while(cumul_pr < 1-mp.mpf(2)**-lam):
        n_iter_FY += 1
        # discrete distribution w/range representing number of succesfully placed
        # positions
        discrete_prob = [mp.mpf(0)]*(successes_to_achieve+1)
        discrete_prob[0] = mp.mpf(1)
        for step in range(n_iter_FY):
            next_discrete_prob = [mp.mpf(0)]*(successes_to_achieve+1)
            for n_succ in range(successes_to_achieve):
                remaining = successes_to_achieve-n_succ
                bits_for_rem = bits_to_represent(remaining)
                # probability that a binary string is an acceptable index in 0...remaining-1
                pr_ok = mp.mpf(remaining)/mp.mpf(2**bits_for_rem)
                pr_ko = 1-pr_ok
                next_discrete_prob[n_succ]   += discrete_prob[n_succ] * mp.mpf(pr_ko)
                next_discrete_prob[n_succ+1] += discrete_prob[n_succ] * mp.mpf(pr_ok)
            # mass of t successes is preserved and added
            next_discrete_prob[successes_to_achieve] += discrete_prob[successes_to_achieve]
            discrete_prob = next_discrete_prob
        cumul_pr = discrete_prob[successes_to_achieve]
        # print(f"natt {n_iter_FY}, prob {mp.log(mp.mpf(1)-cumul_pr,2)}")
        # Overestimate by using the number of bits to represent the highest index in 
        # each iteration.
        # TODO: Rewrite for better estimation
    return n_iter_FY*bits_to_represent(t-1)


# Compute the offsets for the truncated trees required to move between two levels
def tree_offsets_and_nodes(T):

    # Full trees on the left half, so we can already count (i.e. subtract) these values as well as the root node
    missing_nodes_per_level = [2**(i-1) for i in range(1, clog2(T)+1)]
    missing_nodes_per_level.insert(0,0)

    remaining_leaves = T - 2**(clog2(T)-1)
    level = 1

    # Starting from the first level, we construct the tree in a way that the left
    # subtree is always a full binary tree.
    while(remaining_leaves > 0):
        depth = 0
        stree_found = False
        while not stree_found:
            if (remaining_leaves <= 2**depth):
                for i in range(depth, 0, -1):
                    missing_nodes_per_level[level+i] -= 2**(i-1)
                remaining_leaves -= (2**clog2(remaining_leaves)) // 2

                # Subtract root and increase level for next iteration
                missing_nodes_per_level[level] -= 1
                level += 1
                stree_found = True
            else:
                depth += 1
            
    # The offsets are the missing nodes per level subtracted by the missing nodes of all previous levels, as this
    # is already included 
    offsets = [missing_nodes_per_level[i] for i in range(len(missing_nodes_per_level))]
    for i in range(clog2(T), -1, -1):
        for j in range(i):
            offsets[i] -= offsets[j]

    nodes_per_level = [2**i - missing_nodes_per_level[i] for i in range(clog2(T)+1)]
    return offsets, nodes_per_level

# Compute the number of subtrees and corresponding start indices of the leaf nodes within
# the full tree.
def tree_leaves(T, offsets):
    leaves = [0]*T
    leaves_per_level = [0]*(clog2(T)+1)
    start_index_per_level = [0]*(clog2(T)+1)
    ctr = 0

    remaining_leaves = T
    depth = 0
    level = 0
    root_node = 0
    left_child = l_child(root_node) - offsets[level+depth]
    
    while (remaining_leaves > 0):
        depth = 1
        subtree_found = False
        while not subtree_found:
            if (remaining_leaves <= 2**depth):
                for i in range(2**clog2(remaining_leaves)//2):
                    leaves[ctr] = root_node if remaining_leaves==1 else left_child+i
                    if (remaining_leaves==1):
                        leaves_per_level[level] += 1
                        start_index_per_level[level] = root_node if start_index_per_level[level] == 0 else start_index_per_level[level]
                    else:
                        leaves_per_level[level+depth] += 1
                        start_index_per_level[level+depth] = left_child if start_index_per_level[level+depth] == 0 else start_index_per_level[level+depth]
                    ctr += 1
                root_node = r_child(root_node) - offsets[level]
                left_child = l_child(root_node) - offsets[level]
                level += 1
                remaining_leaves -= 2**clog2(remaining_leaves)//2
                subtree_found = True
            else:
                left_child = l_child(left_child) - offsets[level+depth]
                depth += 1

    # Now create array with start idx and number of leaves by removing zeros
    cons_leaves = [i for i in leaves_per_level if i != 0]
    start_index_per_level = [i for i in start_index_per_level if i != 0]

    return leaves_per_level, len(cons_leaves), start_index_per_level[::-1], cons_leaves[::-1]

def hamming_weight_of_value(x):
    hw = 0
    while x > 0:
        hw += x & 1
        x = x >> 1
    return hw

def worst_case_tree_nodes(num_rounds_t,seeds_to_hide):
    num_set_bits_w = num_rounds_t-seeds_to_hide
    weight_t_1 = hamming_weight_of_value(num_rounds_t) - 1
    return floor(seeds_to_hide*log2(num_rounds_t/seeds_to_hide) + weight_t_1)

def print_tree_defines(npl,off,lpl,start_indices,cons_leaves,nodes_to_store):
    c_array_npl = repr(npl).replace("[","{").replace("]","}")
    c_array_off = repr(off).replace("[","{").replace("]","}")
    c_array_lpl = repr(lpl).replace("[","{").replace("]","}")
    c_array_start_indices = repr(start_indices).replace("[","{").replace("]","}")
    c_array_cons_leaves = repr(cons_leaves).replace("[","{").replace("]","}")
    print(f"#define TREE_OFFSETS {c_array_off}")
    print(f"#define TREE_NODES_PER_LEVEL {c_array_npl}")
    print(f"#define TREE_LEAVES_PER_LEVEL {c_array_lpl}")
    print(f"#define TREE_SUBROOTS {subroots}")
    print(f"#define TREE_LEAVES_START_INDICES {c_array_start_indices}")
    print(f"#define TREE_CONSECUTIVE_LEAVES {c_array_cons_leaves}")
    print(f'#define TREE_NODES_TO_STORE {nodes_to_store}')


if __name__=="__main__":
    # arbitrary precision floats mantissa precision (in decimal digits)
    mp.dps = 200
    # check that a security level is provided
    if (len(argv)< 2):
        print("Script computing the auxiliary parameters for CROSS")
        print(f"Usage:{argv[0]} <CSV with code parameters>")
        print(f"Expected CSV header: lambda,opt,p,z,n,m,t")
        print("opt is one out of SPEED,BALANCED,SIG_SIZE")
        print(f"m=0 marks a RSDP instance")
        exit()

    try:
        parameter_file = open(argv[1],"r")
    except:
        print(f"parameter file {argv[1]} not found")

    print("\n/* Derived parameters computed via compute_derived_parameters.py */")
    csv_field_names=["lambda","opt","p","z","n","k","m","t","w"]
    csvreader = csv.DictReader(parameter_file,fieldnames=csv_field_names)
    for line in csvreader:
        if (line["lambda"]=="lambda"):
            continue
        seclevel = int(line["lambda"])
        opt = line["opt"]
        p = int(line["p"])
        z = int(line["z"])
        n = int(line["n"])
        k = int(line["k"])
        m = int(line["m"])
        t = int(line["t"])
        w = int(line["w"])

        category_strings = {128: "CATEGORY_1",192: "CATEGORY_3",256: "CATEGORY_5"}

        if(seclevel == 128) and (m == 0) and (opt == "SPEED"):
            print(f"#if ( defined({category_strings[seclevel]}) ",end='')
        else :
            print(f"\n#elif ( defined({category_strings[seclevel]}) ",end='')

        if( m == 0 ):
            print(f"&& defined(RSDP) ",end='')
        else :
            print(f"&&  defined(RSDPG) ",end='')

        print(f" && defined({opt}) )")

        off, npl = tree_offsets_and_nodes(t)
        lpl, subroots, start_idx, cons_leaves = tree_leaves(t, off)
        if opt == "SPEED":
            nodes_to_store = w
        else:
            nodes_to_store = worst_case_tree_nodes(t,t-w)
        print_tree_defines(npl,off,lpl,start_idx,cons_leaves,nodes_to_store)

        # u' is n elems mod p
        bits_uprime = find_num_csprng_bits(n,p,seclevel)
        print(f"#define BITS_N_FP_CT_RNG {bits_uprime}")
        bits_beta = find_num_csprng_bits(t,p-1,seclevel)
        print(f"#define BITS_BETA_FPSTAR_CT_RNG {bits_beta}")

        bits_V_tr = find_num_csprng_bits((n-k)*k,p,seclevel)
        print(f"#define BITS_V_CT_RNG {bits_V_tr}")
        if (m==0):
            # eta is n elems mod z
            bits_etaprime = find_num_csprng_bits(n,z,seclevel)
            tot=bits_uprime+bits_etaprime
            print(f"#define BITS_N_FZ_CT_RNG {bits_etaprime}")
        else:
            bits_W = find_num_csprng_bits((n-m)*m,z,seclevel)
            print(f"#define BITS_W_CT_RNG {bits_W}")
            # zeta is m elems mod z
            bits_zetaprime = find_num_csprng_bits(m,z,seclevel)
            tot=bits_uprime+bits_zetaprime
            print(f"#define BITS_M_FZ_CT_RNG {bits_zetaprime}")
        num_bits_cwstring = find_num_csprng_bits_cw_string(t,p,seclevel)
        print(f"#define BITS_CWSTR_RNG {num_bits_cwstring}")
    print("\n#endif")

