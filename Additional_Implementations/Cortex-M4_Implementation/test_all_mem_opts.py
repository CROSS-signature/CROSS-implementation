#!/usr/bin/env python3

import subprocess

OPTS_SIGN = ["MEM_OPT_SIGN_INC_CMT_1", "MEM_OPT_SIGN_RECOMP_E_V_U", "MEM_OPT_SIGN_RECOMP_Y", "MEM_OPT_SIGN_INC_MTREE", "MEM_OPT_SIGN_INC_STREE", "MEM_OPT_SIGN_BS_SEEDS_FROM_TREE", "MEM_OPT_HASH_TO_MTREE", "MEM_OPT_IN_PLACE_SAMP"]
OPTS_VERIFY = ["MEM_OPT_VERIFY_INC_CMT_1", "MEM_OPT_VERIFY_INC_Y", "MEM_OPT_VERIFY_INC_STREE", "MEM_OPT_VERIFY_INC_MTREE", "MEM_OPT_VERIFY_BS_SEEDS_FROM_TREE", "MEM_OPT_HASH_TO_MTREE", "MEM_OPT_IN_PLACE_SAMP"]

KAT = False

if KAT is True:
    kat = "-DGEN_KAT=1"
else:
    kat = ""


for i in range(2**len(OPTS_SIGN)):
    result = []
    for j in range(len(OPTS_SIGN)):
        if ((i >> j) & 1):
            result.append(OPTS_SIGN[j])
    cot = "-DCUSTOM_OPT_TARGET= -DDUMMY=1"
    os = "-DOPT_STRING="
    for o in result:
        cot = cot + r' -D' + o + "=1"
        os = os + o + "__"

    # Make all w. current target
    print("Building" + os)
    out = subprocess.run(['rm -rf ./*'], capture_output=True, cwd='build', shell=True)
    out = subprocess.run(['ls'], capture_output=True, cwd='build', text=True)
    out = subprocess.run(['cmake', '..', cot, os, kat], capture_output=True, cwd='build')
    # print(out)
    out = subprocess.run(['make', '-j16'], capture_output=True, cwd='build', check=True)
    # print(out)
    out = subprocess.run(['ls'], capture_output=True, cwd='build/bin', text=True)
    bins = out.stdout.splitlines()
    for binary in bins:
        print("Testing: " + binary + ":", end='')
        out = subprocess.run(['./' + binary], capture_output=True, cwd='build/bin', check=True)
        print(" Pass!")
    if KAT is True:
        # Copy KATs and verify hashes...
        out = subprocess.run(['mv *.req ../../../../KAT/'], capture_output=True, cwd='build/bin', shell=True)
        out = subprocess.run(['mv *.rsp ../../../../KAT/'], capture_output=True, cwd='build/bin', shell=True)
        # print(out)
        print("Testing KATS: ", end='')
        out = subprocess.run(['sha512sum', '-c', 'sha_512_sum_KATs'], capture_output=True, cwd='../../KAT', check=True)
        print("Pass")
        # print(out)

for i in range(2**len(OPTS_VERIFY)):
    result = []
    for j in range(len(OPTS_VERIFY)):
        if ((i >> j) & 1):
            result.append(OPTS_VERIFY[j])
    cot = "-DCUSTOM_OPT_TARGET= -DDUMMY=1"
    os = "-DOPT_STRING="
    for o in result:
        cot = cot + r' -D' + o + "=1"
        os = os + o + "__"

    # Make all w. current target
    print("Building" + os)
    out = subprocess.run(['rm -rf ./*'], capture_output=True, cwd='build', shell=True)
    # print(out)
    out = subprocess.run(['ls'], capture_output=True, cwd='build', text=True)
    # print(out)
    out = subprocess.run(['cmake', '..', cot, os, kat], capture_output=True, cwd='build')
    # print(out)
    out = subprocess.run(['make', '-j16'], capture_output=True, cwd='build', check=True)
    out = subprocess.run(['ls'], capture_output=True, cwd='build/bin', text=True)
    bins = out.stdout.splitlines()
    for binary in bins:
        print("Testing: " + binary + ":", end='')
        out = subprocess.run(['./' + binary], capture_output=True, cwd='build/bin', check=True)
        print(" Pass!")
    if KAT is True:
        # Copy KATs and verify hashes...
        out = subprocess.run(['mv *.req ../../../../KAT/'], capture_output=True, cwd='build/bin', shell=True)
        out = subprocess.run(['mv *.rsp ../../../../KAT/'], capture_output=True, cwd='build/bin', shell=True)
        # print(out)
        print("Testing KATS: ", end='')
        out = subprocess.run(['sha512sum', '-c', 'sha_512_sum_KATs'], capture_output=True, cwd='../../KAT', check=True)
        print("Pass")
        # print(out)
