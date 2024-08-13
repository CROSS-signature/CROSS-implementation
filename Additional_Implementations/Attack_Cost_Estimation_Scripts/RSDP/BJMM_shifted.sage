##########################################################
# FINITE REGIME ESTIMATE of COMPLEXITY of SHIFTED BJMM   #
#               TWO REPRESENTATION LAYERS                #
##########################################################


def cost_shifted_bjmm(p, z, n, k, l, v0, delta1, delta2, Mmax = -1, verb = False):
    """
    computes conservative estimate on complexity of BJMM algorithm for shifted instances

    inputs:
    p:              field size
    z:              size of restricted set
    n:              code length
    k:              code dimension
    l:              redundancy of small instance
    delta1, delta2: overlap with D on level i
    Mmax:           maximum memory restriction, deactivated if < 0
    verb:           verbose

    output:
    lower bound on time complexity C and memory complexity M of BJMM algorithm
    """ 
    if verb:
        print("l =",l)
        print("v0=",v0)
        print("delta1 =",delta1,"delta2 =",delta2)

    assert z == 6,"only shifted z=7 possible!"
    # additive structure
    global alphaE, alphaD, zd
    if verb: print("alphaE=",alphaE,"zd=",zd,"alphaD=",alphaD)

    # logarithmic values of sizes
    Z  = log(z, 2)
    ZD = log(zd,2)
    P  = log(p, 2)

    # set weights for all levels, make sure they satisfy divisibility as needed
    if v0 > k+l: return 10**6+1, 10**6

    if v0%2: return 10**6+1, 10**6
    v1 = v0/2
    if v1%2: return 10**6+2, 10**6
    v2 = v1/2
    if v2%2: return 10**6+3, 10**6
    v3 = v2/2 

    d1 = delta1
    if d1%2: return 10**6+4, 10**6
    d2 = d1/2 + delta2
    if d2%2: return 10**6+5, 10**6
    d3 = d2/2

    # number of solutions
    nsol = N(7**n * p**(k-n))+1
    if verb: print("#sol=", rts(nsol))

    # success probability
    Psucc = binomial(k+l,v0)*z**v0 * (z+1)**(-k-l)
    if verb: print("Psucc=", rts(log(Psucc,2)))

    # number of representation for entry in L0
    r0 = binomial(v0, v0/2)  * binomial(v0/2, delta1)**2 * alphaD**(2*delta1) 
    u0 = floor(log(r0,p))
    u0 = N(min(u0,l))

    # number of representation for entry in L1
    r1 = binomial(v1, v1/2) * binomial(v1/2, delta2)**2 * alphaD**(2*delta2) * binomial(delta1, delta1/2)
    l1 = floor(log(r1,p))
    l1 = N(min(l1,u0))
    l0 = u0 - l1
    if verb: print("l1=",rts(l1),"l0=",rts(l0))

    # Cost of level 3: size of a single base list
    L3 = multinomial(v3, d3, (k+l)/2-v3-d3) * z**v3 * zd**d3
    C3 = 2 * L3 * (l1*P + v3*Z + d3*ZD)
    if verb: print("L3=",rts(log(L3,2)/n),"C3=",rts(log(C3,2)/n))


    # Cost of level 2: size of a list on level 2
    L2 = N(L3^2 * p^(-l1))
    C2 = 2 * L2 * (l0*P + v2*Z + d2*ZD)
    if verb: print("L2=",rts(log(L2,2)/n),"C2=",rts(log(C2,2)/n))

    # Cost of level 1: cost of merging
    L1 = (multinomial(v1, d1, (k+l)-v1-d1) * z**v1 * zd**d1) * p**(-l0-l1)
    C1 = 2 * L2^2 * p**(-l0) * P
    if verb: print("L1=",rts(log(L1,2)/n),"C1=",rts(log(C1,2)/n))

    # Cost of level 0: final number of candidates for small instance
    C0 = L1^2 * p^(-(l-l0-l1)) * P
    if verb: print("C0=",rts(log(C0,2)/n))


    # logarithm of space complexity
    Mbjmm =  max(N(log(max( L3*(v3*Z+d3*ZD), L2*(v2*Z+d2*ZD), L1*(v1*Z+d1*ZD) ),2)),1)
    if verb:  print("space=", rts(Mbjmm), rts(Mbjmm/n))
    if Mbjmm > Mmax and Mmax > 0: return 10**6, Mbjmm

    # conservative estimate of total time complexity in bit
    Cbjmm = N(log((C3 + C2 + C1 + C0) * Mbjmm / nsol /Psucc, 2))
    if verb:  print("Cbjmm=", rts(Cbjmm), rts(Cbjmm/n))

    return Cbjmm, Mbjmm


def opt_shifted_bjmm(p, z, n, k, Mmax=-1, verb=False):
    """
    optimize parameters of BJMM solver on shifted instance
    inputs:
    p:      field size
    z:      size of restricted set
    n:      code length
    k:      code dimension
    Mmax:   maximum memory restriction, deactivated if < 0
    verb:   verbose

    ouput:
    optimized time complexity, memory cost, algoithm parameters
    """
    C_opt = 10**6
    M_opt = 10**6
    P_opt = {'ell': -1, 'v0': -1, 'delta1': -1, 'delta2': -1}
    for l in range(0, n-k+1):
        for v0 in range(0,k+l+1):
            if v0%8: continue
            for delta1 in range(0, min(v0/2, 17)+1): 
                d1 = delta1
                if d1%2: continue
                for delta2 in range(0, min(v0/4, 9)+1): 
                    d2 = d1/2+delta2
                    if d2%2:continue
                    C, M = cost_shifted_bjmm(p, z, n, k, l, v0, delta1, delta2, Mmax, False)
                    if C < C_opt:
                        M_opt       = M
                        C_opt       = C
                        P_opt = {'ell': l, 'v0': v0, 'delta1': delta1, 'delta2': delta2}

    if verb: cost_shifted_bjmm(p, z, n, k, P_opt['ell'], P_opt['v0'], P_opt['delta1'], P_opt['delta2'], Mmax, True)
    return C_opt, M_opt, P_opt

