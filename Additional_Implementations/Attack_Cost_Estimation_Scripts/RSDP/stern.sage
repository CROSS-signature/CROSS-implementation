##########################################################
#     FINITE REGIME ESTIMATE of COMPLEXITY of Stern      #
##########################################################

def cost_stern(p, z, n, k, l,  Mmax, verb):
    """
    inputs:
    p: field size
    z: size of restricted set
    n: code length
    k: code dimension
    l: redundancy of small instance
    verb:   verbose

    output: estimate of time and space complexity of Stern in bits
    """ 
    if verb: 
        print("n=",n,"k=",k)
        print("l=",l)
    # number of solutions
    nsol = N(z**n * p**(k-n))+1
    if verb: print("#sol=", N(nsol))

    # assumes full error weight
    v      = floor((k+l)/2)
    vprime =  ceil((k+l)/2)

    # sizes of base lists
    L      = z**v
    Lprime = z**vprime

    # cost of base lists
    C       = L      * (v*log(z,2)      + l*log(p,2))
    Cprime  = Lprime * (vprime*log(z,2) + l*log(p,2))

    # memory cost in bits
    Mstern = N(log(L * v * log(z,2),2))
    if verb:  print("M_stern=",rts(Mstern), rts(Mstern/n))
    if log(Mstern,2) > Mmax and Mmax > 0: return 10**6, Mstern

    # cost of collision search:
    Ccoll = L * Lprime * p**(-l) *(k+l)*log(p,2)

    if verb: print("C=",N(log(C,2)/n),"C'=",N(log(Cprime,2)/n),"C_coll=",N(log(Ccoll,2)/n))


    # overall cost estimate
    Cstern = N(log((C + Cprime + Ccoll)*Mstern / nsol, 2))
    if verb:  print("C_stern=",Cstern, N(Cstern/n))

    return Cstern, Mstern

def opt_stern(p, z, n, k, Mmax = -1, verb = False):
    """
    optimize parameters of stern

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
    P_opt =  {'ell': -1}
    for l in range(n-k+1):
            C, M = cost_stern(p, z, n, k, l, Mmax, False)
            if C < C_opt:
                C_opt = C
                M_opt = M
                P_opt =  {'ell': l}


    if verb: cost_stern(p, z, n, k, P_opt['ell'],  Mmax, True)

    return C_opt, M_opt, P_opt


