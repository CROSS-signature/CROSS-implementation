##########################################################
#  FINITE REGIME ESTIMATE of COMPLEXITY of SHIFTED Stern #
##########################################################

def cost_shifted_stern(p, zshifted, n, k, l, v0a, v0b, Mmax=-1, verb=False):
    """
    inputs:
    p:      field size
    z:      size of restricted set
    n:      code length
    k:      code dimension
    l:      redundancy of small instance
    Mmax:   maximum memory restriction, deactivated if < 0
    verb:   verbose

    output: 
    estimate of time and space complexity of Stern in bits
    """ 
    if verb:
        print("l =",l)
        print("v0a=",v0a,"v0b=",v0b)

    # number of solutions
    nsol = N((zshifted+1)**n * p**(k-n))+1
    if verb: print("#sol=", N(nsol))

    # assumes full error weight
    sizea = floor((k+l)/2)
    sizeb =  ceil((k+l)/2)

    # success probability
    Psucc = binomial(sizea,v0a)*binomial(sizeb,v0b)*zshifted**(v0a+v0b) * (zshifted+1)**(-k-l)
    if verb: print("Psucc=", rts(log(Psucc,2)))

    # sizes of base lists
    La = binomial(sizea, v0a) * zshifted**v0a
    Lb = binomial(sizeb, v0b) * zshifted**v0b

    # cost of base lists
    Ca  = La * (v0a*log(zshifted,2) + l*log(p,2))
    Cb  = Lb * (v0b*log(zshifted,2) + l*log(p,2))

    # memory cost in bits
    Mstern = N(log(min(La*v0a,Lb*v0b) * log(zshifted,2), 2))
    if verb:  print("M_stern=",N(log(Mstern,2)),N(log(Mstern,2)/n))
    if Mstern > Mmax and Mmax > 0: return 10**6, Mstern

    # cost of collision search:
    Ccoll = La * Lb * p**(-l) *(k+l)*log(p,2)

    if verb: print("C=",N(log(Ca,2)/n),"C'=",N(log(Cb,2)/n),"C_coll=",N(log(Ccoll,2)/n))

    # overall cost estimate
    Cstern = N(log((Ca +  Cb + Ccoll)*Mstern / nsol /Psucc, 2))
    if verb:  print("C_stern=",Cstern, N(Cstern/n))

    return Cstern, Mstern

def opt_shifted_stern(p, zshifted, n, k, Mmax=-1, verb=False):
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
            for v0a in range(floor((k+l)/2)+1):
                for v0b in range(max(0, v0a-5), min(ceil((k+l)/2), v0a+5)+1):
                    C, M = cost_shifted_stern(p, zshifted, n, k, l, v0a, v0b, Mmax, False)
                    if C < C_opt:
                        C_opt = C
                        M_opt = M
                        P_opt =  {'ell': l, 'v0a':v0a, 'v0b':v0b}

    if verb: cost_shifted_stern(p, zshifted, n, k, P_opt['ell'], P_opt['v0a'], P_opt['v0b'], Mmax, True)

    return C_opt, M_opt, P_opt

