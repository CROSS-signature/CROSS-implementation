def rts(number):
    """
    rts = round (float) to string
    """
    s = str(round(number,2))
    return (s + "0" * (3 + s.find(".") -len(s)))

##########################################
# functions for analyzing structure of E #
##########################################
def get_Eset(p,z):
    """
    generate mathbb{E} for (p,z)-pair
    """
    assert is_prime(p),"field size p is not prime!"
    assert (p-1)%z == 0, "z does not divide p-1!"

    F = GF(p)
    alpha = F.primitive_element()
    g = alpha**((p-1)/z)
    return [g^l for l in range(z)]

def shift_E(E, shift):
    return [e-shift for e in E if (e-shift) != 0]

def get_Dset(Eset):
    """
    generates subset D of difference set {a-b|a,b in E}
    """
    F = Eset[0].parent()
    z = len(Eset)
    F_list = F.list()
    p = F.cardinality()
    counts = [0 for _ in range(p)]
    for i in range(p):
        a = F_list[i]
        for b in Eset:
            x = a + b
            if x in Eset:
                counts[i] += 1

    counts[0] = -1
    for i in range(z):
        counts[Eset[i]] = -1
    
    Dset = []
    for i in range(p):
        if counts[i] == max(counts):
            Dset.append(F_list[i])
    return Dset

def get_alphaE(Eset):
    """
    determine additive structure of E
    """
    z = len(Eset)
    counts = [0 for _ in range(z)]
    for a in Eset:
        for b in Eset:
            x = a+b 
            if x in Eset:
                counts[Eset.index(x)] += 1

    irregular = 0
    for l in range(1,z):
        if counts[l] != counts[0]:
            irregular = 1
    
    assert irregular == 0,"additive structure of E not regular!"

    return counts[0]

def get_alphaD(Eset, Dset):
    """
    determine additive structure of E w.r.t. D
    """
    z = len(Eset)
    counts = [0 for _ in range(z)]
    for a in Eset:
        for b in Dset:
            x = a+b 
            if x in Eset:
                counts[Eset.index(x)] += 1

    irregular = 0
    for l in range(1,z):
        if counts[l] != counts[0]:
            irregular = 1
    
    assert irregular == 0,"additive structure of E wrt D not regular!"

    return counts[0]

#########################################################
# functions for finding parameter sets and making plots #
#########################################################


def find_params(p, z, nlow, nhigh, nist, solver, Mmax=-1):
    """
    find parameters which satisfy NIST I, III or V wrt to specified solver

    inputs:
    p:          field size
    z:          size of restricted set
    nlow:       minimum considered code length
    nhigh:      maximum considered code length
    nist:       NIST security level in {1,3,5}
    solver:     algorithm which is considered as solver
                'stern', 'shifted_stern', 'bjmm' or 'shifted_bjmm'
    Mmax:   maximum memory restriction, deactivated if < 0

    output:
    n, k that achieve chosen security level and parameters of solver
    """
    assert nist in [1,3,5],"nist has to be 1, 3 or 5"
    if   nist == 1: seclvl = 143
    elif nist == 3: seclvl = 207
    elif nist == 5: seclvl = 272

    # max rate such that solution unique
    Rmax = 1 - log(z,2)/log(p,2)

    for n in range(nlow, nhigh+1):
        print("checking n=",n," with k in {",max(0, floor(n*(Rmax-0.1))),",...,",min(ceil(n*(Rmax+0.1)), n),"}")
        for k in range(max(0, floor(n*(Rmax-0.1))), min(ceil(n*(Rmax+0.1)), n)+1):

            if solver == 'stern':           C, M, P = opt_stern(p, z, n, k, Mmax, False)
            if solver == 'shifted_stern':   C, M, P = opt_shifted_stern(p, z-1, n, k, Mmax, False)
            if solver == 'bjmm':            C, M, P = opt_bjmm(p, z, n, k, Mmax, False)
            if solver == 'shifted_bjmm':    C, M, P = opt_shifted_bjmm(p, z-1, n, k, Mmax, False)
            print(n,k,C)
            if C > seclvl:
                return C, N(M), P

def make_plot(p, z, R, nlow, nhigh, solver, Mmax=-1):
    """
    reproduce Figure 7 of documentation
    p:          field size
    z:          size of restricted set
    R:          code rate
    nlow:       minimum considered code length
    nhigh:      maximum considered code length
    solver:     algorithm which is considered as solver
                'stern', 'shifted_stern', 'bjmm' or 'shifted_bjmm'
    Mmax:   maximum memory restriction, deactivated if < 0
    """
    C = -1
    for n in range(nlow, nhigh+1):
        k = round(n*R)
        if solver == 'stern':           C = opt_stern(p, z, n, k, Mmax, False)[0]
        if solver == 'shifted_stern':   C = opt_shifted_stern(p, z-1, n, k, Mmax, False)[0]
        if solver == 'bjmm':            C = opt_bjmm(p, z, n, k, Mmax, False)[0]
        if solver == 'shifted_bjmm':    C = opt_shifted_bjmm(p, z-1, n, k, Mmax, False)[0]
        print("(",n,",",rts(C),") %k=",k,"nsol=",rts(z**n/p**(n-k)))
    return 1