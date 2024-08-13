# Created for SageMath version 9.0                  
# Using Python 3.8.10


def get_num_subcodes(z, n, m, w, d): #get_num_subcodes(z,n,n-m,j,j-rho)
    """
    get the expected number of subcodes with given parameters

    inputs:
    z: field size (here z!)
    n: length of C
    m: dimension of code (here m)
    w: support size of subcode
    d: dimension of subcode

    output:
    expected number of subcodes with given parameters
    """
    return binomial(n,w) * (z^d-1)**(w-d) * q_binomial(m,d,z) / q_binomial(n,d,z)

def get_all_subcodes(z, n, m, nist):
    """
    find all subcodes of a certain code over F_z
    that could possibly speed up enumeration

    input:
    z: field size (here z!)
    n: code length
    m: code dimension (here m!)
    nist: NIST security level (only returns relevant subcodes)

    output:
    list of relevant subcodes
    """
    assert nist in [1,3,5],"nist has to be 1, 3 or 5"
    if   nist == 1: seclvl = 143
    elif nist == 3: seclvl = 207
    elif nist == 5: seclvl = 272

    Z = log(z, 2)

    subcodes = []
    # submatrices with width <0 m
    for j in range(1,m+1):
        subcodes += [{'j': j, 'rho': j,'P(j,rho)': 0.0}]
        for rho in range(1, j+1):
            if rho <= seclvl / Z:
                numcodes = get_num_subcodes(z,n,n-m,j,j-rho)
                P = N(log( min(numcodes, 1.0),2))
                if P > -seclvl:
                    subcodes += [{'j': j, 'rho': rho,'P(j,rho)': P}] 

    # submatrices with width > m
    for j in range(m+1, n+1):
        subcodes += [{'j': j, 'rho': m,'P(j,rho)': 0.0}]
        for rho in range(1, m+1):
            if rho <= seclvl / Z:
                numcodes = get_num_subcodes(z,n,m,n-j,m-rho)
                P = N(log( min(numcodes, 1.0),2))

                if P > -seclvl:
                    subcodes += [{'j': j, 'rho': rho,'P(j,rho)': P}] 
    return subcodes


def rts(number):
    """
    rts = round (float) to string
    """
    s = str(round(number,2))
    return (s + "0" * (3 + s.find(".") -len(s)))


def submatrix_stern(p, z, n, k, m, nist, Mmax):
    """
    inputs:
    p:      field size
    z:      size of restricted set
    n:      code length
    k:      code dimension
    m:      order of group G
    nist:   NIST security level [1,3,5]
    Mmax:   maximum memory restriction, deactivated if < 0

    output:
    time complexity, memory cost and parameters of optimized submatrix stern solver
    """
    assert nist in [1,3,5],"nist has to be 1, 3 or 5"
    if   nist == 1: seclvl = 143
    elif nist == 3: seclvl = 207
    elif nist == 5: seclvl = 272

    # find all relevant subcodes
    subcodes = get_all_subcodes(z, n, m, nist)
    numsubc = len(subcodes)

    C_opt = 10**9
    P_max = -1
    Cstern = 0

    Z = log(z,2)
    nsol = N(z**m * p**(k-n))+1

    # try all possible submatrices a 
    for ia in range(numsubc):
        sca   = subcodes[ia]
        j_a   = sca['j']
        rho_a = sca['rho']
        P_a   = sca['P(j,rho)']

        # try submatrices b 
        for ib in range(ia, numsubc):
            scb   = subcodes[ib]
            j_b   = scb['j']
            rho_b = scb['rho']
            P_b   = scb['P(j,rho)']

            if - P_b - P_a > C_opt: 
                continue # more costly already due to probability of occurrence

            l       = j_a   + j_b   - k
            ltilde  = max(rho_a + rho_b - m,0)
            if l < 0:
                continue

            # number of enumerated error vectors
            La = z**rho_a
            Lb = z**rho_b

            # costs associated with list a and b
            Ca = La * (rho_a*Z + ltilde*Z+ l*log(p,2))
            Cb = Lb * (rho_b*Z + ltilde*Z+ l*log(p,2))

            # cost of collision search
            Ccoll = La * Lb * p**(-l) *(k+l)*log(p,2) * z**(-ltilde)
            

            # memory cost
            Mstern = N(log(max(min(La*rho_a, Lb*rho_b) * Z, 2),2))
            if Mstern > Mmax and Mmax > 0:  continue

            # overall time complexity
            Cstern = N(log((Ca + Cb + Ccoll)*Mstern/nsol + 2**(- P_b - P_a) ,2)) #

            if Cstern < C_opt or (Cstern == C_opt and min(P_a, P_b) > P_max):
                P_max =  max(P_a, P_b)
                M_opt = Mstern
                C_opt = Cstern
                P_opt = {'ell': l, 'subcode_a': sca,  'subcode_b': scb}

    return C_opt, M_opt, P_opt

def find_params(p, z, nmin, nmax, nist, Mmax):
    """
    determine valid parameters

    inputs:
    p:      finite field size
    z:      size of restricted set
    nmin:   minimum considered codelength
    nmax:   maximum considered code length
    nist:   targeted NIST security level
    Mmax:   maximum memory that can be used by solver, deactivated if < 0

    outputs:
    params: valid parameter sets
    """
    assert nist in [1,3,5],"nist has to be 1, 3 or 5"
    if   nist == 1: 
        seclvl = 143
        kmin   =  20
    elif nist == 3: 
        seclvl = 207
        kmin   = 29
    elif nist == 5: 
        seclvl = 272
        kmin   = 38

    params =[]



    for n in range(nmin, nmax+1):
        print('doing n=',n)
        for k in range(kmin,n):
            for m in range(ceil(k/2),n):

                C_opt, _, _ = submatrix_stern(p, z, n, k, m, nist, Mmax)

                if C_opt >= seclvl:
                    print('found n=',n,'k=',k,'m=',m,' with C=',C_opt)
                    params += [(n,k,m)]
    return params

# field sizes
p = 509
z = 127

print("\n"
"#######################\n"
"# PARAMS for NIST I   #\n"
"#######################\n")

n_I    =  55
k_I    =  36
m_I    =  25
nist_I =   1
Mmax_I =  -1

C_I, M_I, P_I = submatrix_stern(p, z, n_I, k_I, m_I, nist_I, Mmax_I)

print(f'time complexity submatrix stern >= {rts(C_I)}bit')
print(f'memory cost submatrix stern     >= {rts(M_I)}bit')

print("\n"
"#######################\n"
"# PARAMS for NIST III #\n"
"#######################\n")

n_III    =  79
k_III    =  48
m_III    =  40
nist_III =   3
Mmax_III =  -1

C_III, M_III, P_III = submatrix_stern(p, z, n_III, k_III, m_III, nist_III, Mmax_III)

print(f'time complexity submatrix stern >= {rts(C_III)}bit')
print(f'memory cost submatrix stern     >= {rts(M_III)}bit')

print("\n"
"#######################\n"
"# PARAMS for NIST V   #\n"
"#######################\n")

n_V    =  106
k_V    =   69
m_V    =   48
nist_V =    5
Mmax_V =   -1

C_V, M_V, P_V = submatrix_stern(p, z, n_V, k_V, m_V, nist_V, Mmax_V)

print(f'time complexity submatrix stern >= {rts(C_V)}bit')
print(f'memory cost submatrix stern     >= {rts(M_V)}bit')