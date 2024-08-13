# Created for SageMath version 9.0                  │
# Using Python 3.8.10

load('../macros.sage')
load('stern.sage')
load('stern_shifted.sage')
load('BJMM.sage')
load('BJMM_shifted.sage')


# choice of F_p and E #
p = 127
z = 7

# additive structure of E #
Eset   = get_Eset(p,z)
Dset   = get_Dset(Eset)
zd     = len(Dset)
alphaE = get_alphaE(Eset)
alphaD = get_alphaD(Eset, Dset)

# structure after shifting
Eshift = shift_E(Eset, Eset[0])
zshift = len(Eshift)
Dshift = get_Dset(Eshift)
alphaEshift = get_alphaE(Eshift)
alphaDshift = get_alphaD(Eshift, Dshift)

print(
"#######################\n"
"# PARAMS for NIST I   #\n"
"#######################\n")

n_I = 127
k_I =  76
nsolI = N(z**n_I * p**(k_I-n_I))+1
verb = False
MmaxI = -1 # no memory restriction

C_st_I,     M_st_I,     P_st_I      = opt_stern(p, z, n_I, k_I, MmaxI, verb)
C_shst_I,   M_shst_I,   P_shst_I    = opt_shifted_stern(p, zshift, n_I, k_I, MmaxI, verb)
C_bjmm_I,   M_bjmm_I,   P_bjmm_I    = opt_bjmm(p, z, n_I, k_I, MmaxI, verbose)
C_shbjmm_I, M_shbjmm_I, P_shbjmm_I  = opt_shifted_bjmm(p, zshift, n_I, k_I, MmaxI, verb)

print(f'C: Stern={rts(C_st_I)}bit; shifted Stern={rts(C_shst_I)}bit; BJMM={rts(C_bjmm_I)}bit; shifted BJMM={rts(C_shbjmm_I)}bit')
print(f'M: Stern={rts(M_st_I)}bit; shifted Stern={rts(M_shst_I)}bit; BJMM={rts(M_bjmm_I)}bit; shifted BJMM={rts(M_shbjmm_I)}bit')

print("\n"
"#######################\n"
"# PARAMS for NIST III #\n"
"#######################\n")

n_III = 187
k_III = 111
nsolIII = N(z**n_III * p**(k_III-n_III))+1
verb = False
MmaxIII = -1 # no memory restriction

C_st_III,     M_st_III,     P_st_III      = opt_stern(p, z, n_III, k_III, MmaxIII, verb)
C_shst_III,   M_shst_III,   P_shst_III    = opt_shifted_stern(p, zshift, n_III, k_III, MmaxIII, verb)
C_bjmm_III,   M_bjmm_III,   P_bjmm_III   = opt_bjmm(p, z, n_III, k_III, MmaxIII, verb)
C_shbjmm_III, M_shbjmm_III, P_shbjmm_III  = opt_shifted_bjmm(p, zshift, n_III, k_III, MmaxIII, verb)

print(f'C: Stern={rts(C_st_III)}bit; shifted Stern={rts(C_shst_III)}bit; BJMM={rts(C_bjmm_III)}bit; shifted BJMM={rts(C_shbjmm_III)}bit')
print(f'M: Stern={rts(M_st_III)}bit; shifted Stern={rts(M_shst_III)}bit; BJMM={rts(M_bjmm_III)}bit; shifted BJMM={rts(M_shbjmm_III)}bit')

print("\n"
"#######################\n"
"# PARAMS for NIST V   #\n"
"#######################\n")

n_V = 251
k_V = 150
nsolV = N(z**n_V * p**(k_V-n_V))+1
verb = False
MmaxV = -1

C_st_V,     M_st_V,     P_st_V      = opt_stern(p, z, n_V, k_V, MmaxV, verb)
C_shst_V,   M_shst_V,   P_shst_V    = opt_shifted_stern(p, zshift, n_V, k_V, MmaxV, verb)
C_bjmm_V,   M_bjmm_V,   P_bjmm_V    = opt_bjmm(p, z, n_V, k_V, MmaxV, verb)
C_shbjmm_V, M_shbjmm_V, P_shbjmm_V  = opt_shifted_bjmm(p, zshift, n_V, k_V, MmaxV, verb)

print(f'C: Stern={rts(C_st_V)}bit; shifted Stern={rts(C_shst_V)}bit; BJMM={rts(C_bjmm_V)}bit; shifted BJMM={rts(C_shbjmm_V)}bit')
print(f'M: Stern={rts(M_st_V)}bit; shifted Stern={rts(M_shst_V)}bit; BJMM={rts(M_bjmm_V)}bit; shifted BJMM={rts(M_shbjmm_V)}bit')
