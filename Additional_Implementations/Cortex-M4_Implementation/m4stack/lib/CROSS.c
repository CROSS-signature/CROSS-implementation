/**
 *
 * Cortex-M4 Optimized ISO-C11 Implementation of CROSS.
 *
 * @version 2.2 (July 2025)
 *
 * Authors listed in alphabetical order:
 * 
 * @author: Alessandro Barenghi <alessandro.barenghi@polimi.it>
 * @author: Marco Gianvecchio <marco.gianvecchio@mail.polimi.it>
 * @author: Patrick Karl <patrick.karl@tum.de>
 * @author: Gerardo Pelosi <gerardo.pelosi@polimi.it>
 * @author: Jonas Schupp <jonas.schupp@tum.de>
 * 
 * 
 * This code is hereby placed in the public domain.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 **/

#include <assert.h>
#include <stdalign.h>
#include <stdio.h>
#include <stdbool.h>

#include "randombytes.h"
#include "m4_defs.h"
#include "CROSS.h"
#include "csprng_hash.h"
#include "fp_arith.h"
#include "merkle_tree.h"
#include "pack_unpack.h"
#include "seedtree.h"

#if defined(MEM_OPT_MTREE_IN_CCM) || defined(MEM_OPT_STREE_IN_CCM)
#define CCRAM __attribute__((section(".ccm")))
#endif

#if defined(MEM_OPT_MTREE_IN_CCM)
#if defined(MEM_OPT_HASH_TO_MTREE) && !defined(NO_TREES)
    uint8_t CCRAM mtree[NUM_NODES_MERKLE_TREE * HASH_DIGEST_LENGTH];
#else
    uint8_t CCRAM cmt_0[T][HASH_DIGEST_LENGTH] = {0};
#endif
#endif
#if defined(MEM_OPT_STREE_IN_CCM)
#if defined(NO_TREES)
/* #if !defined(MEM_OPT_VERIFY_INC_STREE) */
unsigned char CCRAM round_seeds[T*SEED_LENGTH_BYTES];
/* #endif */
#else
uint8_t CCRAM seed_tree[SEED_LENGTH_BYTES*NUM_NODES_SEED_TREE];
#if !defined(MEM_OPT_SIGN_BS_SEEDS_FROM_TREE) || !defined(MEM_OPT_VERIFY_BS_SEEDS_FROM_TREE)
unsigned char CCRAM round_seeds[T*SEED_LENGTH_BYTES];
#endif
#endif
#endif

// Please note that these structs might contain some padding bytes as we put instances of
// e.g. uint64_t behind an uint8_t array of arbitrary length.
// It should therefore be avoided to access members of these structs directly by byte offset.
// Margins for randomness are calculated based on the assumption that no padding is in place
// thus having sufficient space in that case.
struct V_tr_csprng_e_bar {
    alignas(M4_REG_BYTES) FP_ELEM_M4 V_tr[V_ROWS][V_COLS];
    CSPRNG_STATE_T csprng_state;
    FZ_ELEM e_bar[N];
#if defined(TIME_OPT_SIMD)
#if defined(RSDPG)
    // In the NO-SIMD case, the available buffer using the entire struct is always large enough. In the SIMD case, we have to expand V_tr transposed and can thus not easily store the randonmess in place. Again for RSDP, this is irrelevant as we expand the matrix as uint8_t and realign it afterwards leading to the same requirements for the random buffer as in the NON-SIMD case. For RSDPG though where the elements of V_tr are directly uint16_t, we need additional space for the random buffer in keygen and verify.
    uint8_t buffer[(BITS_V_CT_RNG+7)/8 - 208 - N];
#endif
#endif
};

struct sign_struct {
    alignas(M4_REG_BYTES) FP_ELEM_M4 V_tr[V_ROWS][V_COLS];
    CSPRNG_STATE_T csprng_state;
#if defined(MEM_OPT_SIGN_RECOMP_E_V_U)
    FZ_ELEM e_bar[N];
    FZ_ELEM e_bar_prime[N];
    FZ_ELEM v_bar[N];
    FP_ELEM u_prime[N];
#if defined(RSDPG)
    alignas(M4_REG_BYTES) FZ_ELEM_M4 e_G_bar_prime[M];
    alignas(M4_REG_BYTES) FZ_ELEM_M4 v_G_bar[M];
#endif
#endif
    alignas(M4_REG_BYTES) FP_ELEM_M4 v_u_shared[N];
    FP_ELEM s_prime[N-K];
#if defined(RSDP)
    uint8_t cmt_0_i_input[DENSELY_PACKED_FP_SYN_SIZE+
                          DENSELY_PACKED_FZ_VEC_SIZE+
                          SALT_LENGTH_BYTES];
#elif defined(RSDPG)
    uint8_t cmt_0_i_input[DENSELY_PACKED_FP_SYN_SIZE+
                          DENSELY_PACKED_FZ_RSDP_G_VEC_SIZE+
                          SALT_LENGTH_BYTES];
#endif
#if defined(MEM_OPT_SIGN_INC_MTREE)
#if defined(NO_TREES)
    CSPRNG_STATE_T cmt_0_inc_state;
#endif
#endif
#if defined(MEM_OPT_SIGN_INC_CMT_1) || defined(MEM_OPT_SIGN_RECOMP_Y)
    CSPRNG_STATE_T cmt_1_inc_state;
#endif
#if (defined(MEM_OPT_SIGN_INC_MTREE) || defined(MEM_OPT_SIGN_INC_STREE))
#if !defined(NO_TREES)
    uint8_t flag_tree[NUM_NODES_MERKLE_TREE];
#endif
#endif
    FP_ELEM chall_1[T];
    uint8_t chall_2[T];
// In the following we calculate the worst case length of our buffers as we need to extend it in certain cases/want to be informed at compile time if things don't work out. This is only relevant for in place sampling when we use some other variable spaces for the random buffers we need. We ignore additional bytes from alignment for all cases as they only increase the available buffer size but might not be always available in the same way. Omitting them is here the same as going for a worst-case assumption on the available buffer space.
// We need the length of the flag tree depending on whether we have it or not
#if defined(NO_TREES)
#define FLAG_TREE_LENGTH 0
#else
#define FLAG_TREE_LENGTH NUM_NODES_MERKLE_TREE
#endif
#if !defined(TIME_OPT_SIMD)
// Case 1: No SIMD
#if defined(RSDP)
#define BUFFER_LENGTH (BITS_V_CT_RNG+7)/8 - ((N-K)*K*7)/8 - 208*3 - N*6 - (N-K) - (DENSELY_PACKED_FP_SYN_SIZE+DENSELY_PACKED_FZ_VEC_SIZE+SALT_LENGTH_BYTES) - FLAG_TREE_LENGTH - T*2
#elif defined(RSDPG)
#define BUFFER_LENGTH (BITS_V_CT_RNG+7)/8 - ((N-K)*K*9)/8 - 208*3 - N*7 - M*2*2 - (N-K)*2 - (DENSELY_PACKED_FP_SYN_SIZE+DENSELY_PACKED_FZ_VEC_SIZE+SALT_LENGTH_BYTES) - FLAG_TREE_LENGTH - T*3
#endif
#elif defined(TIME_OPT_SIMD)
// Case 2: SIMD
#if defined(RSDP)
//For RSDP we allocate the random buffer for the matrix in the upper half of the matrix space, sample the matrix as uint8_t in the lower half and expand it afterwards by copying it once into uint16_t for SIMD. The remaining buffer (i.e. everything after V_tr needs therefore to be large enough to hold the random bits exceeding (N-K)*K bytes
#define BUFFER_LENGTH (BITS_V_CT_RNG+7)/8 - (N-K)*K - 208*3 - N*6 - (N-K) - (DENSELY_PACKED_FP_SYN_SIZE+DENSELY_PACKED_FZ_VEC_SIZE+SALT_LENGTH_BYTES) - FLAG_TREE_LENGTH - T*2
#elif defined(RSDPG)
    // The critical component here is the number of random bits necessary for V_tr sampling is RSDPG
#define BUFFER_LENGTH (BITS_V_CT_RNG+7)/8 - 208*3 - N*7 - M*2*2 - (N-K)*2 - (DENSELY_PACKED_FP_SYN_SIZE+DENSELY_PACKED_FZ_RSDP_G_VEC_SIZE+SALT_LENGTH_BYTES) - FLAG_TREE_LENGTH - T*3
#endif
#endif
#if BUFFER_LENGTH > 0
    uint8_t buffer[BUFFER_LENGTH];
#endif
};

#if defined(RSDP)
static
#if defined(MEM_OPT_IN_PLACE_SAMP)
void expand_pk(FP_ELEM_M4 V_tr[V_ROWS][V_COLS],
               const uint8_t seed_pk[KEYPAIR_SEED_LENGTH_BYTES],
               uint8_t *buffer){
#else
void expand_pk(FP_ELEM_M4 V_tr[V_ROWS][V_COLS],
               const uint8_t seed_pk[KEYPAIR_SEED_LENGTH_BYTES]){
#endif

  /* Expansion of pk->seed, explicit domain separation for CSPRNG as in keygen */
  const uint16_t dsc_csprng_seed_pk = CSPRNG_DOMAIN_SEP_CONST + (3*T+2);

#if defined(MEM_OPT_IN_PLACE_SAMP)
  CSPRNG_STATE_T *csprng_state_mat = (CSPRNG_STATE_T *)buffer;
  csprng_initialize(csprng_state_mat, seed_pk, KEYPAIR_SEED_LENGTH_BYTES, dsc_csprng_seed_pk);
#if !defined(TIME_OPT_SIMD)
#define BUFFER_MAT_OFFSET K*(N-K) - K*(N-K)*7/8
#else
#define BUFFER_MAT_OFFSET K*(N-K)
#endif
  // The following line will generate a compiler warning as we cast an uint8_t pointer on a uint16_t matrix.
  // This is intended as we sample in uint8_t and expand the matrix later on.
  csprng_fp_mat(V_tr, csprng_state_mat, buffer + BUFFER_MAT_OFFSET);
#else
  CSPRNG_STATE_T csprng_state_mat;
  csprng_initialize(&csprng_state_mat, seed_pk, KEYPAIR_SEED_LENGTH_BYTES, dsc_csprng_seed_pk);
  csprng_fp_mat(V_tr, &csprng_state_mat);
#endif
}
#elif defined(RSDPG)
static
#if defined(MEM_OPT_IN_PLACE_SAMP)
void expand_pk(FP_ELEM_M4 V_tr[V_ROWS][V_COLS],
               FZ_ELEM_M4 W_mat[W_ROWS][W_COLS],
               const uint8_t seed_pk[KEYPAIR_SEED_LENGTH_BYTES],
               uint8_t *buffer){
#else
void expand_pk(FP_ELEM_M4 V_tr[V_ROWS][V_COLS],
               FZ_ELEM_M4 W_mat[W_ROWS][W_COLS],
               const uint8_t seed_pk[KEYPAIR_SEED_LENGTH_BYTES]){
#endif

  /* Expansion of pk->seed, explicit domain separation for CSPRNG as in keygen */
  const uint16_t dsc_csprng_seed_pk = CSPRNG_DOMAIN_SEP_CONST + (3*T+2);

#if defined(MEM_OPT_IN_PLACE_SAMP)
  CSPRNG_STATE_T *csprng_state_mat = (CSPRNG_STATE_T *)buffer;
  csprng_initialize(csprng_state_mat, seed_pk, KEYPAIR_SEED_LENGTH_BYTES, dsc_csprng_seed_pk);
#else
  CSPRNG_STATE_T csprng_state_mat;
  csprng_initialize(&csprng_state_mat, seed_pk, KEYPAIR_SEED_LENGTH_BYTES, dsc_csprng_seed_pk);
#endif

#if defined(MEM_OPT_IN_PLACE_SAMP)
#if !defined(TIME_OPT_SIMD)
#define BUFFER_MAT_OFFSET K*(N-K)*2 - K*(N-K)*9/8
#define BUFFER_MAT_OFFSET_W BUFFER_MAT_OFFSET
#else
#define BUFFER_MAT_OFFSET K*(N-K)*2
#define BUFFER_MAT_OFFSET_W 208
#endif
  csprng_fz_mat(W_mat, csprng_state_mat, buffer + BUFFER_MAT_OFFSET_W);
  csprng_fp_mat(V_tr, csprng_state_mat, buffer + BUFFER_MAT_OFFSET);
#else
  csprng_fz_mat(W_mat, &csprng_state_mat);
  csprng_fp_mat(V_tr, &csprng_state_mat);
#endif
}
#endif


#if defined(RSDP)
static
#if defined(MEM_OPT_IN_PLACE_SAMP)
void expand_sk(FZ_ELEM e_bar[N],
               FP_ELEM_M4 V_tr[V_ROWS][V_COLS],
               const uint8_t seed_sk[KEYPAIR_SEED_LENGTH_BYTES],
               uint8_t *buffer){
#else
void expand_sk(FZ_ELEM e_bar[N],
               FP_ELEM_M4 V_tr[V_ROWS][V_COLS],
               const uint8_t seed_sk[KEYPAIR_SEED_LENGTH_BYTES]){
#endif

  uint8_t seed_e_seed_pk[2][KEYPAIR_SEED_LENGTH_BYTES];

  /* Expansion of sk->seed, explicit domain separation for CSPRNG, as in keygen */
  const uint16_t dsc_csprng_seed_sk = CSPRNG_DOMAIN_SEP_CONST + (3*T+1);

#if defined(MEM_OPT_IN_PLACE_SAMP)
  CSPRNG_STATE_T *csprng_state = &(((struct V_tr_csprng_e_bar *)buffer)->csprng_state);
  csprng_initialize(csprng_state, seed_sk, KEYPAIR_SEED_LENGTH_BYTES, dsc_csprng_seed_sk);
  csprng_randombytes((uint8_t *)seed_e_seed_pk,
                     2*KEYPAIR_SEED_LENGTH_BYTES,
                     csprng_state);
#else
  CSPRNG_STATE_T csprng_state;
  csprng_initialize(&csprng_state, seed_sk, KEYPAIR_SEED_LENGTH_BYTES, dsc_csprng_seed_sk);
  csprng_randombytes((uint8_t *)seed_e_seed_pk,
                     2*KEYPAIR_SEED_LENGTH_BYTES,
                     &csprng_state);
#endif

#if defined(MEM_OPT_IN_PLACE_SAMP)
  expand_pk(V_tr,seed_e_seed_pk[1], buffer);
#else
  expand_pk(V_tr,seed_e_seed_pk[1]);
#endif

  /* Expansion of seede, explicit domain separation for CSPRNG as in keygen */
  const uint16_t dsc_csprng_seed_e = CSPRNG_DOMAIN_SEP_CONST + (3*T+3);

#if defined(MEM_OPT_IN_PLACE_SAMP)
  csprng_initialize(csprng_state, seed_e_seed_pk[0], KEYPAIR_SEED_LENGTH_BYTES, dsc_csprng_seed_e);
  csprng_fz_vec(e_bar,csprng_state);
#else
  csprng_initialize(&csprng_state, seed_e_seed_pk[0], KEYPAIR_SEED_LENGTH_BYTES, dsc_csprng_seed_e);
  csprng_fz_vec(e_bar,&csprng_state);
#endif
}
#elif defined(RSDPG)
static
#if defined(MEM_OPT_IN_PLACE_SAMP)
void expand_sk(FZ_ELEM e_bar[N],
               FZ_ELEM_M4 e_G_bar[M],
               FP_ELEM_M4 V_tr[V_ROWS][V_COLS],
               FZ_ELEM_M4 W_mat[W_ROWS][W_COLS],
               const uint8_t seed_sk[KEYPAIR_SEED_LENGTH_BYTES],
               uint8_t *buffer){
#else
void expand_sk(FZ_ELEM e_bar[N],
               FZ_ELEM_M4 e_G_bar[M],
               FP_ELEM_M4 V_tr[V_ROWS][V_COLS],
               FZ_ELEM_M4 W_mat[W_ROWS][W_COLS],
               const uint8_t seed_sk[KEYPAIR_SEED_LENGTH_BYTES]){
#endif

  uint8_t seed_e_seed_pk[2][KEYPAIR_SEED_LENGTH_BYTES];

  /* Expansion of sk->seed, explicit domain separation for CSPRNG, as in keygen */
  const uint16_t dsc_csprng_seed_sk = CSPRNG_DOMAIN_SEP_CONST + (3*T+1);

#if defined(MEM_OPT_IN_PLACE_SAMP)
  CSPRNG_STATE_T *csprng_state = &(((struct V_tr_csprng_e_bar *)buffer)->csprng_state);
  csprng_initialize(csprng_state, seed_sk, KEYPAIR_SEED_LENGTH_BYTES, dsc_csprng_seed_sk);
  csprng_randombytes((uint8_t *)seed_e_seed_pk,
                     2*KEYPAIR_SEED_LENGTH_BYTES,
                     csprng_state);

#else
  CSPRNG_STATE_T csprng_state;

  csprng_initialize(&csprng_state, seed_sk, KEYPAIR_SEED_LENGTH_BYTES, dsc_csprng_seed_sk);
  csprng_randombytes((uint8_t *)seed_e_seed_pk,
                     2*KEYPAIR_SEED_LENGTH_BYTES,
                     &csprng_state);
#endif

#if defined(MEM_OPT_IN_PLACE_SAMP)
  expand_pk(V_tr,W_mat,seed_e_seed_pk[1], buffer);
#else
  expand_pk(V_tr,W_mat,seed_e_seed_pk[1]);
#endif

  /* Expansion of seede, explicit domain separation for CSPRNG as in keygen */
  const uint16_t dsc_csprng_seed_e = CSPRNG_DOMAIN_SEP_CONST + (3*T+3);

#if defined(MEM_OPT_IN_PLACE_SAMP)
  csprng_initialize(csprng_state, seed_e_seed_pk[0], KEYPAIR_SEED_LENGTH_BYTES, dsc_csprng_seed_e);
  csprng_fz_inf_w(e_G_bar,csprng_state);
#else
  csprng_initialize(&csprng_state, seed_e_seed_pk[0], KEYPAIR_SEED_LENGTH_BYTES, dsc_csprng_seed_e);
  csprng_fz_inf_w(e_G_bar,&csprng_state);
#endif
  fz_inf_w_by_fz_matrix(e_bar,e_G_bar,W_mat);
  fz_dz_norm_n(e_bar);
}
#endif

void CROSS_keygen(sk_t *SK,
                  pk_t *PK){

  struct V_tr_csprng_e_bar V_tr_csprng_e_bar_inst;
  /* generation of random material for public and private key */
  randombytes(SK->seed_sk,KEYPAIR_SEED_LENGTH_BYTES);

  uint8_t seed_e_seed_pk[2][KEYPAIR_SEED_LENGTH_BYTES];

  /* Expansion of sk->seed, explicit domain separation for CSPRNG */
  const uint16_t dsc_csprng_seed_sk = CSPRNG_DOMAIN_SEP_CONST + (3*T+1);
  
  /* csprng_state moved to matrix to use as random buffer */
  CSPRNG_STATE_T *csprng_state = &(V_tr_csprng_e_bar_inst.csprng_state);
  csprng_initialize(csprng_state, SK->seed_sk, KEYPAIR_SEED_LENGTH_BYTES, dsc_csprng_seed_sk);
  csprng_randombytes((uint8_t *)seed_e_seed_pk,
                     2*KEYPAIR_SEED_LENGTH_BYTES,
                     csprng_state);
  memcpy(PK->seed_pk,seed_e_seed_pk[1],KEYPAIR_SEED_LENGTH_BYTES);

  /* expansion of matrix/matrices */
  /* Matrix moved to struct to use as buffer for pseudorandomness in combination w. other arrays */
  FP_ELEM_M4 (* V_tr)[V_COLS] = V_tr_csprng_e_bar_inst.V_tr;
#if defined(RSDP)
#if defined(MEM_OPT_IN_PLACE_SAMP)
  expand_pk(V_tr,PK->seed_pk,(uint8_t *)&V_tr_csprng_e_bar_inst);
#else
  expand_pk(V_tr,PK->seed_pk);
#endif
#elif defined(RSDPG)
  alignas(M4_REG_BYTES) FZ_ELEM_M4 W_mat[W_ROWS][W_COLS];
#if defined(MEM_OPT_IN_PLACE_SAMP)
  expand_pk(V_tr,W_mat,PK->seed_pk,(uint8_t *)&V_tr_csprng_e_bar_inst);
#else
  expand_pk(V_tr,W_mat,PK->seed_pk);
#endif
#endif

  /* expansion of secret key material */ 
  /* Expansion of seede, explicit domain separation for CSPRNG */
  const uint16_t dsc_csprng_seed_e = CSPRNG_DOMAIN_SEP_CONST + (3*T+3);

  /* Reuse csprng state above */
  csprng_initialize(csprng_state, seed_e_seed_pk[0], KEYPAIR_SEED_LENGTH_BYTES, dsc_csprng_seed_e);

  FZ_ELEM *e_bar = (V_tr_csprng_e_bar_inst.e_bar);
#if defined(RSDP)
  csprng_fz_vec(e_bar,csprng_state);
#elif defined(RSDPG)
  alignas(M4_REG_BYTES) FZ_ELEM_M4 e_G_bar[M];
  csprng_fz_inf_w(e_G_bar,csprng_state);
  fz_inf_w_by_fz_matrix(e_bar,e_G_bar,W_mat);
  fz_dz_norm_n(e_bar);
#endif
  /* compute public syndrome */
  FP_ELEM s[N-K];
  restr_vec_by_fp_matrix(s,e_bar,V_tr);
  fp_dz_norm_synd(s);
  pack_fp_syn(PK->s,s);
}

/* sign cannot fail */
void CROSS_sign(const sk_t *SK,
               const char *const m,
               const uint64_t mlen,
               CROSS_sig_t *sig){
    
    struct sign_struct sign_struct_i;
    /* Key material expansion */
    FP_ELEM_M4 (* V_tr)[V_COLS] = sign_struct_i.V_tr;
#if defined(MEM_OPT_SIGN_RECOMP_E_V_U)
    FZ_ELEM *e_bar = sign_struct_i.e_bar;
#else
    FZ_ELEM e_bar[N];
#endif
#if defined(RSDP)
#if defined(MEM_OPT_IN_PLACE_SAMP)
    expand_sk(e_bar,V_tr,SK->seed_sk, (uint8_t *)&sign_struct_i);
#else
    expand_sk(e_bar,V_tr,SK->seed_sk);
#endif
#elif defined(RSDPG)
    alignas(M4_REG_BYTES) FZ_ELEM_M4 e_G_bar[M];
    alignas(M4_REG_BYTES) FZ_ELEM_M4 W_mat[W_ROWS][W_COLS];
#if defined(MEM_OPT_IN_PLACE_SAMP)
    expand_sk(e_bar,e_G_bar,V_tr,W_mat,SK->seed_sk, (uint8_t *)&sign_struct_i);
#else
    expand_sk(e_bar,e_G_bar,V_tr,W_mat,SK->seed_sk);
#endif
#endif

    /* Wipe any residual information in the sig structure allocated by the 
     * caller */
    memset(sig,0,sizeof(CROSS_sig_t));

    uint8_t root_seed[SEED_LENGTH_BYTES];
    randombytes(root_seed,SEED_LENGTH_BYTES);
    randombytes(sig->salt,SALT_LENGTH_BYTES);

#if defined(NO_TREES)
#if defined(MEM_OPT_SIGN_INC_STREE)
    CSPRNG_STATE_T csprng_state_seeds;
    uint8_t pos_in_quad = 0;
    unsigned char quad_seed[4*SEED_LENGTH_BYTES] = {0};
    seed_leaves_intermed(quad_seed, root_seed, sig->salt, &csprng_state_seeds);
#else
#if !defined(MEM_OPT_STREE_IN_CCM)
    unsigned char round_seeds[T*SEED_LENGTH_BYTES] = {0};
#endif
    seed_leaves(round_seeds,root_seed,sig->salt);
#endif
#else
#if !defined(MEM_OPT_SIGN_INC_STREE)
#if !defined(MEM_OPT_STREE_IN_CCM)
    uint8_t seed_tree[SEED_LENGTH_BYTES*NUM_NODES_SEED_TREE] = {0};
#if !defined(MEM_OPT_SIGN_BS_SEEDS_FROM_TREE)
    unsigned char round_seeds[T*SEED_LENGTH_BYTES] = {0};
#endif
#endif
    gen_seed_tree(seed_tree,root_seed,sig->salt);
#if !defined(MEM_OPT_SIGN_BS_SEEDS_FROM_TREE)
    seed_leaves(round_seeds, seed_tree);
#endif
#endif
#endif
#if defined(MEM_OPT_SIGN_INC_STREE)
#if !defined(NO_TREES)
    uint8_t stree_inc_seed_pos = 0;
    uint8_t directly_from_sig = 0;
    unsigned char seed_tree_inc[2*LOG2(T)*SEED_LENGTH_BYTES] = {0};
    uint16_t seed_tree_inc_pos[2*LOG2(T)] = {0};
    init_seed_tree_inc(seed_tree_inc, seed_tree_inc_pos, root_seed, sig->salt);
#endif
#endif
#if defined(MEM_OPT_SIGN_RECOMP_E_V_U)
    FZ_ELEM *e_bar_prime = sign_struct_i.e_bar_prime;
    FZ_ELEM *v_bar = sign_struct_i.v_bar;
    FP_ELEM *u_prime = sign_struct_i.u_prime;
#else
    FZ_ELEM e_bar_prime[T][N];
    FZ_ELEM v_bar[T][N];
    FP_ELEM u_prime[T][N];
#endif

#if defined(RSDP)
#define CMT_0_I_INPUT_SIZE DENSELY_PACKED_FP_SYN_SIZE + DENSELY_PACKED_FZ_VEC_SIZE + SALT_LENGTH_BYTES
    uint8_t *cmt_0_i_input = sign_struct_i.cmt_0_i_input;
    const int offset_salt = DENSELY_PACKED_FP_SYN_SIZE+DENSELY_PACKED_FZ_VEC_SIZE;
#elif defined(RSDPG)
#if defined(MEM_OPT_SIGN_RECOMP_E_V_U)
    FZ_ELEM_M4 *e_G_bar_prime = sign_struct_i.e_G_bar_prime;
    FZ_ELEM_M4 *v_G_bar = sign_struct_i.v_G_bar;
#else
    alignas(M4_REG_BYTES) FZ_ELEM_M4 e_G_bar_prime[M];
    alignas(M4_REG_BYTES) FZ_ELEM_M4 v_G_bar[T][M];
#endif
#define CMT_0_I_INPUT_SIZE DENSELY_PACKED_FP_SYN_SIZE + DENSELY_PACKED_FZ_RSDP_G_VEC_SIZE + SALT_LENGTH_BYTES
    uint8_t *cmt_0_i_input = sign_struct_i.cmt_0_i_input;
    const int offset_salt = DENSELY_PACKED_FP_SYN_SIZE+DENSELY_PACKED_FZ_RSDP_G_VEC_SIZE;
#endif
    /* cmt_0_i_input is syndrome || v_bar resp. v_G_bar || salt ; place salt at the end */
    memcpy(cmt_0_i_input+offset_salt, sig->salt, SALT_LENGTH_BYTES);

    uint8_t cmt_1_i_input[SEED_LENGTH_BYTES+
                          SALT_LENGTH_BYTES];
    /* cmt_1_i_input is concat(seed,salt,round index + 2T-1) */
    memcpy(cmt_1_i_input+SEED_LENGTH_BYTES, sig->salt, SALT_LENGTH_BYTES);

#if defined(MEM_OPT_SIGN_INC_MTREE) || defined(MEM_OPT_SIGN_INC_STREE)
#if !defined(NO_TREES)
    uint8_t *flag_tree = sign_struct_i.flag_tree;
#endif
#endif
#if defined(MEM_OPT_SIGN_INC_MTREE)
#if !defined(NO_TREES)
    uint8_t mtree_inc_target_address;
    uint8_t mtree_inc[(LOG2(T) + 2) * HASH_DIGEST_LENGTH];
    uint8_t mtree_inc_flags[LOG2(T) + 1] = {false};
    bool recurse_through_tree = false;
#elif defined(NO_TREES)
    CSPRNG_STATE_T *cmt_0_inc_state = &(sign_struct_i.cmt_0_inc_state);
    uint8_t cmt_0_intermed[4*HASH_DIGEST_LENGTH] = {0};
    uint8_t cmt_0_pos = 0;
#endif
#endif
#if !defined(MEM_OPT_MTREE_IN_CCM) && !defined(MEM_OPT_SIGN_INC_MTREE)
#if defined(MEM_OPT_HASH_TO_MTREE) && !defined(NO_TREES)
    uint8_t mtree[NUM_NODES_MERKLE_TREE * HASH_DIGEST_LENGTH];
#else
    uint8_t cmt_0[T][HASH_DIGEST_LENGTH] = {0};
#endif
#endif
#if defined(MEM_OPT_SIGN_INC_CMT_1) || defined(MEM_OPT_SIGN_RECOMP_Y)
    CSPRNG_STATE_T *cmt_1_inc_state = &(sign_struct_i.cmt_1_inc_state);
#endif
#if defined(MEM_OPT_SIGN_INC_CMT_1)
    hash_init(cmt_1_inc_state);
    uint8_t cmt_1[HASH_DIGEST_LENGTH] = {0};
#else
    uint8_t cmt_1[T*HASH_DIGEST_LENGTH] = {0};
#endif

    /* CSPRNG is fed with concat(seed,salt,round index) represented
        * as a 2 bytes little endian unsigned integer */
    uint8_t csprng_input[SEED_LENGTH_BYTES+SALT_LENGTH_BYTES];
    CSPRNG_STATE_T *csprng_state = &(sign_struct_i.csprng_state);

#if defined(MEM_OPT_HASH_TO_MTREE) || defined(MEM_OPT_SIGN_INC_MTREE) || defined(MEM_OPT_SIGN_BS_SEEDS_FROM_TREE) || defined(MEM_OPT_SIGN_INC_STREE)
#if !defined(NO_TREES)
    uint16_t pos_in_tree = 0;
    uint8_t curr_level = 0;
    uint8_t curr_si = 0;
    uint16_t curr_node_in_level = 0;
#endif
#endif

#if defined(TIME_OPT_SIMD)
#if defined(RSDP)
    /* u and y have different sizes, don't reuse memory */
    FP_ELEM v[N];
    alignas(M4_REG_BYTES) FP_ELEM_M4 u[N];
#else
    /* Reuse memory for v, and u */
    FP_ELEM_M4 *v_u_shared = sign_struct_i.v_u_shared;
    FP_ELEM *v = (FP_ELEM *) v_u_shared;
    FP_ELEM_M4 *u = (FP_ELEM_M4 *) v_u_shared;
#endif
#else
    /* Reuse memory for v, and u */
    FP_ELEM_M4 *v_u_shared = sign_struct_i.v_u_shared;
    FP_ELEM *v = (FP_ELEM *) v_u_shared;
    FP_ELEM_M4 *u = (FP_ELEM_M4 *) v_u_shared;
#endif

    FP_ELEM *s_prime = sign_struct_i.s_prime;

#if defined(MEM_OPT_SIGN_RECOMP_Y)
    /* y vectors are packed before being hashed */
    uint8_t y_digest_chall_1[DENSELY_PACKED_FP_VEC_SIZE];
    FP_ELEM *y = v;
#else
    FP_ELEM y[T][N];
#endif

    for(uint16_t i = 0; i<T; i++){

#if defined(MEM_OPT_HASH_TO_MTREE) || defined(MEM_OPT_SIGN_INC_MTREE) || defined(MEM_OPT_SIGN_BS_SEEDS_FROM_TREE) || defined(MEM_OPT_SIGN_INC_STREE)
#if !defined(NO_TREES)
        get_curr_pos_in_tree(&pos_in_tree, &curr_level, &curr_si, &curr_node_in_level);
#endif
#endif
#if defined(MEM_OPT_SIGN_INC_STREE)
#if !defined(NO_TREES)
        rebuild_tree_inc(&stree_inc_seed_pos,
                         &directly_from_sig,
                         seed_tree_inc,
                         flag_tree,
                         0,
                         seed_tree_inc_pos,
                         &pos_in_tree,
                         &curr_level,
                         sig->path,
                         sig->salt,
                         1);
        memcpy(csprng_input,
               seed_tree_inc + SEED_LENGTH_BYTES * stree_inc_seed_pos,
               SEED_LENGTH_BYTES);
#elif defined(NO_TREES)
        seed_leaves_inc(csprng_input, &csprng_state_seeds, quad_seed, &pos_in_quad, i, sig->salt);
#endif
#elif defined(MEM_OPT_SIGN_BS_SEEDS_FROM_TREE) && !defined(NO_TREES)
        memcpy(csprng_input,seed_tree+SEED_LENGTH_BYTES*pos_in_tree,SEED_LENGTH_BYTES);
#else
        memcpy(csprng_input,round_seeds+SEED_LENGTH_BYTES*i,SEED_LENGTH_BYTES);
#endif
        memcpy(csprng_input+SEED_LENGTH_BYTES,sig->salt,SALT_LENGTH_BYTES);

        uint16_t domain_sep_csprng = CSPRNG_DOMAIN_SEP_CONST+i+(2*T-1);

        /* expand seed[i] into seed_e and seed_u */
        csprng_initialize(csprng_state,
                          csprng_input,
                          SEED_LENGTH_BYTES+SALT_LENGTH_BYTES,
                          domain_sep_csprng);
        /* expand e_bar_prime */
#if defined(RSDP)
#if defined(MEM_OPT_SIGN_RECOMP_E_V_U)
        csprng_fz_vec(e_bar_prime, csprng_state);
#else
        csprng_fz_vec(e_bar_prime[i], csprng_state);
#endif
#elif defined(RSDPG)
        csprng_fz_inf_w(e_G_bar_prime, csprng_state);
#if defined(MEM_OPT_SIGN_RECOMP_E_V_U)
        fz_vec_sub_m(v_G_bar, e_G_bar, e_G_bar_prime);
        fz_dz_norm_m(v_G_bar);
        fz_inf_w_by_fz_matrix(e_bar_prime, e_G_bar_prime, W_mat);
        fz_dz_norm_n(e_bar_prime);
#else
        fz_vec_sub_m(v_G_bar[i], e_G_bar, e_G_bar_prime);
        fz_dz_norm_m(v_G_bar[i]);
        fz_inf_w_by_fz_matrix(e_bar_prime[i], e_G_bar_prime, W_mat);
        fz_dz_norm_n(e_bar_prime[i]);
#endif
#endif
#if defined(MEM_OPT_SIGN_RECOMP_E_V_U)
        fz_vec_sub_n(v_bar, e_bar, e_bar_prime);
#else
        fz_vec_sub_n(v_bar[i], e_bar, e_bar_prime[i]);
#endif

#if defined(MEM_OPT_SIGN_RECOMP_E_V_U)
        convert_restr_vec_to_fp(v, v_bar);
        fz_dz_norm_n(v_bar);
        /* expand u_prime */
        csprng_fp_vec(u_prime, csprng_state);
        fp_vec_by_fp_vec_pointwise(u, v, u_prime);
#else
        convert_restr_vec_to_fp(v, v_bar[i]);
        fz_dz_norm_n(v_bar[i]);
        /* expand u_prime */
        csprng_fp_vec(u_prime[i], csprng_state);
        fp_vec_by_fp_vec_pointwise(u, v, u_prime[i]);
#endif

        fp_vec_by_fp_matrix(s_prime, u, V_tr);
        fp_dz_norm_synd(s_prime);

        /* cmt_0_i_input contains s_prime || v_bar resp. v_G_bar || salt */
        pack_fp_syn(cmt_0_i_input,s_prime);

#if defined(RSDP)
#if defined(MEM_OPT_SIGN_RECOMP_E_V_U)
        pack_fz_vec(cmt_0_i_input + DENSELY_PACKED_FP_SYN_SIZE, v_bar);
#else
        pack_fz_vec(cmt_0_i_input + DENSELY_PACKED_FP_SYN_SIZE, v_bar[i]);
#endif
#elif defined(RSDPG)
#if defined(MEM_OPT_SIGN_RECOMP_E_V_U)
        pack_fz_rsdp_g_vec(cmt_0_i_input + DENSELY_PACKED_FP_SYN_SIZE, v_G_bar);
#else
        pack_fz_rsdp_g_vec(cmt_0_i_input + DENSELY_PACKED_FP_SYN_SIZE, v_G_bar[i]);
#endif
#endif
        /* Fixed endianness marshalling of round counter */
        uint16_t domain_sep_hash = HASH_DOMAIN_SEP_CONST+i+(2*T-1);

#if defined(MEM_OPT_SIGN_INC_MTREE)
#if !defined(NO_TREES)
        recompute_root_mtree_inc_recomp(&mtree_inc_target_address,
                                        &recurse_through_tree,
                                        mtree_inc_flags,
                                        &curr_level);
        hash(mtree_inc + mtree_inc_target_address * HASH_DIGEST_LENGTH, cmt_0_i_input, CMT_0_I_INPUT_SIZE, domain_sep_hash);
        if(recurse_through_tree == true) {
            recompute_root_mtree_inc(mtree_inc,
                                     &recurse_through_tree,
                                     mtree_inc_flags);
        }
#elif defined(NO_TREES)
        hash(cmt_0_intermed+cmt_0_pos*HASH_DIGEST_LENGTH, cmt_0_i_input, CMT_0_I_INPUT_SIZE, domain_sep_hash);
        tree_intermed_inc(cmt_0_inc_state, cmt_0_intermed, &cmt_0_pos, cmt_0_intermed+cmt_0_pos*HASH_DIGEST_LENGTH, i);
#endif
#elif defined(MEM_OPT_HASH_TO_MTREE) && !defined(NO_TREES)
        hash(mtree + pos_in_tree * HASH_DIGEST_LENGTH, cmt_0_i_input, CMT_0_I_INPUT_SIZE, domain_sep_hash);
#else
        hash(cmt_0[i], cmt_0_i_input, CMT_0_I_INPUT_SIZE, domain_sep_hash);
#endif
#if defined(MEM_OPT_SIGN_INC_STREE)
#if !defined(NO_TREES)
        memcpy(cmt_1_i_input,
               seed_tree_inc + SEED_LENGTH_BYTES * stree_inc_seed_pos,
               SEED_LENGTH_BYTES);
#elif defined(NO_TREES)
        memcpy(cmt_1_i_input,
               csprng_input,
               SEED_LENGTH_BYTES);
#endif
#elif defined(MEM_OPT_SIGN_BS_SEEDS_FROM_TREE) && !defined(NO_TREES)
        memcpy(cmt_1_i_input,
               seed_tree+SEED_LENGTH_BYTES*pos_in_tree,
               SEED_LENGTH_BYTES);
#else
        memcpy(cmt_1_i_input,
               round_seeds+SEED_LENGTH_BYTES*i,
               SEED_LENGTH_BYTES);
#endif
#if defined(MEM_OPT_SIGN_INC_CMT_1) 
        hash(cmt_1, cmt_1_i_input, sizeof(cmt_1_i_input), domain_sep_hash);
        hash_update(cmt_1_inc_state, cmt_1, sizeof(cmt_1));
#else
        hash(&cmt_1[i*HASH_DIGEST_LENGTH], cmt_1_i_input, sizeof(cmt_1_i_input), domain_sep_hash);
#endif
    }

    /* vector containing d_0 and d_1 from spec */
    uint8_t digest_cmt0_cmt1[2*HASH_DIGEST_LENGTH];

#if defined(MEM_OPT_SIGN_INC_MTREE)
#if !defined(NO_TREES)
    memcpy(digest_cmt0_cmt1, mtree_inc, HASH_DIGEST_LENGTH);
#elif defined(NO_TREES)
    hash(digest_cmt0_cmt1, cmt_0_intermed, 4*HASH_DIGEST_LENGTH, HASH_DOMAIN_SEP_CONST);
#endif
#elif !defined(MEM_OPT_SIGN_INC_MTREE) && defined(NO_TREES)
    tree_root(digest_cmt0_cmt1, cmt_0);
#else
#if !defined(MEM_OPT_HASH_TO_MTREE)
    uint8_t mtree[NUM_NODES_MERKLE_TREE * HASH_DIGEST_LENGTH];
    tree_root(digest_cmt0_cmt1, mtree, cmt_0);
#else
    tree_root(digest_cmt0_cmt1, mtree);
#endif
#endif
#if defined(MEM_OPT_SIGN_INC_CMT_1) 
    hash_fin_and_squeeze(cmt_1_inc_state, digest_cmt0_cmt1 + HASH_DIGEST_LENGTH, HASH_DOMAIN_SEP_CONST);
#else
    hash(digest_cmt0_cmt1 + HASH_DIGEST_LENGTH, cmt_1, sizeof(cmt_1), HASH_DOMAIN_SEP_CONST);
#endif
    hash(sig->digest_cmt, digest_cmt0_cmt1, sizeof(digest_cmt0_cmt1), HASH_DOMAIN_SEP_CONST);

    /* first challenge extraction */
    uint8_t digest_msg_cmt_salt[2*HASH_DIGEST_LENGTH+SALT_LENGTH_BYTES];

    /* place digest_msg at the beginning of the input of the hash generating digest_chall_1 */
    hash(digest_msg_cmt_salt, (uint8_t*) m, mlen, HASH_DOMAIN_SEP_CONST);
    memcpy(digest_msg_cmt_salt+HASH_DIGEST_LENGTH, sig->digest_cmt, HASH_DIGEST_LENGTH);
    memcpy(digest_msg_cmt_salt+2*HASH_DIGEST_LENGTH, sig->salt, SALT_LENGTH_BYTES);

    uint8_t digest_chall_1[HASH_DIGEST_LENGTH];
    hash(digest_chall_1, digest_msg_cmt_salt, sizeof(digest_msg_cmt_salt), HASH_DOMAIN_SEP_CONST);

    // Domain separation unique for expanding chall_1
    const uint16_t dsc_csprng_chall_1 = CSPRNG_DOMAIN_SEP_CONST + (3*T-1);

    FP_ELEM *chall_1 = sign_struct_i.chall_1;
    csprng_initialize(csprng_state, digest_chall_1, sizeof(digest_chall_1), dsc_csprng_chall_1);
#if defined(MEM_OPT_IN_PLACE_SAMP)
#if defined(MEM_OPT_SIGN_RECOMP_E_V_U)
    csprng_fp_vec_chall_1(chall_1, csprng_state, (uint8_t *) sign_struct_i.e_bar_prime);
#else
    uint8_t aux_buffer[(BITS_CWSTR_RNG+7)/8];
    csprng_fp_vec_chall_1(chall_1, csprng_state, aux_buffer);
#endif
#else
    csprng_fp_vec_chall_1(chall_1, csprng_state);
#endif

#if defined(MEM_OPT_SIGN_RECOMP_Y)
    hash_init(cmt_1_inc_state);
#endif
#if (defined(MEM_OPT_SIGN_BS_SEEDS_FROM_TREE) || defined(MEM_OPT_SIGN_INC_STREE))
#if !defined(NO_TREES)
    pos_in_tree = 0;
    curr_level = 0;
    curr_si = 0;
    curr_node_in_level = 0;
#endif
#endif
#if defined(MEM_OPT_SIGN_INC_STREE)
#if !defined(NO_TREES)
    init_seed_tree_inc(seed_tree_inc, seed_tree_inc_pos, root_seed, sig->salt);
#endif
#endif
#if defined(MEM_OPT_SIGN_INC_STREE)
#if defined(NO_TREES)
    pos_in_quad = 0;
#endif
#endif

    for(int i = 0; i < T; i++){
#if defined(MEM_OPT_SIGN_RECOMP_E_V_U)
#if defined(MEM_OPT_SIGN_BS_SEEDS_FROM_TREE) || defined(MEM_OPT_SIGN_INC_STREE)
#if !defined(NO_TREES)
        get_curr_pos_in_tree(&pos_in_tree, &curr_level, &curr_si, &curr_node_in_level);
#endif
#endif
#if defined(MEM_OPT_SIGN_INC_STREE)
#if !defined(NO_TREES)
        rebuild_tree_inc(&stree_inc_seed_pos,
                         &directly_from_sig,
                         seed_tree_inc,
                         flag_tree,
                         0,
                         seed_tree_inc_pos,
                         &pos_in_tree,
                         &curr_level,
                         sig->path,
                         sig->salt,
                         1);
        memcpy(csprng_input,
               seed_tree_inc + SEED_LENGTH_BYTES * stree_inc_seed_pos,
               SEED_LENGTH_BYTES);
#elif defined(NO_TREES)
        seed_leaves_inc(csprng_input, &csprng_state_seeds, quad_seed, &pos_in_quad, i, sig->salt);
#endif
#elif defined(MEM_OPT_SIGN_BS_SEEDS_FROM_TREE) && !defined(NO_TREES)
        memcpy(csprng_input,seed_tree+SEED_LENGTH_BYTES*pos_in_tree,SEED_LENGTH_BYTES);
#else
        memcpy(csprng_input,round_seeds+SEED_LENGTH_BYTES*i,SEED_LENGTH_BYTES);
#endif
        memcpy(csprng_input+SEED_LENGTH_BYTES,sig->salt,SALT_LENGTH_BYTES);

        uint16_t domain_sep_csprng = CSPRNG_DOMAIN_SEP_CONST+i+(2*T-1);
        /* expand seed[i] into seed_e and seed_u */
        csprng_initialize(csprng_state,
                          csprng_input,
                          SEED_LENGTH_BYTES+SALT_LENGTH_BYTES,
                          domain_sep_csprng);
        /* expand e_bar_prime */
#if defined(RSDP)
        csprng_fz_vec(e_bar_prime, csprng_state);
#elif defined(RSDPG)
        csprng_fz_inf_w(e_G_bar_prime, csprng_state);
        fz_vec_sub_m(v_G_bar, e_G_bar, e_G_bar_prime);
        fz_dz_norm_m(v_G_bar);
        fz_inf_w_by_fz_matrix(e_bar_prime, e_G_bar_prime, W_mat);
        fz_dz_norm_n(e_bar_prime);
#endif
        /* expand u_prime */
        csprng_fp_vec(u_prime, csprng_state);
#endif
#if defined(MEM_OPT_SIGN_RECOMP_E_V_U) && defined(MEM_OPT_SIGN_RECOMP_Y)
        fp_vec_by_restr_vec_scaled(y,
                                   e_bar_prime,
                                   chall_1[i],
                                   u_prime);
        fp_dz_norm(y);
        pack_fp_vec(y_digest_chall_1, y);
        hash_update(cmt_1_inc_state, y_digest_chall_1, sizeof(y_digest_chall_1));
#elif defined(MEM_OPT_SIGN_RECOMP_E_V_U)
        fp_vec_by_restr_vec_scaled(y[i],
                                   e_bar_prime,
                                   chall_1[i],
                                   u_prime);
        fp_dz_norm(y[i]);
#elif defined(MEM_OPT_SIGN_RECOMP_Y)
        fp_vec_by_restr_vec_scaled(y,
                                   e_bar_prime[i],
                                   chall_1[i],
                                   u_prime[i]);
        fp_dz_norm(y);
        pack_fp_vec(y_digest_chall_1, y);
        hash_update(cmt_1_inc_state, y_digest_chall_1, sizeof(y_digest_chall_1));
#else
        fp_vec_by_restr_vec_scaled(y[i],
                                   e_bar_prime[i],
                                   chall_1[i],
                                   u_prime[i]);
        fp_dz_norm(y[i]);
#endif
    }
#if defined(MEM_OPT_SIGN_RECOMP_Y)
    /* Second challenge extraction */
    hash_update(cmt_1_inc_state, digest_chall_1, HASH_DIGEST_LENGTH);
    hash_fin_and_squeeze(cmt_1_inc_state, sig->digest_chall_2, HASH_DOMAIN_SEP_CONST);
#else
    /* y vectors are packed before being hashed */
    uint8_t y_digest_chall_1[T*DENSELY_PACKED_FP_VEC_SIZE+HASH_DIGEST_LENGTH];

    for(int x = 0; x < T; x++){
        pack_fp_vec(y_digest_chall_1+(x*DENSELY_PACKED_FP_VEC_SIZE),y[x]);
    }
    /* Second challenge extraction */
    memcpy(y_digest_chall_1+T*DENSELY_PACKED_FP_VEC_SIZE,digest_chall_1,HASH_DIGEST_LENGTH);

    hash(sig->digest_chall_2, y_digest_chall_1, sizeof(y_digest_chall_1), HASH_DOMAIN_SEP_CONST);
#endif

    uint8_t *chall_2 = sign_struct_i.chall_2;
#if defined(MEM_OPT_IN_PLACE_SAMP)
#if defined(MEM_OPT_SIGN_RECOMP_E_V_U)
    expand_digest_to_fixed_weight(chall_2,sig->digest_chall_2, csprng_state, (uint8_t *) sign_struct_i.e_bar_prime);
#else
    expand_digest_to_fixed_weight(chall_2,sig->digest_chall_2, csprng_state, aux_buffer);
#endif
#else
    expand_digest_to_fixed_weight(chall_2,sig->digest_chall_2);
#endif

    /* Computation of the second round of responses */
#if defined(NO_TREES)
#if !defined(MEM_OPT_SIGN_INC_MTREE)
    tree_proof(sig->proof,cmt_0,chall_2);
#endif
#if !defined(MEM_OPT_SIGN_INC_STREE)
    seed_path(sig->path,round_seeds,chall_2);
#endif
#else
#if !defined(MEM_OPT_SIGN_INC_MTREE)
    tree_proof(sig->proof,mtree,chall_2);
#endif
#if !defined(MEM_OPT_SIGN_INC_STREE)
    seed_path(sig->path,seed_tree,chall_2);
#endif
#endif

#if (defined(MEM_OPT_HASH_TO_MTREE) || defined(MEM_OPT_SIGN_INC_MTREE) || defined(MEM_OPT_SIGN_BS_SEEDS_FROM_TREE) || defined(MEM_OPT_SIGN_INC_STREE))
#if !defined(NO_TREES)
    pos_in_tree = 0;
    curr_level = 0;
    curr_si = 0;
    curr_node_in_level = 0;
#endif
#endif

#if defined(NO_TREES)
#if defined(MEM_OPT_SIGN_INC_STREE)
    pos_in_quad = 0;
#endif
#if (defined(MEM_OPT_SIGN_INC_STREE) || defined(MEM_OPT_SIGN_INC_MTREE))
    uint16_t pos_in_proof = 0;
#endif
#endif

#if !defined(NO_TREES)
#if defined(MEM_OPT_SIGN_INC_MTREE) || defined(MEM_OPT_SIGN_INC_STREE)
    uint8_t total_published;
    memset(flag_tree, 0, NUM_NODES_MERKLE_TREE);
    total_published = recompute_flag_tree(flag_tree, chall_2);
#endif
#if defined(MEM_OPT_SIGN_INC_MTREE)
    for(int i = 0; i < LOG2(T)+1; i++){
        mtree_inc_flags[i] = 0;
    }
    recurse_through_tree = false;
#endif
#if defined(MEM_OPT_SIGN_INC_STREE)
    init_seed_tree_inc(seed_tree_inc, seed_tree_inc_pos, root_seed, sig->salt);
    uint16_t last_copied = 0;
#endif
#endif

    int published_rsps = 0;
    for(int i = 0; i<T; i++){
#if (defined(MEM_OPT_HASH_TO_MTREE) || defined(MEM_OPT_SIGN_INC_MTREE) || defined(MEM_OPT_SIGN_BS_SEEDS_FROM_TREE) || defined(MEM_OPT_SIGN_INC_STREE))
#if !defined(NO_TREES)
        get_curr_pos_in_tree(&pos_in_tree, &curr_level, &curr_si, &curr_node_in_level);
#endif
#endif
#if !(defined(MEM_OPT_SIGN_INC_MTREE) || defined(MEM_OPT_SIGN_INC_STREE))
        if(chall_2[i] == 0){
            assert(published_rsps < T-W);
#else
        {
#endif
            assert(published_rsps <= T-W);
#if defined(MEM_OPT_SIGN_INC_STREE)
#if !defined(NO_TREES)
            rebuild_tree_inc(&stree_inc_seed_pos,
                             &directly_from_sig,
                             seed_tree_inc,
                             flag_tree,
                             0,
                             seed_tree_inc_pos,
                             &pos_in_tree,
                             &curr_level,
                             sig->path,
                             sig->salt,
                             1);
            memcpy(csprng_input,
                   seed_tree_inc + SEED_LENGTH_BYTES * stree_inc_seed_pos,
                   SEED_LENGTH_BYTES);
#elif defined(NO_TREES)
            seed_leaves_inc(csprng_input, &csprng_state_seeds, quad_seed, &pos_in_quad, i, sig->salt);
#endif
#elif defined(MEM_OPT_SIGN_BS_SEEDS_FROM_TREE) && !defined(NO_TREES)
            memcpy(csprng_input,seed_tree+SEED_LENGTH_BYTES*pos_in_tree,SEED_LENGTH_BYTES);
#else
            memcpy(csprng_input,round_seeds+SEED_LENGTH_BYTES*i,SEED_LENGTH_BYTES);
#endif
#if defined(MEM_OPT_SIGN_RECOMP_E_V_U)
            memcpy(csprng_input+SEED_LENGTH_BYTES,sig->salt,SALT_LENGTH_BYTES);

            uint16_t domain_sep_csprng = CSPRNG_DOMAIN_SEP_CONST+i+(2*T-1);

            /* expand seed[i] into seed_e and seed_u */
            csprng_initialize(csprng_state,
                              csprng_input,
                              SEED_LENGTH_BYTES+SALT_LENGTH_BYTES,
                              domain_sep_csprng);
            /* expand e_bar_prime */
#if defined(RSDP)
            csprng_fz_vec(e_bar_prime, csprng_state);
            fz_vec_sub_n(v_bar, e_bar, e_bar_prime);
            fz_dz_norm_n(v_bar);
            if(chall_2[i] == 0){
                pack_fz_vec(sig->resp_0[published_rsps].v_bar, v_bar);
            }
#elif defined(RSDPG)
            csprng_fz_inf_w(e_G_bar_prime, csprng_state);
            fz_vec_sub_m(v_G_bar, e_G_bar, e_G_bar_prime);
            fz_dz_norm_m(v_G_bar);
            fz_inf_w_by_fz_matrix(e_bar_prime, e_G_bar_prime, W_mat);
            fz_dz_norm_n(e_bar_prime);
            fz_vec_sub_n(v_bar, e_bar, e_bar_prime);

            if(chall_2[i] == 0){
                pack_fz_rsdp_g_vec(sig->resp_0[published_rsps].v_G_bar, v_G_bar);
            }
#endif
#else
            if(chall_2[i] == 0){
#if defined(RSDP)
                pack_fz_vec(sig->resp_0[published_rsps].v_bar, v_bar[i]);
#elif defined(RSDPG)
                pack_fz_rsdp_g_vec(sig->resp_0[published_rsps].v_G_bar, v_G_bar[i]);
#endif
            }
#endif
#if defined(MEM_OPT_SIGN_RECOMP_E_V_U) && (defined(MEM_OPT_SIGN_RECOMP_Y) || defined(MEM_OPT_SIGN_INC_MTREE))
            /* expand u_prime */
            csprng_fp_vec(u_prime, csprng_state);
#endif
#if defined(MEM_OPT_SIGN_RECOMP_E_V_U) && defined(MEM_OPT_SIGN_RECOMP_Y)
            fp_vec_by_restr_vec_scaled(y,
                                   e_bar_prime,
                                   chall_1[i],
                                   u_prime);
            fp_dz_norm(y);
            if(chall_2[i] == 0){
                pack_fp_vec(sig->resp_0[published_rsps].y, y);
            }
#elif defined(MEM_OPT_SIGN_RECOMP_Y)
            fp_vec_by_restr_vec_scaled(y,
                                   e_bar_prime[i],
                                   chall_1[i],
                                   u_prime[i]);
            fp_dz_norm(y);
            if(chall_2[i] == 0){
                pack_fp_vec(sig->resp_0[published_rsps].y, y);
            }
#else
            if(chall_2[i] == 0){
                pack_fp_vec(sig->resp_0[published_rsps].y, y[i]);
            }
#endif
#if defined(MEM_OPT_SIGN_INC_CMT_1) || defined(MEM_OPT_SIGN_INC_MTREE)
            /* Fixed endianness marshalling of round counter */
            uint16_t domain_sep_hash = HASH_DOMAIN_SEP_CONST+i+(2*T-1);
#endif
#if defined(MEM_OPT_SIGN_INC_CMT_1) 
#if defined(MEM_OPT_SIGN_INC_STREE)
#if !defined(NO_TREES)
            memcpy(cmt_1_i_input,
                   seed_tree_inc + SEED_LENGTH_BYTES * stree_inc_seed_pos,
                   SEED_LENGTH_BYTES);
#elif defined(NO_TREES)
            memcpy(cmt_1_i_input,
                   csprng_input,
                   SEED_LENGTH_BYTES);
#endif
#elif defined(MEM_OPT_SIGN_BS_SEEDS_FROM_TREE) && !defined(NO_TREES)
            memcpy(cmt_1_i_input,
                seed_tree+SEED_LENGTH_BYTES*pos_in_tree,
                SEED_LENGTH_BYTES);
#else
            memcpy(cmt_1_i_input,
                round_seeds+SEED_LENGTH_BYTES*i,
                SEED_LENGTH_BYTES);
#endif
            if(chall_2[i] == 0){
                hash(sig->resp_1[published_rsps], cmt_1_i_input, sizeof(cmt_1_i_input), domain_sep_hash);
                published_rsps++;
            }
#else
            if(chall_2[i] == 0){
                memcpy(sig->resp_1[published_rsps], &cmt_1[i*HASH_DIGEST_LENGTH], HASH_DIGEST_LENGTH);
                published_rsps++;
            }
#endif
        
#if defined(MEM_OPT_SIGN_RECOMP_E_V_U) && defined(MEM_OPT_SIGN_INC_MTREE)
            convert_restr_vec_to_fp(v, v_bar);
            fp_vec_by_fp_vec_pointwise(u, v, u_prime);
            fp_vec_by_fp_matrix(s_prime, u, V_tr);
            fp_dz_norm_synd(s_prime);

            /* cmt_0_i_input is syndrome || v_bar resp. v_G_bar || salt ; place salt at the end */
            memcpy(cmt_0_i_input+offset_salt, sig->salt, SALT_LENGTH_BYTES);
            /* cmt_0_i_input contains s_prime || v_bar resp. v_G_bar || salt */
            pack_fp_syn(cmt_0_i_input,s_prime);
#if defined(RSDP)
            pack_fz_vec(cmt_0_i_input + DENSELY_PACKED_FP_SYN_SIZE, v_bar);
#elif defined(RSDPG)
            pack_fz_rsdp_g_vec(cmt_0_i_input + DENSELY_PACKED_FP_SYN_SIZE, v_G_bar);
#endif
#elif !defined(MEM_OPT_SIGN_RECOMP_E_V_U) && defined(MEM_OPT_SIGN_INC_MTREE)
            convert_restr_vec_to_fp(v, v_bar[i]);
            fp_vec_by_fp_vec_pointwise(u, v, u_prime[i]);
            fp_vec_by_fp_matrix(s_prime, u, V_tr);
            fp_dz_norm_synd(s_prime);

            /* cmt_0_i_input contains s_prime || v_bar resp. v_G_bar || salt */
            pack_fp_syn(cmt_0_i_input,s_prime);
#if defined(RSDP)
            pack_fz_vec(cmt_0_i_input + DENSELY_PACKED_FP_SYN_SIZE, v_bar[i]);
#elif defined(RSDPG)
            pack_fz_rsdp_g_vec(cmt_0_i_input + DENSELY_PACKED_FP_SYN_SIZE, v_G_bar[i]);
#endif
#endif
#if defined(NO_TREES)
#if defined(MEM_OPT_SIGN_INC_MTREE)
            if(chall_2[i] == 1){
                hash(sig->proof+HASH_DIGEST_LENGTH*pos_in_proof, cmt_0_i_input, CMT_0_I_INPUT_SIZE, domain_sep_hash);
            }
#endif
#if defined(MEM_OPT_SIGN_INC_STREE)
            if(chall_2[i] == 1){
                memcpy(sig->path+SEED_LENGTH_BYTES*pos_in_proof, csprng_input, SEED_LENGTH_BYTES);
            }
#endif
#endif
#if defined(MEM_OPT_SIGN_INC_MTREE) || defined(MEM_OPT_SIGN_INC_STREE)
#if defined(NO_TREES)
            if(chall_2[i] == 1){
                pos_in_proof++;
            }
#endif
#endif
#if !defined(NO_TREES)
#if defined(MEM_OPT_SIGN_INC_MTREE)
            recompute_root_mtree_inc_recomp(&mtree_inc_target_address,
                                            &recurse_through_tree,
                                            mtree_inc_flags,
                                            &curr_level);
            hash(mtree_inc + mtree_inc_target_address * HASH_DIGEST_LENGTH, cmt_0_i_input, CMT_0_I_INPUT_SIZE, domain_sep_hash);
            if(flag_tree[pos_in_tree] != 0){
                /* This is a node that we can directly copy, i.e. a leave that is part of the signature */
                memcpy(sig->proof + (flag_tree[pos_in_tree] - 1) * HASH_DIGEST_LENGTH, mtree_inc + mtree_inc_target_address * HASH_DIGEST_LENGTH, HASH_DIGEST_LENGTH);
            }
            if(recurse_through_tree == true) {
                mtree_inc_proof(mtree_inc,
                                sig->proof,
                                flag_tree,
                                &pos_in_tree,
                                &curr_level,
                                &recurse_through_tree,
                                mtree_inc_flags);
            }
#endif
#if defined(MEM_OPT_SIGN_INC_STREE)
            seed_path_inc(seed_tree_inc,
                          sig->path,
                          flag_tree,
                          total_published,
                          &last_copied,
                          seed_tree_inc_pos,
                          &pos_in_tree,
                          &curr_level);
#endif
#endif
        }
    }
}

/* verify returns 1 if signature is ok, 0 otherwise */
int CROSS_verify(const pk_t *const PK,
                 const char *const m,
                 const uint64_t mlen,
                 const CROSS_sig_t *const sig){

    struct V_tr_csprng_e_bar V_tr_csprng_e_bar_inst;
    CSPRNG_STATE_T *csprng_state = &(V_tr_csprng_e_bar_inst.csprng_state);

    FP_ELEM_M4 (* V_tr)[V_COLS] = V_tr_csprng_e_bar_inst.V_tr;

    FP_ELEM s[N-K];
    uint8_t is_padd_key_ok;
    is_padd_key_ok = unpack_fp_syn(s,PK->s);

    uint8_t digest_msg_cmt_salt[2*HASH_DIGEST_LENGTH+SALT_LENGTH_BYTES];
    hash(digest_msg_cmt_salt, (uint8_t*) m, mlen, HASH_DOMAIN_SEP_CONST);
    memcpy(digest_msg_cmt_salt+HASH_DIGEST_LENGTH, sig->digest_cmt, HASH_DIGEST_LENGTH);
    memcpy(digest_msg_cmt_salt+2*HASH_DIGEST_LENGTH, sig->salt, SALT_LENGTH_BYTES);

    uint8_t digest_chall_1[HASH_DIGEST_LENGTH];
    hash(digest_chall_1, digest_msg_cmt_salt, sizeof(digest_msg_cmt_salt), HASH_DOMAIN_SEP_CONST);

    // Domain separation unique for expanding digest_chall_1
    const uint16_t dsc_csprng_chall_1 = CSPRNG_DOMAIN_SEP_CONST + (3*T-1);
    csprng_initialize(csprng_state, digest_chall_1, sizeof(digest_chall_1), dsc_csprng_chall_1);

    FP_ELEM chall_1[T];
#if defined(MEM_OPT_IN_PLACE_SAMP)
    csprng_fp_vec_chall_1(chall_1, csprng_state, (uint8_t *)&V_tr_csprng_e_bar_inst);
#else
    csprng_fp_vec_chall_1(chall_1, csprng_state);
#endif

    uint8_t chall_2[T]={0};
#if defined(MEM_OPT_IN_PLACE_SAMP)
    expand_digest_to_fixed_weight(chall_2, sig->digest_chall_2, csprng_state, (uint8_t *)&V_tr_csprng_e_bar_inst);
#else
    expand_digest_to_fixed_weight(chall_2, sig->digest_chall_2);
#endif

#if defined(RSDP)
#if defined(MEM_OPT_IN_PLACE_SAMP)
    expand_pk(V_tr,PK->seed_pk,(uint8_t *)&V_tr_csprng_e_bar_inst);
#else
    expand_pk(V_tr,PK->seed_pk);
#endif
#elif defined(RSDPG)
    FZ_ELEM_M4 W_mat[W_ROWS][W_COLS];
#if defined(MEM_OPT_IN_PLACE_SAMP)
    expand_pk(V_tr,W_mat,PK->seed_pk,(uint8_t *)&V_tr_csprng_e_bar_inst);
#else
    expand_pk(V_tr,W_mat,PK->seed_pk);
#endif
#endif

    uint8_t is_stree_padding_ok = 0;
#if defined(NO_TREES)
#if defined(MEM_OPT_VERIFY_INC_STREE)
    is_stree_padding_ok = 1;
#else
#if !defined(MEM_OPT_STREE_IN_CCM)
    uint8_t round_seeds[T*SEED_LENGTH_BYTES] = {0};
#endif
    is_stree_padding_ok = rebuild_leaves(round_seeds, chall_2, sig->path);
#endif
#else
#if !defined(MEM_OPT_VERIFY_INC_STREE)
#if !defined(MEM_OPT_STREE_IN_CCM)
    uint8_t seed_tree[SEED_LENGTH_BYTES*NUM_NODES_SEED_TREE] = {0};
#if !defined(MEM_OPT_VERIFY_BS_SEEDS_FROM_TREE)
    unsigned char round_seeds[T*SEED_LENGTH_BYTES] = {0};
#endif
#endif
    is_stree_padding_ok = rebuild_tree(seed_tree, chall_2, sig->path, sig->salt);
#if !defined(MEM_OPT_VERIFY_BS_SEEDS_FROM_TREE)
    seed_leaves(round_seeds, seed_tree);
#endif
#endif
#endif

#if defined(RSDP)
    uint8_t cmt_0_i_input[DENSELY_PACKED_FP_SYN_SIZE+
                          DENSELY_PACKED_FZ_VEC_SIZE+
                          SALT_LENGTH_BYTES];
    const int offset_salt = DENSELY_PACKED_FP_SYN_SIZE+DENSELY_PACKED_FZ_VEC_SIZE;
#elif defined(RSDPG)
    uint8_t cmt_0_i_input[DENSELY_PACKED_FP_SYN_SIZE+
                          DENSELY_PACKED_FZ_RSDP_G_VEC_SIZE+
                          SALT_LENGTH_BYTES];
    const int offset_salt = DENSELY_PACKED_FP_SYN_SIZE+DENSELY_PACKED_FZ_RSDP_G_VEC_SIZE;
#endif
    /* cmt_0_i_input is syndrome || v_bar resp. v_G_bar || salt */
    memcpy(cmt_0_i_input+offset_salt, sig->salt, SALT_LENGTH_BYTES);

    /* cmt_1_i_input is concat(seed,salt,round index) */
    uint8_t cmt_1_i_input[SEED_LENGTH_BYTES+SALT_LENGTH_BYTES];
    memcpy(cmt_1_i_input+SEED_LENGTH_BYTES, sig->salt, SALT_LENGTH_BYTES);

#if defined(MEM_OPT_VERIFY_INC_MTREE) || defined(MEM_OPT_VERIFY_INC_STREE)
#if !defined(NO_TREES)
    uint8_t flag_tree[NUM_NODES_MERKLE_TREE] = {0};
    uint8_t num_published = recompute_flag_tree(flag_tree, chall_2);
#endif
#endif

#if defined(MEM_OPT_VERIFY_INC_MTREE)
#if !defined(NO_TREES)
    uint8_t mtree_inc_target_address;
    uint8_t mtree_inc[(LOG2(T) + 2) * HASH_DIGEST_LENGTH];
    uint8_t mtree_inc_flags[LOG2(T) + 1] = {false};
    bool recurse_through_tree = false;
    uint16_t flag_tree_address;
    uint8_t mtree_last_taken = 0;
#elif defined(NO_TREES)
    CSPRNG_STATE_T cmt_0_inc_state;
    uint8_t cmt_0_intermed[4*HASH_DIGEST_LENGTH] = {0};
    uint8_t cmt_0_pos = 0;
#endif
#else
#if !defined(MEM_OPT_MTREE_IN_CCM)
#if defined(MEM_OPT_HASH_TO_MTREE) && !defined(NO_TREES)
    uint8_t mtree[NUM_NODES_MERKLE_TREE * HASH_DIGEST_LENGTH];
#else
    uint8_t cmt_0[T][HASH_DIGEST_LENGTH] = {0};
#endif
#endif
#endif

#if defined(MEM_OPT_VERIFY_INC_STREE)
#if !defined(NO_TREES)
    is_stree_padding_ok = 1;
    for(int i = num_published; i < TREE_NODES_TO_STORE; i++){
        for(int j = 0; j < SEED_LENGTH_BYTES; j++){
            is_stree_padding_ok &= (sig->path[i*SEED_LENGTH_BYTES + j] == 0);
        }
    }
#endif
#endif

#ifdef MEM_OPT_VERIFY_INC_CMT_1
    CSPRNG_STATE_T cmt_1_inc_state;
    hash_init(&cmt_1_inc_state);
    uint8_t cmt_1[HASH_DIGEST_LENGTH] = {0};
#else
    uint8_t cmt_1[T*HASH_DIGEST_LENGTH] = {0};
#endif

    FZ_ELEM *e_bar_prime = V_tr_csprng_e_bar_inst.e_bar;

#if defined(TIME_OPT_SIMD)
#if defined(RSDP)
    FP_ELEM u_prime[N];
    alignas(M4_REG_BYTES) FP_ELEM_M4 y_prime[N];
#else
    alignas(M4_REG_BYTES) FP_ELEM_M4 u_y_shared[N];
    FP_ELEM *u_prime = (FP_ELEM *) u_y_shared;
    FP_ELEM_M4 *y_prime = (FP_ELEM_M4 *) u_y_shared;
#endif
#else
    alignas(M4_REG_BYTES) FP_ELEM_M4 u_y_shared[N];
    FP_ELEM *u_prime = (FP_ELEM *) u_y_shared;
    FP_ELEM_M4 *y_prime = (FP_ELEM_M4 *) u_y_shared;
#endif

    FP_ELEM y_prime_H[N-K] = {0};
    FP_ELEM *s_prime = y_prime_H;

#if defined(MEM_OPT_VERIFY_INC_Y)
    FP_ELEM *y = u_prime;
    CSPRNG_STATE_T y_digest_chall_1_state;
    hash_init(&y_digest_chall_1_state);
    uint8_t y_packed[DENSELY_PACKED_FP_VEC_SIZE];
#else
    FP_ELEM y[T][N];
#endif

#if defined(NO_TREES)
#if defined(MEM_OPT_VERIFY_INC_STREE)
    int used_seeds = 0;
#endif
#if defined(MEM_OPT_VERIFY_INC_MTREE)
    uint16_t pos_in_proof = 0;
#endif
#endif

#if !defined(NO_TREES)
#if defined(MEM_OPT_HASH_TO_MTREE) || defined(MEM_OPT_VERIFY_INC_MTREE) || defined(MEM_OPT_VERIFY_BS_SEEDS_FROM_TREE) || defined(MEM_OPT_VERIFY_INC_STREE)
    uint16_t pos_in_tree = 0;
    uint8_t curr_level = 0;
    uint8_t curr_si = 0;
    uint16_t curr_node_in_level = 0;
#endif
#if defined(MEM_OPT_VERIFY_INC_STREE)
    uint8_t stree_inc_seed_pos = 0;
    uint8_t directly_from_sig = 0;
    unsigned char seed_tree_inc[2*LOG2(T)*SEED_LENGTH_BYTES] = {0};
    uint16_t seed_tree_inc_pos[2*LOG2(T)] = {0};
#endif
#endif

    int used_rsps = 0;
    int is_signature_ok = 1;
    uint8_t is_packed_padd_ok = 1;

    for(uint16_t i = 0; i< T; i++){
        /* CSPRNG is fed with concat(seed,salt,round index) represented
        * as a 2 bytes little endian unsigned integer */
        const int csprng_input_length = SALT_LENGTH_BYTES+SEED_LENGTH_BYTES;
        uint8_t csprng_input[csprng_input_length];

#if defined(MEM_OPT_HASH_TO_MTREE) || defined(MEM_OPT_VERIFY_INC_MTREE) || defined(MEM_OPT_VERIFY_INC_STREE) || defined(MEM_OPT_VERIFY_BS_SEEDS_FROM_TREE)
#if !defined(NO_TREES)
        get_curr_pos_in_tree(&pos_in_tree, &curr_level, &curr_si, &curr_node_in_level);
#endif
#endif

        uint16_t domain_sep_csprng = CSPRNG_DOMAIN_SEP_CONST + i + (2*T-1);
        uint16_t domain_sep_hash = HASH_DOMAIN_SEP_CONST + i + (2*T-1);

        if(chall_2[i] == 1){

#if defined(MEM_OPT_VERIFY_INC_STREE)
#if !defined(NO_TREES)
        rebuild_tree_inc(&stree_inc_seed_pos,
                         &directly_from_sig,
                         seed_tree_inc,
                         flag_tree,
                         num_published,
                         seed_tree_inc_pos,
                         &pos_in_tree,
                         &curr_level,
                         sig->path,
                         sig->salt,
                         0);
#endif
#endif
#if defined(MEM_OPT_VERIFY_INC_STREE)
#if defined(NO_TREES)
            memcpy(cmt_1_i_input,
                   sig->path+used_seeds*SEED_LENGTH_BYTES,
                   SEED_LENGTH_BYTES);
#elif !defined(NO_TREES)
            uint8_t const *seed_ptr;
            if(directly_from_sig == 0){
                seed_ptr = seed_tree_inc + SEED_LENGTH_BYTES*stree_inc_seed_pos;
            }else{
                seed_ptr = sig->path + SEED_LENGTH_BYTES*stree_inc_seed_pos;
            }
            memcpy(cmt_1_i_input,
                   seed_ptr,
                   SEED_LENGTH_BYTES);
#endif
#elif defined(MEM_OPT_VERIFY_BS_SEEDS_FROM_TREE) && !defined(NO_TREES)
            memcpy(cmt_1_i_input,
                   seed_tree+SEED_LENGTH_BYTES*pos_in_tree,
                   SEED_LENGTH_BYTES);
#else
            memcpy(cmt_1_i_input,
                   round_seeds+SEED_LENGTH_BYTES*i,
                   SEED_LENGTH_BYTES);
#endif
#ifdef MEM_OPT_VERIFY_INC_CMT_1
            hash(cmt_1, cmt_1_i_input, sizeof(cmt_1_i_input), domain_sep_hash);
            hash_update(&cmt_1_inc_state, cmt_1, sizeof(cmt_1));
#else
            hash(&cmt_1[i*HASH_DIGEST_LENGTH],cmt_1_i_input,sizeof(cmt_1_i_input), domain_sep_hash);
#endif

            memcpy(csprng_input+SEED_LENGTH_BYTES,sig->salt,SALT_LENGTH_BYTES);
#if defined(MEM_OPT_VERIFY_INC_STREE)
#if defined(NO_TREES)
            memcpy(csprng_input,sig->path+used_seeds*SEED_LENGTH_BYTES,SEED_LENGTH_BYTES);
            used_seeds++;
#elif !defined(NO_TREES)
            memcpy(csprng_input,seed_ptr,SEED_LENGTH_BYTES);
#endif
#elif defined(MEM_OPT_VERIFY_BS_SEEDS_FROM_TREE) && !defined(NO_TREES)
            memcpy(csprng_input,seed_tree+SEED_LENGTH_BYTES*pos_in_tree,SEED_LENGTH_BYTES);
#else
            memcpy(csprng_input,round_seeds+SEED_LENGTH_BYTES*i,SEED_LENGTH_BYTES);
#endif

            /* expand seed[i] into seed_e and seed_u */
            csprng_initialize(csprng_state, csprng_input, csprng_input_length, domain_sep_csprng);
#if defined(RSDP)
            /* expand e_bar_prime */
            csprng_fz_vec(e_bar_prime, csprng_state);
#elif defined(RSDPG)
            alignas(M4_REG_BYTES) FZ_ELEM_M4 e_G_bar_prime[M];
            csprng_fz_inf_w(e_G_bar_prime, csprng_state);
            fz_inf_w_by_fz_matrix(e_bar_prime,e_G_bar_prime,W_mat);
            fz_dz_norm_n(e_bar_prime);
#endif
            /* expand u_prime */
            csprng_fp_vec(u_prime, csprng_state);
#if defined(MEM_OPT_VERIFY_INC_Y)
            fp_vec_by_restr_vec_scaled(y,
                                       e_bar_prime,
                                       chall_1[i],
                                       u_prime);
            fp_dz_norm(y);
            pack_fp_vec(y_packed, y);
            hash_update(&y_digest_chall_1_state, y_packed, DENSELY_PACKED_FP_VEC_SIZE);
#else
            fp_vec_by_restr_vec_scaled(y[i],
                                       e_bar_prime,
                                       chall_1[i],
                                       u_prime);
            fp_dz_norm(y[i]);
#endif
#if defined(MEM_OPT_VERIFY_INC_MTREE)
#if !defined(NO_TREES)
            recompute_root_mtree_inc_proof(&mtree_inc_target_address,
                                           &flag_tree_address,
                                           flag_tree,
                                           &mtree_last_taken,
                                           mtree_inc_flags,
                                           &pos_in_tree,
                                           &curr_level,
                                           &recurse_through_tree);
            if(mtree_inc_target_address != 0){
                memcpy(mtree_inc + mtree_inc_target_address * HASH_DIGEST_LENGTH, sig->proof + (flag_tree[flag_tree_address] - 1) * HASH_DIGEST_LENGTH, HASH_DIGEST_LENGTH);
            }
#elif defined(NO_TREES)
        memcpy(cmt_0_intermed+cmt_0_pos*HASH_DIGEST_LENGTH, sig->proof + HASH_DIGEST_LENGTH*pos_in_proof, HASH_DIGEST_LENGTH);
        pos_in_proof++;
#endif
#endif
        } else {
#if defined(MEM_OPT_VERIFY_INC_Y)
            is_packed_padd_ok = is_packed_padd_ok &&
                                unpack_fp_vec(y, sig->resp_0[used_rsps].y);
            hash_update(&y_digest_chall_1_state, sig->resp_0[used_rsps].y, DENSELY_PACKED_FP_VEC_SIZE);
#else
            /* place y[i] in the buffer for later on hashing */
            is_packed_padd_ok = is_packed_padd_ok &&
                                unpack_fp_vec(y[i], sig->resp_0[used_rsps].y);
#endif

            FZ_ELEM v_bar[N];
#if defined(RSDP)
            /*v_bar is memcpy'ed directly into cmt_0 input buffer */
            FZ_ELEM* v_bar_ptr = cmt_0_i_input+DENSELY_PACKED_FP_SYN_SIZE;
            is_packed_padd_ok = is_packed_padd_ok &&
                                unpack_fz_vec(v_bar, sig->resp_0[used_rsps].v_bar);
            memcpy(v_bar_ptr,
                   &sig->resp_0[used_rsps].v_bar,
                   DENSELY_PACKED_FZ_VEC_SIZE);
            is_signature_ok = is_signature_ok &&
                              is_fz_vec_in_restr_group_n(v_bar);
#elif defined(RSDPG)
            /*v_G_bar is memcpy'ed directly into cmt_0 input buffer */
            FZ_ELEM* v_G_bar_ptr = cmt_0_i_input+DENSELY_PACKED_FP_SYN_SIZE;
            memcpy(v_G_bar_ptr,
                   &sig->resp_0[used_rsps].v_G_bar,
                   DENSELY_PACKED_FZ_RSDP_G_VEC_SIZE);
            alignas(M4_REG_BYTES) FZ_ELEM_M4 v_G_bar[M];
            is_packed_padd_ok = is_packed_padd_ok &&
                                unpack_fz_rsdp_g_vec(v_G_bar, sig->resp_0[used_rsps].v_G_bar);
            is_signature_ok = is_signature_ok &&
                              is_fz_vec_in_restr_group_m(v_G_bar);
            fz_inf_w_by_fz_matrix(v_bar,v_G_bar,W_mat);

#endif
#ifdef MEM_OPT_VERIFY_INC_CMT_1
            hash_update(&cmt_1_inc_state, sig->resp_1[used_rsps], HASH_DIGEST_LENGTH);
#else
            memcpy(&cmt_1[i*HASH_DIGEST_LENGTH], sig->resp_1[used_rsps], HASH_DIGEST_LENGTH);
#endif

            used_rsps++;

            FP_ELEM v[N];
            convert_restr_vec_to_fp(v,v_bar);
#if defined(MEM_OPT_VERIFY_INC_Y)
            fp_vec_by_fp_vec_pointwise(y_prime,v,y);
#else
            fp_vec_by_fp_vec_pointwise(y_prime,v,y[i]);
#endif
            fp_vec_by_fp_matrix(y_prime_H,y_prime,V_tr);
            fp_dz_norm_synd(y_prime_H);
            fp_synd_minus_fp_vec_scaled(s_prime,
                                        y_prime_H,
                                        chall_1[i],
                                        s);
            fp_dz_norm_synd(s_prime);
            pack_fp_syn(cmt_0_i_input, s_prime);
#if defined(MEM_OPT_VERIFY_INC_MTREE)
#if !defined(NO_TREES)
            recompute_root_mtree_inc_recomp(&mtree_inc_target_address,
                                            &recurse_through_tree,
                                            mtree_inc_flags,
                                            &curr_level);
            hash(mtree_inc + mtree_inc_target_address * HASH_DIGEST_LENGTH, cmt_0_i_input, CMT_0_I_INPUT_SIZE, domain_sep_hash);
#elif defined(NO_TREES)
            hash(cmt_0_intermed+cmt_0_pos*HASH_DIGEST_LENGTH, cmt_0_i_input, CMT_0_I_INPUT_SIZE, domain_sep_hash);
#endif
#elif defined(MEM_OPT_HASH_TO_MTREE) && !defined(NO_TREES)
            hash(mtree + pos_in_tree * HASH_DIGEST_LENGTH, cmt_0_i_input, CMT_0_I_INPUT_SIZE, domain_sep_hash);
#else
            hash(cmt_0[i], cmt_0_i_input, CMT_0_I_INPUT_SIZE, domain_sep_hash);
#endif
        }
#if defined(MEM_OPT_VERIFY_INC_MTREE)
#if !defined(NO_TREES)
        if(recurse_through_tree == true) {
            recompute_root_mtree_inc(mtree_inc,
                                     &recurse_through_tree,
                                     mtree_inc_flags);
        }
#elif defined(NO_TREES)
        tree_intermed_inc(&cmt_0_inc_state, cmt_0_intermed, &cmt_0_pos, cmt_0_intermed+cmt_0_pos*HASH_DIGEST_LENGTH, i);
#endif
#endif
    } /* end for iterating on ZKID iterations */

    #ifndef SKIP_ASSERT
    assert(is_signature_ok);
    #endif

    uint8_t digest_cmt0_cmt1[2*HASH_DIGEST_LENGTH];

#if defined(MEM_OPT_VERIFY_INC_MTREE)
#if !defined(NO_TREES)
    memcpy(digest_cmt0_cmt1, mtree_inc, HASH_DIGEST_LENGTH);
    uint8_t is_mtree_padding_ok = 1;
    for(int i = num_published; i < TREE_NODES_TO_STORE; i++){
        for(int j = 0; j < HASH_DIGEST_LENGTH; j++){
            is_mtree_padding_ok &= (sig->proof[i*HASH_DIGEST_LENGTH + j] == 0);
        }
    }
#elif defined(NO_TREES)
    hash(digest_cmt0_cmt1, cmt_0_intermed, 4*HASH_DIGEST_LENGTH, HASH_DOMAIN_SEP_CONST);
    uint8_t is_mtree_padding_ok = 1;
#endif
#elif defined(MEM_OPT_HASH_TO_MTREE) && !defined(NO_TREES)
    uint8_t is_mtree_padding_ok = recompute_root(digest_cmt0_cmt1,
                                                 mtree,
                                                 sig->proof,
                                                 chall_2);
#else
    uint8_t is_mtree_padding_ok = recompute_root(digest_cmt0_cmt1,
                                                 cmt_0,
                                                 sig->proof,
                                                 chall_2);
#endif

#ifdef MEM_OPT_VERIFY_INC_CMT_1
    hash_fin_and_squeeze(&cmt_1_inc_state, digest_cmt0_cmt1 + HASH_DIGEST_LENGTH, HASH_DOMAIN_SEP_CONST);
#else
    hash(digest_cmt0_cmt1 + HASH_DIGEST_LENGTH, cmt_1, sizeof(cmt_1), HASH_DOMAIN_SEP_CONST);
#endif

    uint8_t digest_cmt_prime[HASH_DIGEST_LENGTH];
    hash(digest_cmt_prime, digest_cmt0_cmt1 ,sizeof(digest_cmt0_cmt1), HASH_DOMAIN_SEP_CONST);

#if defined(MEM_OPT_VERIFY_INC_Y)
    hash_update(&y_digest_chall_1_state, digest_chall_1, HASH_DIGEST_LENGTH);

    uint8_t digest_chall_2_prime[HASH_DIGEST_LENGTH];
    hash_fin_and_squeeze(&y_digest_chall_1_state, digest_chall_2_prime, HASH_DOMAIN_SEP_CONST);
#else
    uint8_t y_digest_chall_1[T*DENSELY_PACKED_FP_VEC_SIZE+HASH_DIGEST_LENGTH];

    for(int x = 0; x < T; x++){
        pack_fp_vec(y_digest_chall_1+(x*DENSELY_PACKED_FP_VEC_SIZE), y[x]);
    }
    memcpy(y_digest_chall_1+T*DENSELY_PACKED_FP_VEC_SIZE, digest_chall_1, HASH_DIGEST_LENGTH);

    uint8_t digest_chall_2_prime[HASH_DIGEST_LENGTH];
    hash(digest_chall_2_prime, y_digest_chall_1, sizeof(y_digest_chall_1), HASH_DOMAIN_SEP_CONST);
#endif

    int does_digest_cmt_match = ( memcmp(digest_cmt_prime,
                                        sig->digest_cmt,
                                        HASH_DIGEST_LENGTH) == 0);

    #ifndef SKIP_ASSERT
    assert(does_digest_cmt_match);
    #endif

    int does_digest_chall_2_match = ( memcmp(digest_chall_2_prime,
                                        sig->digest_chall_2,
                                        HASH_DIGEST_LENGTH) == 0);
    #ifndef SKIP_ASSERT
    assert(does_digest_chall_2_match);
    #endif

    is_signature_ok = is_signature_ok &&
                      does_digest_cmt_match &&
                      does_digest_chall_2_match && 
                      is_mtree_padding_ok &&
                      is_stree_padding_ok &&
                      is_padd_key_ok &&
                      is_packed_padd_ok;
    return is_signature_ok;
}
