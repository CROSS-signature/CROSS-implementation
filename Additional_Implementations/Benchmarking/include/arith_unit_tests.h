/**
 *
 * Reference ISO-C11 Implementation of CROSS.
 *
 * @version 2.0 (February 2025)
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

#pragma once



static inline
int trivial_is_restricted(FZ_ELEM* v, int len){
    int is_ok = 1;
    for(int i =0; i<len; i++){
        is_ok = is_ok ||
                ( (v[i] == 1) || (v[i] == 2) || (v[i] == 4) || (v[i] == 8) ||
                  (v[i] == 16) || (v[i] == 32) || (v[i] == 64));
    }
    return is_ok;
}

int restr_belonging_test(void){
    uint8_t seed[SEED_LENGTH_BYTES];
    randombytes(seed, SEED_LENGTH_BYTES);
    CSPRNG_STATE_T csprng_state;
    csprng_initialize(&csprng_state, seed, SEED_LENGTH_BYTES, 0);

#if defined(RSDP)
    FZ_ELEM a[N]={0};
    csprng_fz_vec(a, &csprng_state);
    return trivial_is_restricted(a,N) == is_fz_vec_in_restr_group_n(a);
#else
    FZ_ELEM a[M]={0};
    csprng_fz_inf_w(a,&csprng_state);
    return trivial_is_restricted(a,M) == is_fz_vec_in_restr_group_m(a);
#endif
}

void prettyprint(char* name, FP_ELEM *v, int l){
    fprintf(stderr,"%s: [",name);
    for(int i=0; i<l; i++){
        fprintf(stderr," %d, ", v[i]);
    }
    fprintf(stderr,"]\n");
}

#if defined(RSDP)
int signature_invariant_test(void){
    uint8_t seed[SEED_LENGTH_BYTES];
    randombytes(seed,SEED_LENGTH_BYTES);
    CSPRNG_STATE_T csprng_state;
    csprng_initialize(&csprng_state, seed, SEED_LENGTH_BYTES, 0);

    /* keygen */
    FZ_ELEM e_bar[N];
    FP_ELEM V_tr[K][N-K];
    csprng_fz_vec(e_bar, &csprng_state);
    csprng_fp_mat(V_tr,&csprng_state);
    FP_ELEM s[N-K] = {0};
    restr_vec_by_fp_matrix(s, e_bar, V_tr);
    fp_dz_norm_synd(s);

#if defined(HIGH_PERFORMANCE_X86_64)
    /* Expanded */
    alignas(EPI8_PER_REG) FP_DOUBLEPREC V_tr_avx[K][ROUND_UP(N-K,EPI16_PER_REG)] = {{0}};
    for(int i = 0; i < K; i++){
      for (int j = 0; j < N-K; j++){
         V_tr_avx[i][j] = V_tr[i][j];
      }
    }
#endif
    FZ_ELEM e_bar_prime[N], v_bar[N];

    csprng_fz_vec(e_bar_prime, &csprng_state);

    fz_vec_sub_n(v_bar, e_bar, e_bar_prime);
    FP_ELEM v[N];
    convert_restr_vec_to_fp(v, v_bar);

    FP_ELEM u[N], u_prime[N];
    csprng_fp_vec(u_prime, &csprng_state);

    fp_vec_by_fp_vec_pointwise(u, v, u_prime);

    FP_ELEM s_prime_sig[N-K] = {0},
            s_prime_ver[N-K] = {0};
#if defined(HIGH_PERFORMANCE_X86_64)
    fp_vec_by_fp_matrix(s_prime_sig, u, V_tr_avx);
#else
    fp_vec_by_fp_matrix(s_prime_sig, u, V_tr);
#endif
    fp_dz_norm_synd(s_prime_sig);

    FP_ELEM y_prime[N], y[N];
    FP_ELEM chall_1 = 1;
    /* y = u_prime + chall_1 * e_bar_prime */
    fp_vec_by_restr_vec_scaled(y, e_bar_prime, chall_1, u_prime);
    fp_dz_norm(y);

    /* Verify starts here */
    fp_vec_by_fp_vec_pointwise(y_prime, v, y);

    FP_ELEM s_prime_tmp_ver[N-K] = {0};
#if defined(HIGH_PERFORMANCE_X86_64)
    fp_vec_by_fp_matrix(s_prime_tmp_ver, y_prime, V_tr_avx);
#else
    fp_vec_by_fp_matrix(s_prime_tmp_ver, y_prime, V_tr);
#endif
    fp_dz_norm_synd(s_prime_tmp_ver);

    fp_synd_minus_fp_vec_scaled(s_prime_ver,
                                s_prime_tmp_ver,
                                chall_1,
                                s);
    fp_dz_norm_synd(s_prime_ver);

    int outcome = (memcmp(s_prime_ver,s_prime_sig,N-K) == 0);
    if(!outcome){
        fprintf(stderr,"syn_tilde_sig: [");
        for(int i=0;i<N-K;i++){
            fprintf(stderr," %d, ", s_prime_sig[i]);
        }
        fprintf(stderr,"]\n");
        fprintf(stderr,"syn_tilde_ver: [");
        for(int i=0;i<N-K;i++){
            fprintf(stderr," %d, ", s_prime_ver[i]);
        }
        fprintf(stderr,"]\n");

    }
    return outcome;
}
#elif defined(RSDPG)

static
void expand_pk(FP_ELEM V_tr[K][N-K],
                        FZ_ELEM W_mat[M][N-M],
                        const uint8_t seed_pub[SEED_LENGTH_BYTES]){
  uint8_t seedV_seed_W[2][SEED_LENGTH_BYTES];
  CSPRNG_STATE_T csprng_state_mat;
  csprng_initialize(&csprng_state_mat, seed_pub, SEED_LENGTH_BYTES,0);
  csprng_randombytes((uint8_t *)seedV_seed_W,
                     2*SEED_LENGTH_BYTES,
                     &csprng_state_mat);

  csprng_initialize(&csprng_state_mat, seedV_seed_W[0], SEED_LENGTH_BYTES,0);
  csprng_fp_mat(V_tr,&csprng_state_mat);

  csprng_initialize(&csprng_state_mat, seedV_seed_W[1], SEED_LENGTH_BYTES,0);
  csprng_fz_mat(W_mat,&csprng_state_mat);
}

static
void expand_sk(FZ_ELEM eta[N],
                         FZ_ELEM zeta[M],
                         FP_ELEM V_tr[K][N-K],
                         FZ_ELEM W_mat[M][N-M],
                         const uint8_t seed[SEED_LENGTH_BYTES]){
  uint8_t seede_seed_pub[2][SEED_LENGTH_BYTES];
  CSPRNG_STATE_T csprng_state;
  csprng_initialize(&csprng_state,seed,SEED_LENGTH_BYTES,0);
  csprng_randombytes((uint8_t *)seede_seed_pub,
                     2*SEED_LENGTH_BYTES,
                     &csprng_state);

  expand_pk(V_tr, W_mat, seede_seed_pub[1]);
  CSPRNG_STATE_T csprng_state_eta;
  csprng_initialize(&csprng_state_eta, seede_seed_pub[0], SEED_LENGTH_BYTES,0);
  csprng_fz_inf_w(zeta, &csprng_state_eta);
#if (defined(HIGH_PERFORMANCE_X86_64) && defined(RSDPG) )
  alignas(EPI8_PER_REG) uint16_t W_mat_avx[M][ROUND_UP(N-M,EPI16_PER_REG)] = {{0}};
  for(int i = 0; i < M; i++){
    for (int j = 0; j < N-M; j++){
       W_mat_avx[i][j] = W_mat[i][j];
    }
  }
  fz_inf_w_by_fz_matrix(eta,zeta,W_mat_avx);
#else
  fz_inf_w_by_fz_matrix(eta,zeta,W_mat);
#endif
  fz_dz_norm_n(eta);
}

int signature_invariant_test(void){
    uint8_t seed[SEED_LENGTH_BYTES];
    randombytes(seed,SEED_LENGTH_BYTES);
    CSPRNG_STATE_T csprng_state;
    csprng_initialize(&csprng_state, seed, SEED_LENGTH_BYTES, 0);

    FZ_ELEM e_bar[N];
    FZ_ELEM e_G_bar[M];
    FP_ELEM V_tr[K][N-K];
    FZ_ELEM W_mat[M][N-M];
    expand_sk(e_bar, e_G_bar, V_tr, W_mat, seed);
#if (defined(HIGH_PERFORMANCE_X86_64) && defined(RSDPG) )
    alignas(EPI8_PER_REG) uint16_t W_mat_avx[M][ROUND_UP(N-M,EPI16_PER_REG)] = {{0}};
    for(int i = 0; i < M; i++){
      for (int j = 0; j < N-M; j++){
         W_mat_avx[i][j] = W_mat[i][j];
      }
    }
#endif
    /* keygen */
    FP_ELEM s[N-K] = {0};
    restr_vec_by_fp_matrix(s, e_bar, V_tr);
    fp_dz_norm_synd(s);

#if (defined(HIGH_PERFORMANCE_X86_64) && defined(RSDPG) )
    alignas(EPI8_PER_REG) FP_DOUBLEPREC V_tr_avx[K][ROUND_UP(N-K,EPI32_PER_REG)] = {{0}};
    for(int i = 0; i < K; i++){
      for (int j = 0; j < N-K; j++){
         V_tr_avx[i][j] = V_tr[i][j];
      }
    }
#endif

    /* sign */
    FZ_ELEM e_G_bar_prime[M], e_bar_prime[N], v_G_bar[M], v_bar[N];

    csprng_fz_inf_w(e_G_bar_prime, &csprng_state);
    fz_vec_sub_m(v_G_bar, e_G_bar, e_G_bar_prime);
#if defined(HIGH_PERFORMANCE_X86_64)
    fz_inf_w_by_fz_matrix(e_bar_prime, e_G_bar_prime, W_mat_avx);
#else
    fz_inf_w_by_fz_matrix(e_bar_prime, e_G_bar_prime, W_mat);
#endif

    fz_vec_sub_n(v_bar, e_bar, e_bar_prime);

    FP_ELEM v[N];
    convert_restr_vec_to_fp(v, v_bar);

    FP_ELEM u[N], u_prime[N];
    csprng_fp_vec(u_prime, &csprng_state);

    fp_vec_by_fp_vec_pointwise(u, v, u_prime);

    FP_ELEM s_prime_sig[N-K] = {0},
            s_prime_ver[N-K] = {0};

#if defined(HIGH_PERFORMANCE_X86_64)
    fp_vec_by_fp_matrix(s_prime_sig, u, V_tr_avx);
#else
    fp_vec_by_fp_matrix(s_prime_sig, u, V_tr);
#endif
    fp_dz_norm_synd(s_prime_sig);

    FP_ELEM y_prime[N], y[N];
    FP_ELEM chall_1 = 1;
    /* y = u_tilde + chall_1 * e_bar_prime */
    fp_vec_by_restr_vec_scaled(y, e_bar_prime, chall_1, u_prime);
    fp_dz_norm(y);

    /* Verify starts here */
    FZ_ELEM v_bar_ver[N];
    FP_ELEM v_ver[N];
#if defined(HIGH_PERFORMANCE_X86_64)
    fz_inf_w_by_fz_matrix(v_bar_ver, v_G_bar, W_mat_avx);
#else
    fz_inf_w_by_fz_matrix(v_bar_ver, v_G_bar, W_mat);
#endif
    convert_restr_vec_to_fp(v_ver, v_bar_ver);

    fp_vec_by_fp_vec_pointwise(y_prime, v_ver, y);

    FP_ELEM s_prime_tmp_ver[N-K] = {0};
#if defined(HIGH_PERFORMANCE_X86_64)
    fp_vec_by_fp_matrix(s_prime_tmp_ver, y_prime, V_tr_avx);
#else
    fp_vec_by_fp_matrix(s_prime_tmp_ver, y_prime, V_tr);
#endif
    fp_dz_norm_synd(s_prime_tmp_ver);

    fp_synd_minus_fp_vec_scaled(s_prime_ver,
                                s_prime_tmp_ver,
                                chall_1,
                                s);
    fp_dz_norm_synd(s_prime_ver);

    int outcome = (memcmp(s_prime_ver, s_prime_sig, N-K) == 0);
    if(!outcome){
        fprintf(stderr,"syn_tilde_sig: [");
        for(int i=0;i<N-K;i++){
            fprintf(stderr," %d, ", s_prime_sig[i]);
        }
        fprintf(stderr,"]\n");
        fprintf(stderr,"syn_tilde_ver: [");
        for(int i=0;i<N-K;i++){
            fprintf(stderr," %d, ", s_prime_ver[i]);
        }
        fprintf(stderr,"]\n");

    }
    return outcome;
}
#endif

int fp_arith_testing(void){
  /*Testing fast constant time reduction against regular one
   * Barrett's method for q = 509
   * k s.t. 2^k > q is k = 9
   * r = floor( 2^(2k)/q)) = 515 */
#if defined(RSDPG)
#define BARRET_R 515
  for (unsigned int i = 0; i < P*P; i++ ){
      uint32_t t = i - ((i*BARRET_R) >> 18)*P;
      int32_t mask = ((int32_t) P - ((int32_t)t+1)) >> 31;
      /*bitmask of all ones iff t >= P*/
      t = ( (t - P) & mask ) | (t & ~mask);
      if(t != (i%P)){
        fprintf(stderr, "(%d,%d,%d,%08X)\n", i, t, i%P, mask);
        return 0;
      }
  }
#endif

  /* test restr_to_val against junior grade school exponentiation */
  for (FZ_ELEM exp = 0; exp < Z; exp++){
    uint32_t check = 1;
    for (int i = exp; i>0; i--){
        check = (check*RESTR_G_GEN) % P;
    }

    FP_ELEM fun;
    fun = RESTR_TO_VAL(exp);
    if( fun != check){
            fprintf(stderr,"mismatch: 16**%u = %u != %u\n",exp,check,fun);
    }
  }
  return 1;
}
