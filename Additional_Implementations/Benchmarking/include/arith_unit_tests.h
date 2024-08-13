/**
 *
 * Reference ISO-C11 Implementation of CROSS.
 *
 * @version 1.1 (March 2023)
 *
 * @author Alessandro Barenghi <alessandro.barenghi@polimi.it>
 * @author Gerardo Pelosi <gerardo.pelosi@polimi.it>
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
    randombytes(seed,SEED_LENGTH_BYTES);
    CSPRNG_STATE_T CSPRNG_state;
    initialize_csprng(&CSPRNG_state,seed,SEED_LENGTH_BYTES);

#if defined(RSDP)
    FZ_ELEM a[N]={0};
    CSPRNG_zz_vec(a,&CSPRNG_state);
    return trivial_is_restricted(a,N) == is_zz_vec_in_restr_group(a);
#else
    FZ_ELEM a[M]={0};
    CSPRNG_zz_inf_w(a,&CSPRNG_state);
    return trivial_is_restricted(a,M) == is_zz_inf_w_valid(a);
#endif
}

void prettyprint(char* name, FQ_ELEM *v, int l){
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
    CSPRNG_STATE_T CSPRNG_state;
    initialize_csprng(&CSPRNG_state,seed,SEED_LENGTH_BYTES);

    /* keygen */
    FZ_ELEM e[N];
    FQ_ELEM V_tr[K][N-K];
    CSPRNG_zz_vec(e, &CSPRNG_state);
    CSPRNG_fq_mat(V_tr,&CSPRNG_state);
    FQ_ELEM public_key_synd[N-K] = {0};
    restr_vec_by_fq_matrix(public_key_synd,e,V_tr);
    fq_dz_norm_synd(public_key_synd);

#if defined(HIGH_PERFORMANCE_X86_64)
    /* Expanded */
    alignas(EPI8_PER_REG) FQ_DOUBLEPREC V_tr_avx[K][ROUND_UP(N-K,EPI16_PER_REG)] = {{0}};
    for(int i = 0; i < K; i++){
      for (int j = 0; j < N-K; j++){
         V_tr_avx[i][j] = V_tr[i][j];
      }
    }
#endif
    FZ_ELEM e_tilde[N],sigma[N];

    CSPRNG_zz_vec(e_tilde, &CSPRNG_state);

    restr_vec_sub(sigma,e,e_tilde);
    FQ_ELEM tau[N];
    convert_restr_vec_to_fq(tau,sigma);

    FQ_ELEM u[N],u_tilde[N];
    CSPRNG_fq_vec(u_tilde, &CSPRNG_state);

    fq_vec_by_fq_vec_pointwise(u,tau,u_tilde);

    FQ_ELEM syn_tilde_sig[N-K] = {0},
            syn_tilde_ver[N-K] = {0};
#if defined(HIGH_PERFORMANCE_X86_64)
    fq_vec_by_fq_matrix(syn_tilde_sig,u,V_tr_avx);
#else
    fq_vec_by_fq_matrix(syn_tilde_sig,u,V_tr);
#endif
    fq_dz_norm_synd(syn_tilde_sig);

    FQ_ELEM y_tilde[N], y_sig[N];
    FQ_ELEM beta = 1;
    /* y = u_tilde + beta * e*/
    fq_vec_by_restr_vec_scaled(y_sig, e_tilde, beta, u_tilde);
    fq_dz_norm(y_sig);

    /* Verify starts here */
    fq_vec_by_fq_vec_pointwise(y_tilde,tau,y_sig);

    FQ_ELEM syn_tilde_tmp_ver[N-K] = {0};
#if defined(HIGH_PERFORMANCE_X86_64)
    fq_vec_by_fq_matrix(syn_tilde_tmp_ver,y_tilde,V_tr_avx);
#else
    fq_vec_by_fq_matrix(syn_tilde_tmp_ver,y_tilde,V_tr);
#endif
    fq_dz_norm_synd(syn_tilde_tmp_ver);

    fq_synd_minus_fq_vec_scaled(syn_tilde_ver,
                                syn_tilde_tmp_ver,
                                beta,
                                public_key_synd);
    fq_dz_norm_synd(syn_tilde_ver);

    int outcome = (memcmp(syn_tilde_ver,syn_tilde_sig,N-K) == 0);
    if(!outcome){
        fprintf(stderr,"syn_tilde_sig: [");
        for(int i=0;i<N-K;i++){
            fprintf(stderr," %d, ", syn_tilde_sig[i]);
        }
        fprintf(stderr,"]\n");
        fprintf(stderr,"syn_tilde_ver: [");
        for(int i=0;i<N-K;i++){
            fprintf(stderr," %d, ", syn_tilde_ver[i]);
        }
        fprintf(stderr,"]\n");

    }
    return outcome;
}
#elif defined(RSDPG)

static
void expand_public_seed(FQ_ELEM V_tr[K][N-K],
                        FZ_ELEM W_mat[M][N-M],
                        const uint8_t seed_pub[SEED_LENGTH_BYTES]){
  uint8_t seedV_seed_W[2][SEED_LENGTH_BYTES];
  CSPRNG_STATE_T CSPRNG_state_mat;
  initialize_csprng(&CSPRNG_state_mat, seed_pub, SEED_LENGTH_BYTES);
  csprng_randombytes((uint8_t *)seedV_seed_W,
                     2*SEED_LENGTH_BYTES,
                     &CSPRNG_state_mat);

  initialize_csprng(&CSPRNG_state_mat, seedV_seed_W[0], SEED_LENGTH_BYTES);
  CSPRNG_fq_mat(V_tr,&CSPRNG_state_mat);

  initialize_csprng(&CSPRNG_state_mat, seedV_seed_W[1], SEED_LENGTH_BYTES);
  CSPRNG_fz_mat(W_mat,&CSPRNG_state_mat);
}

static
void expand_private_seed(FZ_ELEM eta[N],
                         FZ_ELEM zeta[M],
                         FQ_ELEM V_tr[K][N-K],
                         FZ_ELEM W_mat[M][N-M],
                         const uint8_t seed[SEED_LENGTH_BYTES]){
  uint8_t seede_seed_pub[2][SEED_LENGTH_BYTES];
  CSPRNG_STATE_T CSPRNG_state;
  initialize_csprng(&CSPRNG_state,seed,SEED_LENGTH_BYTES);
  csprng_randombytes((uint8_t *)seede_seed_pub,
                     2*SEED_LENGTH_BYTES,
                     &CSPRNG_state);

  expand_public_seed(V_tr,W_mat,seede_seed_pub[1]);
  CSPRNG_STATE_T CSPRNG_state_eta;
  initialize_csprng(&CSPRNG_state_eta, seede_seed_pub[0], SEED_LENGTH_BYTES);
  CSPRNG_zz_inf_w(zeta,&CSPRNG_state_eta);
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
  fz_dz_norm_sigma(eta);
}

int signature_invariant_test(void){
    uint8_t seed[SEED_LENGTH_BYTES];
    randombytes(seed,SEED_LENGTH_BYTES);
    CSPRNG_STATE_T CSPRNG_state;
    initialize_csprng(&CSPRNG_state,seed,SEED_LENGTH_BYTES);

    FZ_ELEM eta[N];
    FZ_ELEM zeta[M];
    FQ_ELEM V_tr[K][N-K];
    FZ_ELEM W_mat[M][N-M];
    expand_private_seed(eta, zeta, V_tr, W_mat,seed);
#if (defined(HIGH_PERFORMANCE_X86_64) && defined(RSDPG) )
    alignas(EPI8_PER_REG) uint16_t W_mat_avx[M][ROUND_UP(N-M,EPI16_PER_REG)] = {{0}};
    for(int i = 0; i < M; i++){
      for (int j = 0; j < N-M; j++){
         W_mat_avx[i][j] = W_mat[i][j];
      }
    }
#endif
    /* keygen */
    FQ_ELEM public_key_synd[N-K] = {0};
    restr_vec_by_fq_matrix(public_key_synd,eta,V_tr);
    fq_dz_norm_synd(public_key_synd);

    /* sign */
    FZ_ELEM zeta_tilde[M],eta_tilde[N],delta[M],sigma[N];

    // CSPRNG_zz_vec(e_tilde, &CSPRNG_state);
    CSPRNG_zz_inf_w(zeta_tilde, &CSPRNG_state);
    restr_inf_w_sub(delta, zeta,zeta_tilde);
#if defined(HIGH_PERFORMANCE_X86_64)
    fz_inf_w_by_fz_matrix(eta_tilde,zeta_tilde,W_mat_avx);
#else
    fz_inf_w_by_fz_matrix(eta_tilde,zeta_tilde,W_mat);
#endif

    restr_vec_sub(sigma,eta,eta_tilde);

    FQ_ELEM tau[N];
    convert_restr_vec_to_fq(tau,sigma);

    FQ_ELEM u[N],u_tilde[N];
    CSPRNG_fq_vec(u_tilde, &CSPRNG_state);

    fq_vec_by_fq_vec_pointwise(u,tau,u_tilde);

    FQ_ELEM syn_tilde_sig[N-K] = {0},
            syn_tilde_ver[N-K] = {0};

    fq_vec_by_fq_matrix(syn_tilde_sig,u,V_tr);
    fq_dz_norm_synd(syn_tilde_sig);

    FQ_ELEM y_tilde[N], y_sig[N];
    FQ_ELEM beta = 1;
    /* y = u_tilde + beta * e*/
    fq_vec_by_restr_vec_scaled(y_sig, eta_tilde, beta, u_tilde);
    fq_dz_norm(y_sig);

    /* Verify starts here */
    FZ_ELEM sigma_ver[N];
    FQ_ELEM tau_ver[N];
#if defined(HIGH_PERFORMANCE_X86_64)
    fz_inf_w_by_fz_matrix(sigma_ver,delta,W_mat_avx);
#else
    fz_inf_w_by_fz_matrix(sigma_ver,delta,W_mat);
#endif
    convert_restr_vec_to_fq(tau_ver,sigma);

    fq_vec_by_fq_vec_pointwise(y_tilde,tau_ver,y_sig);

    FQ_ELEM syn_tilde_tmp_ver[N-K] = {0};
    fq_vec_by_fq_matrix(syn_tilde_tmp_ver,y_tilde,V_tr);
    fq_dz_norm_synd(syn_tilde_tmp_ver);

    fq_synd_minus_fq_vec_scaled(syn_tilde_ver,
                                syn_tilde_tmp_ver,
                                beta,
                                public_key_synd);
    fq_dz_norm_synd(syn_tilde_ver);

    int outcome = (memcmp(syn_tilde_ver,syn_tilde_sig,N-K) == 0);
    if(!outcome){
        fprintf(stderr,"syn_tilde_sig: [");
        for(int i=0;i<N-K;i++){
            fprintf(stderr," %d, ", syn_tilde_sig[i]);
        }
        fprintf(stderr,"]\n");
        fprintf(stderr,"syn_tilde_ver: [");
        for(int i=0;i<N-K;i++){
            fprintf(stderr," %d, ", syn_tilde_ver[i]);
        }
        fprintf(stderr,"]\n");

    }
    return outcome;
}
#endif

int fq_arith_testing(void){
  /*Testing fast constant time reduction against regular one
   * Barrett's method for q = 509
   * k s.t. 2^k > q is k = 9
   * r = floor( 2^(2k)/q)) = 515 */
#if defined(RSDPG)
#define BARRET_R 515
  for (unsigned int i = 0; i < Q*Q; i++ ){
      uint32_t t = i - ((i*BARRET_R) >> 18)*Q;
      int32_t mask = ((int32_t) Q - ((int32_t)t+1)) >> 31;
      /*bitmask of all ones iff t >= Q*/
      t = ( (t - Q) & mask ) | (t & ~mask);
      if(t != (i%Q)){
        fprintf(stderr,"(%d,%d,%d,%08X)\n",i,t,i%Q,mask);
        return 0;
      }
  }
#endif

  /* test restr_to_val against junior grade school exponentiation */
  for (FZ_ELEM exp = 0; exp < Z; exp++){
    uint32_t check = 1;
    for (int i = exp; i>0; i--){
        check = (check*RESTR_G_GEN) % Q;
    }

    FQ_ELEM fun;
    fun = RESTR_TO_VAL(exp);
    if( fun != check){
            fprintf(stderr,"mismatch: 16**%u = %u != %u\n",exp,check,fun);
    }
  }
  return 1;
}
