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

#include "parameters.h"
#include "sha3.h"

/************************* CSPRNG ********************************/

#define CSPRNG_STATE_T SHAKE_STATE_STRUCT
/* initializes a CSPRNG, given the seed and a state pointer */
static inline
void initialize_csprng(CSPRNG_STATE_T * const csprng_state,
                       const unsigned char * const seed,
                       const uint32_t seed_len_bytes) {
   // the second parameter is the security level of the SHAKE instance
   xof_shake_init(csprng_state, SEED_LENGTH_BYTES*8);
   xof_shake_update(csprng_state,seed,seed_len_bytes);
   xof_shake_final(csprng_state);
} /* end initialize_csprng */

/* extracts xlen bytes from the CSPRNG, given the state */
static inline
void csprng_randombytes(unsigned char * const x,
                        unsigned long long xlen,
                        CSPRNG_STATE_T * const csprng_state){
   xof_shake_extract(csprng_state,x,xlen);
}

/******************************************************************************/

/* global csprng state employed to have deterministic randombytes for testing */
extern CSPRNG_STATE_T platform_csprng_state;
/* extracts xlen bytes from the global CSPRNG */
static inline
void randombytes(unsigned char * x,
                 unsigned long long xlen) {
   csprng_randombytes(x,xlen,&platform_csprng_state);
}

/************************* HASH functions ********************************/

/* Opaque algorithm agnostic hash call */
static inline
void hash(uint8_t digest[HASH_DIGEST_LENGTH],
          const unsigned char *const m,
          const uint64_t mlen){
   /* SHAKE with a 2*lambda bit digest is employed also for hashing */
   CSPRNG_STATE_T csprng_state;    
   xof_shake_init(&csprng_state, SEED_LENGTH_BYTES*8);
   xof_shake_update(&csprng_state,m,mlen);
   xof_shake_final(&csprng_state);    
   xof_shake_extract(&csprng_state,digest,HASH_DIGEST_LENGTH);
}


/********************** CSPRNG Sampling functions helpers ********************/

static inline
FQ_ELEM fq_star_rnd_state(CSPRNG_STATE_T * const csprng_state)
{
   const FQ_ELEM mask = ( (FQ_ELEM) 1 << BITS_TO_REPRESENT(Q-2) ) - 1;
   FQ_ELEM rnd_value;
   do {
      csprng_randombytes((unsigned char *) &rnd_value,
                         sizeof(FQ_ELEM),
                         csprng_state);
      rnd_value = mask & rnd_value;
   } while (rnd_value > Q-2);

   return rnd_value+1;
} /* end fq_star_rnd_state */


/***************** Specialized CSPRNGs for non binary domains *****************/

/* CSPRNG sampling fixed weight strings */
void expand_digest_to_fixed_weight(uint8_t fixed_weight_string[T],
                                   const uint8_t digest[HASH_DIGEST_LENGTH]);

#define BITS_FOR_Q BITS_TO_REPRESENT(Q-1) 
#define BITS_FOR_Z BITS_TO_REPRESENT(Z-1) 

static inline
void CSPRNG_fq_vec(FQ_ELEM res[N],
                   CSPRNG_STATE_T * const csprng_state){
    const FQ_ELEM mask = ( (FQ_ELEM) 1 << BITS_FOR_Q) - 1;
    uint8_t CSPRNG_buffer[ROUND_UP(BITS_N_ZQ_CT_RNG,8)/8];
    /* To facilitate hardware implementations, the uint64_t 
     * sub-buffer is consumed starting from the least significant byte 
     * i.e., from the first being output by SHAKE. Bits in the byte are 
     * discarded shifting them out to the right, shifting fresh ones
     * in from the left end */
    csprng_randombytes(CSPRNG_buffer,sizeof(CSPRNG_buffer),csprng_state);    
    int placed = 0;
    uint64_t sub_buffer = *(uint64_t*)CSPRNG_buffer;
    int bits_in_sub_buf = 64;
    /* position of the next fresh byte in CSPRNG_buffer*/
    int pos_in_buf = 8;
    while(placed < N){
        if(bits_in_sub_buf <= 32){
            /* get 32 fresh bits from main buffer with a single load */
            uint32_t refresh_buf = *(uint32_t*) (CSPRNG_buffer+pos_in_buf);
            pos_in_buf += 4;
            sub_buffer |=  ((uint64_t) refresh_buf) << bits_in_sub_buf;
            bits_in_sub_buf += 32; 
        }
        res[placed] = sub_buffer & mask;
        if(res[placed] < Q) {
           placed++;
           sub_buffer = sub_buffer >> BITS_FOR_Q;
           bits_in_sub_buf -= BITS_FOR_Q;
            
        } else {
           sub_buffer = sub_buffer >> 1;
           bits_in_sub_buf -= 1;
        }
    }
}

#define BITS_FOR_Q_M_ONE BITS_TO_REPRESENT(Q-2) 

static inline
void CSPRNG_fq_vec_beta(FQ_ELEM res[T],
                   CSPRNG_STATE_T * const csprng_state){
    const FQ_ELEM mask = ( (FQ_ELEM) 1 << BITS_FOR_Q_M_ONE) - 1;
    uint8_t CSPRNG_buffer[ROUND_UP(BITS_BETA_ZQSTAR_CT_RNG,8)/8];
    /* To facilitate hardware implementations, the uint64_t 
     * sub-buffer is consumed starting from the least significant byte 
     * i.e., from the first being output by SHAKE. Bits in the byte are 
     * discarded shifting them out to the right , shifting fresh ones
     * in from the left end */
    csprng_randombytes(CSPRNG_buffer,sizeof(CSPRNG_buffer),csprng_state);    
    int placed = 0;
    uint64_t sub_buffer = *(uint64_t*)CSPRNG_buffer;
    int bits_in_sub_buf = 64;
    /* position of the next fresh byte in CSPRNG_buffer*/
    int pos_in_buf = 8;
    while(placed < T){
        if(bits_in_sub_buf <= 32){
            /* get 32 fresh bits from main buffer with a single load */
            uint32_t refresh_buf = *(uint32_t*) (CSPRNG_buffer+pos_in_buf);
            pos_in_buf += 4;
            sub_buffer |=  ((uint64_t) refresh_buf) << bits_in_sub_buf;
            bits_in_sub_buf += 32; 
        }
        /* draw from 0 ... Q-2, then add 1*/
        res[placed] = (sub_buffer & mask)+1;
        if(res[placed] < Q) {
           placed++;
           sub_buffer = sub_buffer >> BITS_FOR_Q_M_ONE;
           bits_in_sub_buf -= BITS_FOR_Q_M_ONE;
        } else {
           sub_buffer = sub_buffer >> 1;
           bits_in_sub_buf -= 1;
        }
    }
}

static inline
void CSPRNG_fq_mat(FQ_ELEM res[K][N-K],
                   CSPRNG_STATE_T * const csprng_state){
    const FQ_ELEM mask = ( (FQ_ELEM) 1 << BITS_TO_REPRESENT(Q-1)) - 1;
    uint8_t CSPRNG_buffer[ROUND_UP(BITS_V_CT_RNG,8)/8];
    /* To facilitate hardware implementations, the uint64_t 
     * sub-buffer is consumed starting from the least significant byte 
     * i.e., from the first being output by SHAKE. Bits in the byte are 
     * discarded shifting them out to the right , shifting fresh ones
     * in from the left end */
    csprng_randombytes(CSPRNG_buffer,sizeof(CSPRNG_buffer),csprng_state);    
    int placed = 0;
    uint64_t sub_buffer = *(uint64_t*)CSPRNG_buffer;
    int bits_in_sub_buf = 64;
    /* position of the next fresh byte in CSPRNG_buffer*/
    int pos_in_buf = 8;
    while(placed < K*(N-K)){
        if(bits_in_sub_buf <= 32){
            /* get 32 fresh bits from main buffer with a single load */
            uint32_t refresh_buf = *(uint32_t*) (CSPRNG_buffer+pos_in_buf);
            pos_in_buf += 4;
            sub_buffer |=  ((uint64_t) refresh_buf) << bits_in_sub_buf;
            bits_in_sub_buf += 32; 
        }
        *((FQ_ELEM*)res+placed) = sub_buffer & mask;
        if(*((FQ_ELEM*)res+placed) < Q) {
           placed++;
           sub_buffer = sub_buffer >> BITS_FOR_Q;
           bits_in_sub_buf -= BITS_FOR_Q;
            
        } else {
           sub_buffer = sub_buffer >> 1;
           bits_in_sub_buf -= 1;
        }
    }   
}

#if defined(RSDP)
static inline
void CSPRNG_zz_vec(FZ_ELEM res[N],
                   CSPRNG_STATE_T * const csprng_state){
    const FZ_ELEM mask = ( (FZ_ELEM) 1 << BITS_TO_REPRESENT(Z-1)) - 1;
    uint8_t CSPRNG_buffer[ROUND_UP(BITS_N_ZZ_CT_RNG,8)/8];
    /* To facilitate hardware implementations, the uint64_t 
     * sub-buffer is consumed starting from the least significant byte 
     * i.e., from the first being output by SHAKE. Bits in the byte are 
     * discarded shifting them out to the right , shifting fresh ones
     * in from the left end */
    csprng_randombytes(CSPRNG_buffer,sizeof(CSPRNG_buffer),csprng_state);    
    int placed = 0;
    uint64_t sub_buffer = *(uint64_t*)CSPRNG_buffer;
    int bits_in_sub_buf = 64;
    /* position of the next fresh byte in CSPRNG_buffer*/
    int pos_in_buf = 8;
    while(placed < N){
        if(bits_in_sub_buf <= 32){
            /* get 32 fresh bits from main buffer with a single load */
            uint32_t refresh_buf = *(uint32_t*) (CSPRNG_buffer+pos_in_buf);
            pos_in_buf += 4;
            sub_buffer |=  ((uint64_t) refresh_buf) << bits_in_sub_buf;
            bits_in_sub_buf += 32; 
      }
        /* get */
        res[placed] = sub_buffer & mask;
        if(res[placed] < Z) {
           placed++;
           sub_buffer = sub_buffer >> BITS_FOR_Z;
           bits_in_sub_buf -= BITS_FOR_Z;
            
        } else {
           sub_buffer = sub_buffer >> 1;
           bits_in_sub_buf -= 1;
        }
    }
}
#elif defined(RSDPG)
static inline
void CSPRNG_zz_inf_w(FZ_ELEM res[M],
                   CSPRNG_STATE_T * const csprng_state){
    const FZ_ELEM mask = ( (FZ_ELEM) 1 << BITS_TO_REPRESENT(Z-1)) - 1;
    uint8_t CSPRNG_buffer[ROUND_UP(BITS_M_ZZ_CT_RNG,8)/8];
    /* To facilitate hardware implementations, the uint64_t 
     * sub-buffer is consumed starting from the least significant byte 
     * i.e., from the first being output by SHAKE. Bits in the byte are 
     * discarded shifting them out to the right , shifting fresh ones
     * in from the left end */
    csprng_randombytes(CSPRNG_buffer,sizeof(CSPRNG_buffer),csprng_state);    
    int placed = 0;
    uint64_t sub_buffer = *(uint64_t*)CSPRNG_buffer;
    int bits_in_sub_buf = 64;
    /* position of the next fresh byte in CSPRNG_buffer*/
    int pos_in_buf = 8;
    while(placed < M){
        if(bits_in_sub_buf <= 32){
            /* get 32 fresh bits from main buffer with a single load */
            uint32_t refresh_buf = *(uint32_t*) (CSPRNG_buffer+pos_in_buf);
            pos_in_buf += 4;
            sub_buffer |=  ((uint64_t) refresh_buf) << bits_in_sub_buf;
            bits_in_sub_buf += 32;             
        }
        res[placed] = sub_buffer & mask;
        if(res[placed] < Z) {
           placed++;
           sub_buffer = sub_buffer >> BITS_FOR_Z;
           bits_in_sub_buf -= BITS_FOR_Z;
            
        } else {
           sub_buffer = sub_buffer >> 1;
           bits_in_sub_buf -= 1;
        }
    }
}

static inline
void CSPRNG_fz_mat(FZ_ELEM res[M][N-M],
                   CSPRNG_STATE_T * const csprng_state){
    const FZ_ELEM mask = ( (FZ_ELEM) 1 << BITS_TO_REPRESENT(Z-1)) - 1;
    uint8_t CSPRNG_buffer[ROUND_UP(BITS_W_CT_RNG,8)/8];
    /* To facilitate hardware implementations, the uint64_t 
     * sub-buffer is consumed starting from the least significant byte 
     * i.e., from the first being output by SHAKE. Bits in the byte are 
     * discarded shifting them out to the right , shifting fresh ones
     * in from the left end */
    csprng_randombytes(CSPRNG_buffer,sizeof(CSPRNG_buffer),csprng_state);    
    int placed = 0;
    uint64_t sub_buffer = *(uint64_t*)CSPRNG_buffer;
    int bits_in_sub_buf = 64;
    /* position of the next fresh byte in CSPRNG_buffer*/
    int pos_in_buf = 8;
    while(placed < M*(N-M)){
        if(bits_in_sub_buf <= 32){
            /* get 32 fresh bits from main buffer with a single load */
            uint32_t refresh_buf = *(uint32_t*) (CSPRNG_buffer+pos_in_buf);
            pos_in_buf += 4;
            sub_buffer |=  ((uint64_t) refresh_buf) << bits_in_sub_buf;
            bits_in_sub_buf += 32;             
        }
        *((FZ_ELEM*)res+placed) = sub_buffer & mask;
        if(*((FZ_ELEM*)res+placed) < Z) {
           placed++;
           sub_buffer = sub_buffer >> BITS_FOR_Z;
           bits_in_sub_buf -= BITS_FOR_Z;
            
        } else {
           sub_buffer = sub_buffer >> 1;
           bits_in_sub_buf -= 1;
        }
    }    
}
#endif
