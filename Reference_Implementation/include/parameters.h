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
#include <stdint.h>

/******************************************************************************/
/*************************** Base Fields Parameters ***************************/
/******************************************************************************/
#if defined(RSDP)

/* The same base field and restriction are employed for all categories of RSDP */
#define   Q (127)
#define   Z (  7)
/* single-register table representation of E, the value of g^7=1 is also
 * represented to avoid exponent renormalization*/
#define RESTR_G_TABLE ((uint64_t) (0x0140201008040201))
#define RESTR_G_GEN 2
#define FQ_ELEM uint8_t
#define FZ_ELEM uint8_t
#define FQ_DOUBLEPREC uint16_t
#define FQ_TRIPLEPREC uint32_t
#elif defined(RSDPG)

/* The same base field and restriction are employed for all categories of RSDP */
#define   Q (509)
#define   Z (127)
/* Restricted subgroup generator */
#define RESTR_G_GEN 16
#define FZ_ELEM uint8_t
#define FZ_DOUBLEPREC uint16_t
#define FQ_ELEM uint16_t
#define FQ_DOUBLEPREC uint32_t
#define FQ_TRIPLEPREC uint32_t

#else
#error define either RSDP or RSDPG
#endif

/******************************************************************************/
/****************************** RSDP Parameters *******************************/
/******************************************************************************/
#if defined(RSDP)
/********************************* Category 1 *********************************/
#if defined(CATEGORY_1)
#define SEC_MARGIN_LAMBDA (128)
#define   N (127)
#define   K ( 76)

#if defined(SPEED)
#define   T (163)
#define   W ( 85)
#define POSITION_IN_FW_STRING_T uint16_t
#elif defined(BALANCED)
#define   T (252)
#define   W (212)
#define POSITION_IN_FW_STRING_T uint16_t
#elif defined(SIG_SIZE)
#define   T (960)
#define   W (938)
#define POSITION_IN_FW_STRING_T uint16_t

#else
#error define optimization corner in Cmakelist
#endif

/********************************* Category 3 *********************************/
#elif defined(CATEGORY_3)
#define SEC_MARGIN_LAMBDA (192)
#define   N (187)
#define   K (111)

#if defined(SPEED)
#define   T (245)
#define   W (127)
#define POSITION_IN_FW_STRING_T uint16_t
#elif defined(BALANCED)
#define   T (398)
#define   W (340)
#define POSITION_IN_FW_STRING_T uint16_t
#elif defined(SIG_SIZE)
#define   T (945)
#define   W (907)
#define POSITION_IN_FW_STRING_T uint16_t

#else
#error define optimization corner in Cmakelist
#endif

/********************************* Category 5 *********************************/
#elif defined(CATEGORY_5)
#define SEC_MARGIN_LAMBDA (256)
#define   N (251)
#define   K (150)

#if defined(SPEED)
#define   T (327)
#define   W (169)
#define POSITION_IN_FW_STRING_T uint16_t
#elif defined(BALANCED)
#define   T (507)
#define   W (427)
#define POSITION_IN_FW_STRING_T uint16_t
#elif defined(SIG_SIZE)
#define   T (968)
#define   W (912)
#define POSITION_IN_FW_STRING_T uint16_t

#else
#error define optimization corner in Cmakelist
#endif

#else
#error define category for parameters
#endif

/******************************************************************************/
/****************************** RSDP(G) Parameters ****************************/
/******************************************************************************/
#elif defined(RSDPG)
/********************************* Category 1 *********************************/
#if defined(CATEGORY_1)
#define SEC_MARGIN_LAMBDA (128)
#define   N ( 55)
#define   K ( 36)
#define   M ( 25)

#if defined(SPEED)
#define   T (153)
#define   W (79)
#define POSITION_IN_FW_STRING_T uint8_t
#elif defined(BALANCED)
#define   T (243)
#define   W (206)
#define POSITION_IN_FW_STRING_T uint8_t
#elif defined(SIG_SIZE)
#define   T (871)
#define   W (850)
#define POSITION_IN_FW_STRING_T uint16_t

#else
#error define optimization corner in Cmakelist
#endif

/********************************* Category 3 *********************************/
#elif defined(CATEGORY_3)
#define SEC_MARGIN_LAMBDA (192)
#define   N ( 79)
#define   K ( 48)
#define   M ( 40)

#if defined(SPEED)
#define   T (230)
#define   W (123)
#define POSITION_IN_FW_STRING_T uint8_t
#elif defined(BALANCED)
#define   T (255)
#define   W (176)
#define POSITION_IN_FW_STRING_T uint8_t
#elif defined(SIG_SIZE)
#define   T (949)
#define   W (914)
#define POSITION_IN_FW_STRING_T uint16_t

#else
#error define optimization corner in Cmakelist
#endif

/********************************* Category 5 *********************************/
#elif defined(CATEGORY_5)
#define SEC_MARGIN_LAMBDA (256)
#define   N (106)
#define   K ( 69)
#define   M ( 48)

#if defined(SPEED)
#define   T (306)
#define   W (157)
#define POSITION_IN_FW_STRING_T uint16_t
#elif defined(BALANCED)
#define   T (356)
#define   W (257)
#define POSITION_IN_FW_STRING_T uint16_t
#elif defined(SIG_SIZE)
#define   T (996)
#define   W (945)
#define POSITION_IN_FW_STRING_T uint16_t

#else
#error define optimization corner in Cmakelist
#endif

#else
#error define category for parameters
#endif

#else
#error define either RSDP or RSDPG
#endif


#define HASH_CSPRNG_DOMAIN_SEP_CONST ((uint16_t)32768)

/************* Helper macros for derived parameter computation ***************/

#define ROUND_UP(amount, round_amt) ( ((amount+round_amt-1)/round_amt)*round_amt )

#define IS_REPRESENTABLE_IN_D_BITS(D, N)                \
  (((unsigned long) N>=(1UL << (D-1)) && (unsigned long) N<(1UL << D)) ? D : -1)

#define BITS_TO_REPRESENT(N)                            \
  (N == 0 ? 1 : (15                                     \
                 + IS_REPRESENTABLE_IN_D_BITS( 1, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 2, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 3, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 4, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 5, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 6, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 7, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 8, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 9, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(10, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(11, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(12, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(13, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(14, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(15, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(16, N)    \
                 )                                      \
   )

#define LOG2(L) ( (BITS_TO_REPRESENT(L) > BITS_TO_REPRESENT(L-1)) ? (BITS_TO_REPRESENT(L-1)) : (BITS_TO_REPRESENT(L)) )

/***************** Derived parameters *****************************************/
#define SEED_LENGTH_BYTES (SEC_MARGIN_LAMBDA/8)
#define KEYPAIR_SEED_LENGTH_BYTES (2*(SEC_MARGIN_LAMBDA/8))
#define HASH_DIGEST_LENGTH (2*(SEC_MARGIN_LAMBDA/8))
#define SALT_LENGTH_BYTES (2*(SEC_MARGIN_LAMBDA/8))

#define NUM_LEAVES_MERKLE_TREE (T)
#define NUM_NODES_MERKLE_TREE (2*NUM_LEAVES_MERKLE_TREE-1)

/*to be derived via script for each T/W*/
#define NUM_LEAVES_SEED_TREE ( T )
// #define NUM_NODES_SEED_TREE ( 2*NUM_LEAVES_SEED_TREE-1 )
#define NUM_INNER_NODES_SEED_TREE ( NUM_NODES_SEED_TREE-NUM_LEAVES_SEED_TREE )

/* Sizes of bitpacked field element vectors
 * Bitpacking an n-elements vector of num_bits_for_q-1 bits long values
 * will pack 8 values in num_bits_for_q-1 bytes exactly, leaving the remaining
 * N % 8 as a tail */
#define DENSELY_PACKED_FQ_VEC_SIZE ((N/8)*BITS_TO_REPRESENT(Q-1) + \
                                   ROUND_UP( ((N%8)*BITS_TO_REPRESENT(Q-1)),8)/8)
#define DENSELY_PACKED_FQ_SYN_SIZE (((N-K)/8)*BITS_TO_REPRESENT(Q-1) + \
                                   ROUND_UP( (((N-K)%8)*BITS_TO_REPRESENT(Q-1)),8)/8)
#define DENSELY_PACKED_FZ_VEC_SIZE ((N/8)*BITS_TO_REPRESENT(Z-1) + \
                                   ROUND_UP( ((N%8)*BITS_TO_REPRESENT(Z-1)),8)/8)
#ifdef RSDPG
#define DENSELY_PACKED_FZ_RSDP_G_VEC_SIZE ((M/8)*BITS_TO_REPRESENT(Z-1) + \
                                          ROUND_UP( ((M%8)*BITS_TO_REPRESENT(Z-1)),8)/8)
#endif


/* Derived parameters computed via compute_derived_parameters.py */
#if ( defined(CATEGORY_1) && defined(RSDP)  && defined(SPEED) )
#define TREE_NODES_TO_STORE 83
#define NUM_NODES_SEED_TREE 330
#define NODES_PER_LEVEL_ARRAY {1, 2, 3, 6, 11, 21, 41, 82, 163}
#define MISSING_NODES_BEFORE_LEVEL_ARRAY {0, 0, 0, 1, 3, 8, 19, 42, 88}
#define BITS_N_ZQ_CT_RNG 923
#define BITS_BETA_ZQSTAR_CT_RNG 1187
#define BITS_V_CT_RNG 27260
#define BITS_N_ZZ_CT_RNG 493
#define BITS_CWSTR_RNG 1355

#elif ( defined(CATEGORY_1) && defined(RSDP)  && defined(BALANCED) )
#define TREE_NODES_TO_STORE 107
#define NUM_NODES_SEED_TREE 504
#define NODES_PER_LEVEL_ARRAY {1, 2, 4, 8, 16, 32, 63, 126, 252}
#define MISSING_NODES_BEFORE_LEVEL_ARRAY {0, 0, 0, 0, 0, 0, 0, 1, 3}
#define BITS_N_ZQ_CT_RNG 923
#define BITS_BETA_ZQSTAR_CT_RNG 1817
#define BITS_V_CT_RNG 27260
#define BITS_N_ZZ_CT_RNG 493
#define BITS_CWSTR_RNG 2102

#elif ( defined(CATEGORY_1) && defined(RSDP)  && defined(SIG_SIZE) )
#define TREE_NODES_TO_STORE 120
#define NUM_NODES_SEED_TREE 1920
#define NODES_PER_LEVEL_ARRAY {1, 2, 4, 8, 15, 30, 60, 120, 240, 480, 960}
#define MISSING_NODES_BEFORE_LEVEL_ARRAY {0, 0, 0, 0, 0, 1, 3, 7, 15, 31, 63}
#define BITS_N_ZQ_CT_RNG 923
#define BITS_BETA_ZQSTAR_CT_RNG 6811
#define BITS_V_CT_RNG 27260
#define BITS_N_ZZ_CT_RNG 493
#define BITS_CWSTR_RNG 9371

#elif ( defined(CATEGORY_3) && defined(RSDP)  && defined(SPEED) )
#define TREE_NODES_TO_STORE 125
#define NUM_NODES_SEED_TREE 492
#define NODES_PER_LEVEL_ARRAY {1, 2, 4, 8, 16, 31, 62, 123, 245}
#define MISSING_NODES_BEFORE_LEVEL_ARRAY {0, 0, 0, 0, 0, 0, 1, 3, 8}
#define BITS_N_ZQ_CT_RNG 1361
#define BITS_BETA_ZQSTAR_CT_RNG 1785
#define BITS_V_CT_RNG 59289
#define BITS_N_ZZ_CT_RNG 729
#define BITS_CWSTR_RNG 2125

#elif ( defined(CATEGORY_3) && defined(RSDP)  && defined(BALANCED) )
#define TREE_NODES_TO_STORE 162
#define NUM_NODES_SEED_TREE 799
#define NODES_PER_LEVEL_ARRAY {1, 2, 4, 7, 13, 25, 50, 100, 199, 398}
#define MISSING_NODES_BEFORE_LEVEL_ARRAY {0, 0, 0, 0, 1, 4, 11, 25, 53, 110}
#define BITS_N_ZQ_CT_RNG 1361
#define BITS_BETA_ZQSTAR_CT_RNG 2868
#define BITS_V_CT_RNG 59289
#define BITS_N_ZZ_CT_RNG 729
#define BITS_CWSTR_RNG 3647

#elif ( defined(CATEGORY_3) && defined(RSDP)  && defined(SIG_SIZE) )
#define TREE_NODES_TO_STORE 177
#define NUM_NODES_SEED_TREE 1894
#define NODES_PER_LEVEL_ARRAY {1, 2, 4, 8, 15, 30, 60, 119, 237, 473, 945}
#define MISSING_NODES_BEFORE_LEVEL_ARRAY {0, 0, 0, 0, 0, 1, 3, 7, 16, 35, 74}
#define BITS_N_ZQ_CT_RNG 1361
#define BITS_BETA_ZQSTAR_CT_RNG 6730
#define BITS_V_CT_RNG 59289
#define BITS_N_ZZ_CT_RNG 729
#define BITS_CWSTR_RNG 9332

#elif ( defined(CATEGORY_5) && defined(RSDP)  && defined(SPEED) )
#define TREE_NODES_TO_STORE 166
#define NUM_NODES_SEED_TREE 658
#define NODES_PER_LEVEL_ARRAY {1, 2, 3, 6, 11, 21, 41, 82, 164, 327}
#define MISSING_NODES_BEFORE_LEVEL_ARRAY {0, 0, 0, 1, 3, 8, 19, 42, 88, 180}
#define BITS_N_ZQ_CT_RNG 1827
#define BITS_BETA_ZQSTAR_CT_RNG 2383
#define BITS_V_CT_RNG 106427
#define BITS_N_ZZ_CT_RNG 979
#define BITS_CWSTR_RNG 3044

#elif ( defined(CATEGORY_5) && defined(RSDP)  && defined(BALANCED) )
#define TREE_NODES_TO_STORE 214
#define NUM_NODES_SEED_TREE 1015
#define NODES_PER_LEVEL_ARRAY {1, 2, 4, 8, 16, 32, 64, 127, 254, 507}
#define MISSING_NODES_BEFORE_LEVEL_ARRAY {0, 0, 0, 0, 0, 0, 0, 0, 1, 3}
#define BITS_N_ZQ_CT_RNG 1827
#define BITS_BETA_ZQSTAR_CT_RNG 3658
#define BITS_V_CT_RNG 106427
#define BITS_N_ZZ_CT_RNG 979
#define BITS_CWSTR_RNG 4734

#elif ( defined(CATEGORY_5) && defined(RSDP)  && defined(SIG_SIZE) )
#define TREE_NODES_TO_STORE 231
#define NUM_NODES_SEED_TREE 1938
#define NODES_PER_LEVEL_ARRAY {1, 2, 4, 8, 16, 31, 61, 121, 242, 484, 968}
#define MISSING_NODES_BEFORE_LEVEL_ARRAY {0, 0, 0, 0, 0, 0, 1, 4, 11, 25, 53}
#define BITS_N_ZQ_CT_RNG 1827
#define BITS_BETA_ZQSTAR_CT_RNG 6914
#define BITS_V_CT_RNG 106427
#define BITS_N_ZZ_CT_RNG 979
#define BITS_CWSTR_RNG 9665

#elif ( defined(CATEGORY_1) &&  defined(RSDPG)  && defined(SPEED) )
#define TREE_NODES_TO_STORE 78
#define NUM_NODES_SEED_TREE 310
#define NODES_PER_LEVEL_ARRAY {1, 2, 3, 5, 10, 20, 39, 77, 153}
#define MISSING_NODES_BEFORE_LEVEL_ARRAY {0, 0, 0, 1, 4, 10, 22, 47, 98}
#define BITS_N_ZQ_CT_RNG 521
#define BITS_BETA_ZQSTAR_CT_RNG 1413
#define BITS_V_CT_RNG 6208
#define BITS_W_CT_RNG 5311
#define BITS_M_ZZ_CT_RNG 199
#define BITS_CWSTR_RNG 1264

#elif ( defined(CATEGORY_1) &&  defined(RSDPG)  && defined(BALANCED) )
#define TREE_NODES_TO_STORE 101
#define NUM_NODES_SEED_TREE 488
#define NODES_PER_LEVEL_ARRAY {1, 2, 4, 8, 16, 31, 61, 122, 243}
#define MISSING_NODES_BEFORE_LEVEL_ARRAY {0, 0, 0, 0, 0, 0, 1, 4, 10}
#define BITS_N_ZQ_CT_RNG 521
#define BITS_BETA_ZQSTAR_CT_RNG 2228
#define BITS_V_CT_RNG 6208
#define BITS_W_CT_RNG 5311
#define BITS_M_ZZ_CT_RNG 199
#define BITS_CWSTR_RNG 2030

#elif ( defined(CATEGORY_1) &&  defined(RSDPG)  && defined(SIG_SIZE) )
#define TREE_NODES_TO_STORE 113
#define NUM_NODES_SEED_TREE 1745
#define NODES_PER_LEVEL_ARRAY {1, 2, 4, 7, 14, 28, 55, 109, 218, 436, 871}
#define MISSING_NODES_BEFORE_LEVEL_ARRAY {0, 0, 0, 0, 1, 3, 7, 16, 35, 73, 149}
#define BITS_N_ZQ_CT_RNG 521
#define BITS_BETA_ZQSTAR_CT_RNG 7903
#define BITS_V_CT_RNG 6208
#define BITS_W_CT_RNG 5311
#define BITS_M_ZZ_CT_RNG 199
#define BITS_CWSTR_RNG 8468

#elif ( defined(CATEGORY_3) &&  defined(RSDPG)  && defined(SPEED) )
#define TREE_NODES_TO_STORE 119
#define NUM_NODES_SEED_TREE 462
#define NODES_PER_LEVEL_ARRAY {1, 2, 4, 8, 15, 29, 58, 115, 230}
#define MISSING_NODES_BEFORE_LEVEL_ARRAY {0, 0, 0, 0, 0, 1, 4, 10, 23}
#define BITS_N_ZQ_CT_RNG 751
#define BITS_BETA_ZQSTAR_CT_RNG 2125
#define BITS_V_CT_RNG 13483
#define BITS_W_CT_RNG 11025
#define BITS_M_ZZ_CT_RNG 317
#define BITS_CWSTR_RNG 2003

#elif ( defined(CATEGORY_3) &&  defined(RSDPG)  && defined(BALANCED) )
#define TREE_NODES_TO_STORE 134
#define NUM_NODES_SEED_TREE 510
#define NODES_PER_LEVEL_ARRAY {1, 2, 4, 8, 16, 32, 64, 128, 255}
#define MISSING_NODES_BEFORE_LEVEL_ARRAY {0, 0, 0, 0, 0, 0, 0, 0, 0}
#define BITS_N_ZQ_CT_RNG 751
#define BITS_BETA_ZQSTAR_CT_RNG 2351
#define BITS_V_CT_RNG 13483
#define BITS_W_CT_RNG 11025
#define BITS_M_ZZ_CT_RNG 317
#define BITS_CWSTR_RNG 2205

#elif ( defined(CATEGORY_3) &&  defined(RSDPG)  && defined(SIG_SIZE) )
#define TREE_NODES_TO_STORE 167
#define NUM_NODES_SEED_TREE 1901
#define NODES_PER_LEVEL_ARRAY {1, 2, 4, 8, 15, 30, 60, 119, 238, 475, 949}
#define MISSING_NODES_BEFORE_LEVEL_ARRAY {0, 0, 0, 0, 0, 1, 3, 7, 16, 34, 71}
#define BITS_N_ZQ_CT_RNG 751
#define BITS_BETA_ZQSTAR_CT_RNG 8627
#define BITS_V_CT_RNG 13483
#define BITS_W_CT_RNG 11025
#define BITS_M_ZZ_CT_RNG 317
#define BITS_CWSTR_RNG 9373

#elif ( defined(CATEGORY_5) &&  defined(RSDPG)  && defined(SPEED) )
#define TREE_NODES_TO_STORE 155
#define NUM_NODES_SEED_TREE 616
#define NODES_PER_LEVEL_ARRAY {1, 2, 3, 5, 10, 20, 39, 77, 153, 306}
#define MISSING_NODES_BEFORE_LEVEL_ARRAY {0, 0, 0, 1, 4, 10, 22, 47, 98, 201}
#define BITS_N_ZQ_CT_RNG 1007
#define BITS_BETA_ZQSTAR_CT_RNG 2828
#define BITS_V_CT_RNG 23112
#define BITS_W_CT_RNG 19646
#define BITS_M_ZZ_CT_RNG 385
#define BITS_CWSTR_RNG 2832

#elif ( defined(CATEGORY_5) &&  defined(RSDPG)  && defined(BALANCED) )
#define TREE_NODES_TO_STORE 183
#define NUM_NODES_SEED_TREE 715
#define NODES_PER_LEVEL_ARRAY {1, 2, 3, 6, 12, 23, 45, 89, 178, 356}
#define MISSING_NODES_BEFORE_LEVEL_ARRAY {0, 0, 0, 1, 3, 7, 16, 35, 74, 152}
#define BITS_N_ZQ_CT_RNG 1007
#define BITS_BETA_ZQSTAR_CT_RNG 3281
#define BITS_V_CT_RNG 23112
#define BITS_W_CT_RNG 19646
#define BITS_M_ZZ_CT_RNG 385
#define BITS_CWSTR_RNG 3330

#elif ( defined(CATEGORY_5) &&  defined(RSDPG)  && defined(SIG_SIZE) )
#define TREE_NODES_TO_STORE 219
#define NUM_NODES_SEED_TREE 1994
#define NODES_PER_LEVEL_ARRAY {1, 2, 4, 8, 16, 32, 63, 125, 249, 498, 996}
#define MISSING_NODES_BEFORE_LEVEL_ARRAY {0, 0, 0, 0, 0, 0, 0, 1, 4, 11, 25}
#define BITS_N_ZQ_CT_RNG 1007
#define BITS_BETA_ZQSTAR_CT_RNG 9070
#define BITS_V_CT_RNG 23112
#define BITS_W_CT_RNG 19646
#define BITS_M_ZZ_CT_RNG 385
#define BITS_CWSTR_RNG 9947

#endif
