/**
 *
 * Reference ISO-C11 Implementation of CROSS.
 *
 * @version 1.1 (March 2023)
 *
 * @author Alessandro Barenghi <alessandro.barenghi@polimi.it>
 * @author: Patrick Karl <patrick.karl@tum.de>
 * @author Gerardo Pelosi <gerardo.pelosi@polimi.it>
 *
 * Authors listed in lexicographic order.
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

#include "csprng_hash.h"
#if defined(NO_TREES)

#define TO_PUBLISH 1
#define NOT_TO_PUBLISH 0

void merkle_tree_root_compute(uint8_t root[HASH_DIGEST_LENGTH],
                              uint8_t leaves[T][HASH_DIGEST_LENGTH]){
   uint8_t quad_hash[4][HASH_DIGEST_LENGTH];
   int remainders[4] = {0};
   if(T%4 > 0){ remainders[0] = 1; } 
   if(T%4 > 1){ remainders[1] = 1; } 
   if(T%4 > 2){ remainders[2] = 1; }
   int offset = 0;
   for (int i = 0; i < 4; i++){
       hash(quad_hash[i],
            leaves[(T/4)*i+offset],
            (T/4+remainders[i])*HASH_DIGEST_LENGTH);
       offset += remainders[i];  
   }
   hash(root,(const unsigned char *)quad_hash,4*HASH_DIGEST_LENGTH);
}

uint16_t merkle_tree_proof_compute(uint8_t mtp[W*HASH_DIGEST_LENGTH],
                                   uint8_t leaves[T][HASH_DIGEST_LENGTH],
                                   const uint8_t leaves_to_reveal[T]){
    uint16_t published = 0;
    for(int i=0; i<T; i++){
       if(leaves_to_reveal[i] == TO_PUBLISH){
          memcpy(&mtp[HASH_DIGEST_LENGTH*published],
                 &leaves[i],
                 HASH_DIGEST_LENGTH);
          published++;
       }
    }
    return published;    
}

void merkle_tree_root_recompute(uint8_t root[HASH_DIGEST_LENGTH],
                                uint8_t recomputed_leaves[T][HASH_DIGEST_LENGTH],
                                const uint8_t mtp[W*HASH_DIGEST_LENGTH],
                                const uint8_t leaves_to_reveal[T]){
    uint16_t published = 0;
    for(int i=0; i<T; i++){
       if(leaves_to_reveal[i] == TO_PUBLISH){
          memcpy(&recomputed_leaves[i],
                 &mtp[HASH_DIGEST_LENGTH*published],
                 HASH_DIGEST_LENGTH);
          published++;
       }
    }
    merkle_tree_root_compute(root,recomputed_leaves);
}
#else

#include "merkle.h"
/* Stub of the interface to Merkle tree root computer from all leaves */
void merkle_tree_root_compute(uint8_t root[HASH_DIGEST_LENGTH],
                              uint8_t tree[NUM_NODES_MERKLE_TREE * HASH_DIGEST_LENGTH],
                              /* input, although mutable in caller, having as const is non
                               * tolerated in strict ISO C */
                              uint8_t leaves[T][HASH_DIGEST_LENGTH]){
    generate_merkle_tree(tree, leaves);
    /* Root is at first position of the tree */
    memcpy(root, tree, HASH_DIGEST_LENGTH);
}

/* Stub interface to the function computing the Merkle tree proof, storing it
 * in the signature. Returns the number of digests in the merkle tree proof */
uint16_t merkle_tree_proof_compute(uint8_t mtp[HASH_DIGEST_LENGTH*TREE_NODES_TO_STORE],
                               const uint8_t tree[NUM_NODES_MERKLE_TREE*HASH_DIGEST_LENGTH],
                               uint8_t leaves[T][HASH_DIGEST_LENGTH],
                               const uint8_t leaves_to_reveal[T]){
    uint16_t mtp_len; 
    uint16_t merkle_proof_indices[TREE_NODES_TO_STORE];

    generate_merkle_proof(merkle_proof_indices, &mtp_len, leaves_to_reveal);

    for (size_t i=0; i< mtp_len; i++) {
        memcpy(mtp + i*HASH_DIGEST_LENGTH, tree + merkle_proof_indices[i]*HASH_DIGEST_LENGTH,
                HASH_DIGEST_LENGTH);
    }
    return mtp_len;
}

/* stub of the interface to Merkle tree recomputation given the proof and
 * the computed leaves */
void merkle_tree_root_recompute(uint8_t root[HASH_DIGEST_LENGTH],
                                uint8_t recomputed_leaves[T][HASH_DIGEST_LENGTH],
                                const uint8_t mtp[HASH_DIGEST_LENGTH*TREE_NODES_TO_STORE],
                                const uint8_t leaves_to_reveal[T]){

    unsigned char tree[NUM_NODES_MERKLE_TREE * HASH_DIGEST_LENGTH];

    rebuild_merkle_tree(tree, mtp, recomputed_leaves, leaves_to_reveal);
    memcpy(root, tree, HASH_DIGEST_LENGTH);
}
#endif
