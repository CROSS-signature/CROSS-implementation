#pragma once

#include "m4_defs.h"

/* fallback functions to test M4-specific features (e.g. CMSIS) on x86-64 */

uint32_t __smlad(const uint32_t a, const uint32_t b, const uint32_t c){
    return c + ((a & 0xFFFF) * (b & 0xFFFF)) + (((a >> 16) & 0xFFFF) * ((b >> 16) & 0xFFFF));
}

uint32_t smlad_batch_V(const uint16_t *a, const uint16_t *b){
    uint32_t res = 0;
    uint32_t rep = SMLAD_BATCH_SIZE_V;
    for(int i = 0; i < rep; i++){
        res += a[i] * b[i];
    }
    return res;
}

#if defined(RSDPG)
uint32_t smlad_batch_W(const uint16_t *a, const uint16_t *b){
    uint32_t res = 0;
    uint32_t rep = SMLAD_BATCH_SIZE_W;
    for(int i = 0; i < rep; i++){
        res += a[i] * b[i];
    }
    return res;
}
#endif
