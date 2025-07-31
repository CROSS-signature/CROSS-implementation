#include "csprng_hash.h"

/* global csprng state employed to have deterministic randombytes for testing */
CSPRNG_STATE_T platform_csprng_state;
/* extracts xlen bytes from the global CSPRNG */
void randombytes(unsigned char * x,
                 unsigned long long xlen) {
   csprng_randombytes(x,xlen,&platform_csprng_state);
}
