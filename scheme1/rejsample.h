#ifndef REJSAMPLE_H
#define REJSAMPLE_H

#include <stdint.h>
#include "params.h"
#include "symmetric.h"

#define REJ_UNIFORM_AVX_NBLOCKS ((12*OSKR_N/8*(1 << 12)/OSKR_Q + XOF_BLOCKBYTES)/XOF_BLOCKBYTES)//6
#define REJ_UNIFORM_AVX_BUFLEN (REJ_UNIFORM_AVX_NBLOCKS*XOF_BLOCKBYTES)//1008

unsigned int rej_uniform_avx(int16_t *r, const uint8_t *buf);

#endif
