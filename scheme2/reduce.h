#ifndef REDUCE_H
#define REDUCE_H

#include <stdint.h>

#define MONT 4088 // 2^16 % Q
#define QINV 57857 // q^(-1) mod 2^16

int16_t reduce_avx(int16_t *r);

#endif
