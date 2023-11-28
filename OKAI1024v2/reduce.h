#ifndef REDUCE_H
#define REDUCE_H

#include <stdint.h>

#define MONT 4088 // 2^16 mod q
#define QINV 57857 // q^-1 mod 2^16

int16_t freeze(int16_t x);

int16_t montgomery_reduce(int32_t a);

int16_t barrett_reduce(int16_t a);

int16_t csubq(int16_t a);

#endif
