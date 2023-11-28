#ifndef SYMMETRIC_H
#define SYMMETRIC_H

#include "fips202.h"
#include "params.h"
#include <stddef.h>

void OSKR_shake128_absorb(shake128ctx *s, const unsigned char *input, unsigned char x, unsigned char y);
void OSKR_shake128_squeezeblocks(unsigned char *output, size_t nblocks, shake128ctx *s);
void shake256_prf(unsigned char *output, size_t outlen, const unsigned char *key, unsigned char nonce);

#define hash_h(OUT, IN, INBYTES) sha3_512(OUT, IN, INBYTES)
#define hash_g(OUT, IN, INBYTES) shake256(OUT,128, IN, INBYTES)
#define xof_absorb(STATE, IN, X, Y) OSKR_shake128_absorb(STATE, IN, X, Y)
#define xof_squeezeblocks(OUT, OUTBLOCKS, STATE) OSKR_shake128_squeezeblocks(OUT, OUTBLOCKS, STATE)
#define prf(OUT, OUTBYTES, KEY, NONCE) shake256_prf(OUT, OUTBYTES, KEY, NONCE)
#define kdf(OUT, IN, INBYTES) shake256(OUT, 2*OSKR_SSBYTES, IN, INBYTES)
#define kdf1(OUT, IN, INBYTES) shake256(OUT, OSKR_SSBYTES, IN, INBYTES)
#define PREFIXLEN 65
#define XOF_BLOCKBYTES 168

typedef shake128ctx xof_state;

#endif /* SYMMETRIC_H */
