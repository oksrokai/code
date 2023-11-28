#ifndef SYMMETRIC_H
#define SYMMETRIC_H

#include "params.h"
#include "fips202.h"

typedef struct {
  uint64_t s[25];
} keccak_state;

void OKAI_shake128_absorb(keccak_state *s, const unsigned char *input, unsigned char x, unsigned char y);
void OKAI_shake128_squeezeblocks(unsigned char *output, unsigned long long nblocks, keccak_state *s);
void shake256_prf(unsigned char *output, unsigned long long outlen, const unsigned char *key, const unsigned char nonce);

#define hash_h(OUT, IN, INBYTES) sha3_512(OUT, IN, INBYTES)
#define hash_g(OUT, IN, INBYTES) shake256(OUT, 128, IN, INBYTES)
#define xof_absorb(STATE, IN, X, Y) OKAI_shake128_absorb(STATE, IN, X, Y)
#define xof_squeezeblocks(OUT, OUTBLOCKS, STATE) OKAI_shake128_squeezeblocks(OUT, OUTBLOCKS, STATE)
#define prf(OUT, OUTBYTES, KEY, NONCE) shake256_prf(OUT, OUTBYTES, KEY, NONCE)
#define kdf(OUT, IN, INBYTES) shake256(OUT, 2*OKAI_SYMBYTES, IN, INBYTES)
#define kdf1(OUT, IN, INBYTES) shake256(OUT, OKAI_SYMBYTES, IN, INBYTES)
#define PREFIXLEN 65
#define XOF_BLOCKBYTES 168

typedef keccak_state xof_state;

#endif /* SYMMETRIC_H */
