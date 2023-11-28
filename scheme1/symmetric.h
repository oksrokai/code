#ifndef SYMMETRIC_H
#define SYMMETRIC_H

#include <stddef.h>
#include <stdint.h>
#include "params.h"

#include "fips202.h"
#include "fips202x4.h"

typedef keccak_state xof_state;

void OSKR_shake128_absorb(keccak_state *s,
                           const uint8_t seed[OSKR_SYMBYTES],
                           uint8_t x,
                           uint8_t y);

void OSKR_shake256_prf(uint8_t *out,
                        size_t outlen,
                        const uint8_t key[OSKR_SYMBYTES],
                        uint8_t nonce);

#define XOF_BLOCKBYTES SHAKE128_RATE

#if (OSKR_N == 256)
#define hash_h(OUT, IN, INBYTES) sha3_256(OUT, IN, INBYTES)
#define hash_g(OUT, IN, INBYTES) sha3_512(OUT, IN, INBYTES)
#define PREFIXLEN 33
#ifdef CCS21
#define kdf(OUT, IN, INBYTES) shake256(OUT, 2*OSKR_SSBYTES, IN, INBYTES)
#define kdf1(OUT, IN, INBYTES) shake256(OUT, OSKR_SSBYTES, IN, INBYTES)
#else
#define kdf(OUT, IN, INBYTES) shake256(OUT, OSKR_SSBYTES, IN, INBYTES)
#endif
#elif (OSKR_N == 512)
#define hash_h(OUT, IN, INBYTES) sha3_512(OUT, IN, INBYTES)
#define hash_g(OUT, IN, INBYTES) shake256(OUT, 128, IN, INBYTES)
#define PREFIXLEN 65
#ifdef CCS21
#define kdf(OUT, IN, INBYTES) shake256(OUT, 2*OSKR_SSBYTES, IN, INBYTES)
#define kdf1(OUT, IN, INBYTES) shake256(OUT, OSKR_SSBYTES, IN, INBYTES)
#else
#define kdf(OUT, IN, INBYTES) sha3_512(OUT, IN, INBYTES)
#endif
#endif

#define xof_absorb(STATE, SEED, X, Y) OSKR_shake128_absorb(STATE, SEED, X, Y)
#define xof_squeezeblocks(OUT, OUTBLOCKS, STATE) \
        shake128_squeezeblocks(OUT, OUTBLOCKS, STATE)
#define prf(OUT, OUTBYTES, KEY, NONCE) \
        OSKR_shake256_prf(OUT, OUTBYTES, KEY, NONCE)

#endif /* SYMMETRIC_H */