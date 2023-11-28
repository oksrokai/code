#ifndef KEM_H
#define KEM_H

#include <stdint.h>
#include "params.h"

#define CRYPTO_SECRETKEYBYTES  OSKR_SECRETKEYBYTES
#define CRYPTO_PUBLICKEYBYTES  OSKR_PUBLICKEYBYTES
#define CRYPTO_CIPHERTEXTBYTES OSKR_CIPHERTEXTBYTES
#define CRYPTO_BYTES           OSKR_SSBYTES

#if   (OSKR_K == 2)
#define CRYPTO_ALGNAME "OSKR512"
#elif (OSKR_K == 3)
#define CRYPTO_ALGNAME "OSKR768"
#elif (OSKR_K == 4)
#define CRYPTO_ALGNAME "OSKR1024"
#endif

int crypto_kem_keypair(uint8_t *pk, uint8_t *sk);

int crypto_kem_enc(uint8_t *ct, uint8_t *ss, const uint8_t *pk);

int crypto_kem_dec(uint8_t *ss, const uint8_t *ct, const uint8_t *sk);

#endif
