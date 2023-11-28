#ifndef API_H
#define API_H

#include "params.h"

#define CRYPTO_SECRETKEYBYTES  OSKR_SECRETKEYBYTES
#define CRYPTO_PUBLICKEYBYTES  OSKR_PUBLICKEYBYTES
#define CRYPTO_CIPHERTEXTBYTES OSKR_CIPHERTEXTBYTES
#define CRYPTO_BYTES           OSKR_SSBYTES

#define CRYPTO_ALGNAME "OSKR768"

int crypto_kem_keypair(unsigned char *pk, unsigned char *sk);

int crypto_kem_enc(unsigned char *ct, unsigned char *ss, const unsigned char *pk);

int crypto_kem_dec(unsigned char *ss, const unsigned char *ct, const unsigned char *sk);


#endif
