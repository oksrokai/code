#ifndef API_H
#define API_H

#include "params.h"

#define CRYPTO_SECRETKEYBYTES  OKAI_SECRETKEYBYTES
#define CRYPTO_PUBLICKEYBYTES  OKAI_PUBLICKEYBYTES
#define CRYPTO_CIPHERTEXTBYTES OKAI_CIPHERTEXTBYTES
#define CRYPTO_BYTES           OKAI_SYMBYTES

#define CRYPTO_ALGNAME "OKAI-MLWE"

int crypto_kem_keypair(unsigned char *pk, unsigned char *sk);

int crypto_kem_enc(unsigned char *ct, unsigned char *ss, const unsigned char *pk);

int crypto_kem_dec(unsigned char *ss, const unsigned char *ct, const unsigned char *sk);


#endif
