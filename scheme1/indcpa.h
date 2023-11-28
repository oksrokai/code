#ifndef INDCPA_H
#define INDCPA_H

#include <stdint.h>
#include "params.h"
#include "polyvec.h"

void gen_matrix(polyvec *a, const uint8_t seed[OSKR_SYMBYTES], int transposed);
void indcpa_keypair(uint8_t pk[OSKR_INDCPA_PUBLICKEYBYTES],
                    uint8_t sk[OSKR_INDCPA_SECRETKEYBYTES]);

void indcpa_enc(uint8_t c[OSKR_INDCPA_BYTES],
                const uint8_t m[OSKR_INDCPA_MSGBYTES],
                const uint8_t pk[OSKR_INDCPA_PUBLICKEYBYTES],
                const uint8_t coins[OSKR_SYMBYTES]);

void indcpa_dec(uint8_t m[OSKR_INDCPA_MSGBYTES],
                const uint8_t c[OSKR_INDCPA_BYTES],
                const uint8_t sk[OSKR_INDCPA_SECRETKEYBYTES]);

void unpack_sk(polyvec *sk, const uint8_t packedsk[OSKR_INDCPA_SECRETKEYBYTES]);
void unpack_pk(polyvec *pk,
                      uint8_t seed[OSKR_SYMBYTES],
                      const uint8_t packedpk[OSKR_INDCPA_PUBLICKEYBYTES]);

#endif
