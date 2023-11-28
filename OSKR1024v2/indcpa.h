#ifndef INDCPA_H
#define INDCPA_H
#include "params.h"
#include "polyvec.h"
#include <stdint.h>
void indcpa_keypair(unsigned char *pk,
                    unsigned char *sk);

void indcpa_enc(unsigned char *c,
                const unsigned char *m,
                const unsigned char *pk,
                const unsigned char *coins);

void indcpa_dec(unsigned char *m,
                const unsigned char *c,
                const unsigned char *sk);

void unpack_sk(polyvec *sk, const uint8_t packedsk[OSKR_INDCPA_SECRETKEYBYTES]);

#endif
