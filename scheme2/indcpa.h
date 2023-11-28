#ifndef INDCPA_H
#define INDCPA_H

#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "polyvec.h"

void pack_pk(unsigned char *r, const polyvec *pk, const unsigned char *seed);
void unpack_pk(polyvec *pk, unsigned char *seed, const unsigned char *packedpk);

void pack_ciphertext(unsigned char *r, const polyvec *b, const poly *v);
void unpack_ciphertext(polyvec *b, poly *v, const unsigned char *c);

void pack_sk(unsigned char *r, polyvec *sk);
void unpack_sk(polyvec *sk, const unsigned char *packedsk);

void gen_matrix(polyvec *a, const uint8_t seed[OKAI_SYMBYTES], int transposed);

void indcpa_keypair(unsigned char *pk, 
                   unsigned char *sk);

void indcpa_enc(unsigned char *c,
               const unsigned char *m,
               const unsigned char *pk,
               const unsigned char *coins);

void pke_enc(unsigned char *c,
             const unsigned char *m,
             const unsigned char *pk);

void indcpa_dec(unsigned char *m,
               const unsigned char *c,
               const unsigned char *sk);

#endif
