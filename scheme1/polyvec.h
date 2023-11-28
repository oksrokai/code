#ifndef POLYVEC_H
#define POLYVEC_H

#include <stdint.h>
#include "params.h"
#include "poly.h"

typedef struct{
  poly vec[OSKR_K];
} polyvec;

void polyvec_compress(uint8_t r[OSKR_POLYVECCOMPRESSEDBYTES+2], const polyvec *a);
void polyvec_decompress(polyvec *r, const uint8_t a[OSKR_POLYVECCOMPRESSEDBYTES+12]);

void polyvec_tobytes(uint8_t r[OSKR_POLYVECBYTES], const polyvec *a);
void polyvec_frombytes(polyvec *r, const uint8_t a[OSKR_POLYVECBYTES]);

void polyvec_ntt(polyvec *r);
void polyvec_invntt_tomont(polyvec *r);
void polyvec_basemul_acc_montgomery(poly *r, const polyvec *a, const polyvec *b);
void polyvec_reduce(polyvec *r);
void polyvec_add(polyvec *r, const polyvec *a, const polyvec *b);
#endif
