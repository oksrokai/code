#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include "params.h"

/*
 * Elements of R_q = Z_q[X]/(X^n + 1). Represents polynomial
 * coeffs[0] + X*coeffs[1] + X^2*xoeffs[2] + ... + X^{n-1}*coeffs[n-1]
 */
typedef struct{
  int16_t coeffs[OKAI_N];
} poly;

void poly_con(uint8_t c[OKAI_POLYCOMPRESSEDBYTES], const uint8_t m[OKAI_INDCPA_MSGBYTES], poly *s);
void poly_rec(uint8_t m[OKAI_INDCPA_MSGBYTES], const uint8_t c[OKAI_POLYCOMPRESSEDBYTES], poly *s);

void poly_compress(unsigned char *r, poly *a);
void poly_decompress(poly *r, const unsigned char *a);

void poly_tobytes(unsigned char *r, const poly *a);
void poly_frombytes(poly *r, const unsigned char *a);

void poly_frommsg(poly *r, const unsigned char msg[OKAI_SYMBYTES]);
void poly_tomsg(unsigned char msg[OKAI_SYMBYTES], poly *r);

void poly_getnoise(poly *r,const unsigned char *seed, unsigned char nonce);
void poly_getsecret(poly *r,const unsigned char *seed, unsigned char nonce);

void poly_ntt(poly *r);
void poly_invntt(poly *r);
void poly_basemul(poly *r, const poly *a, const poly *b);
void poly_basemul_acc(poly *r, const poly *a, const poly *b);
void poly_reduce(poly *r);

void poly_add(poly *r, const poly *a, const poly *b);
void poly_sub(poly *r, const poly *a, const poly *b);
void poly_csubq(poly *r);

void poly_compresspk(uint8_t *r, const poly *a);
void poly_decompresspk(poly *r, const uint8_t *a);
void poly_compressc(uint8_t *r, const poly *a);
void poly_decompressc(poly *r, const uint8_t *a);
void poly_zeroize(poly *p);
void poly_con_4(uint8_t c[OKAI_POLYCOMPRESSEDBYTES],
                       const uint8_t m[OKAI_INDCPA_MSGBYTES],
                       poly *s);
void poly_rec_4(uint8_t m[OKAI_INDCPA_MSGBYTES],
                       const uint8_t c[OKAI_POLYCOMPRESSEDBYTES],
                       poly *s);
#endif
