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
void poly_basemul_256(poly *r, const poly *a, const poly *b);
void poly_basemul(poly *r, const poly *a, const poly *b);
void poly_basemul_acc(poly *r, const poly *a, const poly *b);
void poly_reduce(poly *r);
void poly_multi_basemul(poly *r, const poly *a, const poly *b);

void poly_add(poly *r, const poly *a, const poly *b);
void poly_sub(poly *r, const poly *a, const poly *b);
void poly_csubq(poly *r);

void poly_compress_5(unsigned char *r, poly *a);
void poly_decompress_5(poly *r, unsigned char *a);
void poly_zeroize(poly *p);

void poly_pointwiseyhat(int16_t r[256], const int16_t a[256]);
void poly_pointwise(int16_t r[256], const int16_t a[256], const int16_t b[256]);
void basemulyhat(int16_t r[2], const int16_t a[2], int16_t zeta);
void basemul(int16_t r[2], const int16_t a[2], const int16_t b[2], int16_t zeta);
void poly_single_basemul(int16_t *r,int16_t *a,int16_t *b,int16_t zeta);

void poly_con_4(uint8_t c[128],
              const uint8_t m[32],
              int16_t *s);
void poly_rec_4(uint8_t m[32],
              const uint8_t c[128],
              int16_t *s);
void poly_con(uint8_t *c,
              const uint8_t *m,
              poly *s);
void poly_rec(uint8_t *m,
              const uint8_t *c,
              poly *s);
void poly_compress_10(uint8_t r[320], const int16_t *acoeffs);
void poly_decompress_10(int16_t *rcoeffs, const uint8_t a[320]);
void poly_compress512(uint8_t r[320*2], const poly *a);
void poly_decompress512(poly *r, const uint8_t a[320*2]);
#endif
