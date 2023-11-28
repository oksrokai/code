#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include <immintrin.h>
#include "params.h"

/*
 * Elements of R_q = Z_q[X]/(X^n + 1). Represents polynomial
 * coeffs[0] + X*coeffs[1] + X^2*xoeffs[2] + ... + X^{n-1}*coeffs[n-1]
 */
typedef struct{
  int16_t __attribute__((aligned(32))) coeffs[OKAI_N];
} poly;

void poly_con(uint8_t c[OKAI_POLYCOMPRESSEDBYTES], const uint8_t m[OKAI_INDCPA_MSGBYTES], poly *s);
void poly_rec(uint8_t m[OKAI_INDCPA_MSGBYTES], const uint8_t c[OKAI_POLYCOMPRESSEDBYTES], poly *s);

void con_avx_4(uint8_t c[OKAI_POLYCOMPRESSEDBYTES], const uint8_t m[OKAI_INDCPA_MSGBYTES], int16_t *s);
void rec_avx_4(uint8_t m[OKAI_INDCPA_MSGBYTES], const uint8_t c[OKAI_POLYCOMPRESSEDBYTES], int16_t *s);

void poly_compress(uint8_t *r, const poly *a);
void poly_decompress(poly *r, const uint8_t *a);

void poly_compress_to_8_avx2(uint8_t *r, const poly *a);
void poly_decompress_from_8_avx2(poly *r, const uint8_t *a);

void poly_compress_to_9_avx2(uint8_t *r, const poly *a);
void poly_decompress_from_9_avx2(poly *r, const uint8_t *a);

void poly_compress_to_10_avx2(uint8_t *r, const int16_t *a);
void poly_decompress_from_10_avx2(int16_t *r, const uint8_t *a);

void poly_tobytes(unsigned char *r, const poly *a);
void poly_tobytes_avx2(uint8_t *r, const int16_t *a);

void poly_frombytes(poly *r, const unsigned char *a);
void poly_frombytes_avx2(int16_t *r, const uint8_t *a);

void poly_frommsg(poly *r, const uint8_t msg[OKAI_SYMBYTES]);
void poly_tomsg(uint8_t msg[OKAI_SYMBYTES], poly *r);

void poly_getnoises(poly *r, const uint8_t seed[OKAI_SYMBYTES], uint8_t nonce);
void poly_getnoisee(poly *r, const uint8_t seed[OKAI_SYMBYTES], uint8_t nonce);

void poly_ntt(poly *r);
void poly_invntt(poly *r);

void poly_nttpack(poly *r);
void poly_nttunpack(poly *r);

void poly_reduce(poly *r);

void poly_add(poly *r, const poly *a, const poly *b);
void poly_sub(poly *r, const poly *a, const poly *b);

void poly_basemul_montgomery_multi(poly *r, const poly *a, const poly *b);
#endif
