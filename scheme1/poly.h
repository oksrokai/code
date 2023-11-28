#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include "align.h"
#include "params.h"

typedef ALIGNED_INT16(OSKR_N) poly;

void poly_con(uint8_t c[OSKR_POLYCOMPRESSEDBYTES], const uint8_t m[OSKR_INDCPA_MSGBYTES], poly *s);
void poly_rec(uint8_t m[OSKR_INDCPA_MSGBYTES], const uint8_t c[OSKR_POLYCOMPRESSEDBYTES], poly *s);

void con_avx_3(uint8_t c[OSKR_POLYCOMPRESSEDBYTES], const uint8_t m[OSKR_INDCPA_MSGBYTES], poly *s, const __m256i *qdata);
void rec_avx_3(uint8_t m[OSKR_INDCPA_MSGBYTES], const uint8_t c[OSKR_POLYCOMPRESSEDBYTES], poly *s, const __m256i *qdata);

void con_avx_4(uint8_t c[OSKR_POLYCOMPRESSEDBYTES], const uint8_t m[OSKR_INDCPA_MSGBYTES], poly *s, const __m256i *qdata);
void rec_avx_4(uint8_t m[OSKR_INDCPA_MSGBYTES], const uint8_t c[OSKR_POLYCOMPRESSEDBYTES], poly *s, const __m256i *qdata);

void con_avx_5(uint8_t c[OSKR_POLYCOMPRESSEDBYTES], const uint8_t m[OSKR_INDCPA_MSGBYTES], __m256i *s, const __m256i *qdata);
void rec_avx_5(uint8_t m[OSKR_INDCPA_MSGBYTES], const uint8_t c[OSKR_POLYCOMPRESSEDBYTES], __m256i *s, const __m256i *qdata);

void poly_compress_to_10_avx2(uint8_t r[OSKR_POLYCOMPRESSEDBYTES], poly *a, const __m256i *qdata);
void poly_decompress_from_10_avx2(poly *r, const uint8_t a[OSKR_POLYCOMPRESSEDBYTES], const __m256i *qdata);

void poly_compress_to_11_avx2(uint8_t r[OSKR_POLYCOMPRESSEDBYTES], int16_t *a, const __m256i *qdata);
void poly_decompress_from_11_avx2(int16_t *r, const uint8_t a[OSKR_POLYCOMPRESSEDBYTES], const __m256i *qdata);

void poly_tobytes(uint8_t r[OSKR_POLYBYTES], const poly *a);
void poly_frombytes(poly *r, const uint8_t a[OSKR_POLYBYTES]);

void poly_getnoise_eta1(poly *r, const uint8_t seed[OSKR_SYMBYTES], uint8_t nonce);

void poly_getnoise_eta2(poly *r, const uint8_t seed[OSKR_SYMBYTES], uint8_t nonce);

#ifndef OSKR_90S
void poly_getnoise_eta1_4x(poly *r0,
                           poly *r1,
                           poly *r2,
                           poly *r3,
                           const uint8_t seed[32],
                           uint8_t nonce0,
                           uint8_t nonce1,
                           uint8_t nonce2,
                           uint8_t nonce3);

#if OSKR_K == 2
void poly_getnoise_eta1122_4x(poly *r0,
                              poly *r1,
                              poly *r2,
                              poly *r3,
                              const uint8_t seed[32],
                              uint8_t nonce0,
                              uint8_t nonce1,
                              uint8_t nonce2,
                              uint8_t nonce3);
#endif
#endif

void poly_ntt(poly *r);
void poly_invntt_tomont(poly *r);
void poly_nttunpack(poly *r);
void poly_basemul_montgomery(poly *r, const poly *a, const poly *b);
void poly_tomont(poly *r);
void poly_reduce(poly *r);
void poly_add(poly *r, const poly *a, const poly *b);
void poly_sub(poly *r, const poly *a, const poly *b);
void poly_basemul_montgomery_multi(poly *r, const poly *a, const poly *b);

#endif
