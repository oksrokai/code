#include <stdio.h>
#include <stdint.h>
#include <immintrin.h>
#include <string.h>
#include "align.h"
#include "params.h"
#include "poly.h"
#include "ntt.h"
#include "consts.h"
#include "reduce.h"
#include "cbd.h"
#include "symmetric.h"

void poly_con(uint8_t c[OSKR_POLYCOMPRESSEDBYTES],
              const uint8_t m[OSKR_INDCPA_MSGBYTES],
              poly *s)
{
#if (OSKR_N == 256)
  con_avx_4(c, m, s, qdata.vec);
#elif (OSKR_N == 512)
  con_avx_5(c, m, &s->vec[0], qdata.vec);
  con_avx_5(c + OSKR_POLYCOMPRESSEDBYTES/2, m + OSKR_INDCPA_MSGBYTES/2, &s->vec[16], qdata.vec);
#endif
}

void poly_rec(uint8_t m[OSKR_INDCPA_MSGBYTES],
              const uint8_t c[OSKR_POLYCOMPRESSEDBYTES],
              poly *s)
{
#if (OSKR_N == 256)
  rec_avx_4(m, c, s, qdata.vec);
#elif (OSKR_N == 512)
  rec_avx_5(m, c, &s->vec[0], qdata.vec);
  rec_avx_5(m + OSKR_INDCPA_MSGBYTES/2, c + OSKR_POLYCOMPRESSEDBYTES/2, &s->vec[16], qdata.vec);
#endif
}

/*************************************************
* Name:        poly_tobytes
*
* Description: Serialization of a polynomial in NTT representation.
*              The coefficients of the input polynomial are assumed to
*              lie in the invertal [0,q], i.e. the polynomial must be reduced
*              by poly_reduce(). The coefficients are orderd as output by
*              poly_ntt(); the serialized output coefficients are in bitreversed
*              order.
*
* Arguments:   - uint8_t *r: pointer to output byte array
*                            (needs space for OSKR_POLYBYTES bytes)
*              - poly *a: pointer to input polynomial
**************************************************/
void poly_tobytes(uint8_t r[OSKR_POLYBYTES], const poly *a)
{
  ntttobytes_avx(r, a->vec, qdata.vec);
#if (OSKR_N == 512)
  ntttobytes_avx(r + OSKR_POLYBYTES/2, &a->vec[16], qdata.vec);
#endif
}

/*************************************************
* Name:        poly_frombytes
*
* Description: De-serialization of a polynomial;
*              inverse of poly_tobytes
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *a: pointer to input byte array
*                                  (of OSKR_POLYBYTES bytes)
**************************************************/
void poly_frombytes(poly *r, const uint8_t a[OSKR_POLYBYTES])
{
  nttfrombytes_avx(r->vec, a, qdata.vec);
#if (OSKR_N == 512)
  nttfrombytes_avx(&r->vec[16], a + OSKR_POLYBYTES/2, qdata.vec);
#endif
}

/*************************************************
* Name:        poly_getnoise_eta1
*
* Description: Sample a polynomial deterministically from a seed and a nonce,
*              with output polynomial close to centered binomial distribution
*              with parameter OSKR_ETA1
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *seed: pointer to input seed
*                                     (of length OSKR_SYMBYTES bytes)
*              - uint8_t nonce: one-byte input nonce
**************************************************/
void poly_getnoise_eta1(poly *r, const uint8_t seed[OSKR_SYMBYTES], uint8_t nonce)
{
  ALIGNED_UINT8(OSKR_ETA1*OSKR_N/4+32) buf;
  prf(buf.coeffs, OSKR_ETA1*OSKR_N/4, seed, nonce);
  poly_cbd_eta1(r, buf.vec);
}

/*************************************************
* Name:        poly_getnoise_eta2
*
* Description: Sample a polynomial deterministically from a seed and a nonce,
*              with output polynomial close to centered binomial distribution
*              with parameter OSKR_ETA2
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *seed: pointer to input seed
*                                     (of length OSKR_SYMBYTES bytes)
*              - uint8_t nonce: one-byte input nonce
**************************************************/
void poly_getnoise_eta2(poly *r, const uint8_t seed[OSKR_SYMBYTES], uint8_t nonce)
{
  ALIGNED_UINT8(OSKR_ETA2*OSKR_N/4) buf;
  prf(buf.coeffs, OSKR_ETA2*OSKR_N/4, seed, nonce);
  poly_cbd_eta2(r, buf.vec);
}

#define NOISE_NBLOCKS ((OSKR_ETA1*OSKR_N/4+SHAKE256_RATE-1)/SHAKE256_RATE)
void poly_getnoise_eta1_4x(poly *r0,
                           poly *r1,
                           poly *r2,
                           poly *r3,
                           const uint8_t seed[OSKR_SYMBYTES],
                           uint8_t nonce0,
                           uint8_t nonce1,
                           uint8_t nonce2,
                           uint8_t nonce3)
{
  ALIGNED_UINT8(NOISE_NBLOCKS*SHAKE256_RATE) buf[4];
  __m256i f;
  keccakx4_state state;

  f = _mm256_loadu_si256((__m256i *)seed);
  _mm256_store_si256(buf[0].vec, f);
  _mm256_store_si256(buf[1].vec, f);
  _mm256_store_si256(buf[2].vec, f);
  _mm256_store_si256(buf[3].vec, f);
#if (OSKR_N == 512)
  f = _mm256_loadu_si256((__m256i *)(&seed[32]));
  _mm256_store_si256(&buf[0].vec[1], f);
  _mm256_store_si256(&buf[1].vec[1], f);
  _mm256_store_si256(&buf[2].vec[1], f);
  _mm256_store_si256(&buf[3].vec[1], f);
#endif
  buf[0].coeffs[OSKR_SYMBYTES] = nonce0;
  buf[1].coeffs[OSKR_SYMBYTES] = nonce1;
  buf[2].coeffs[OSKR_SYMBYTES] = nonce2;
  buf[3].coeffs[OSKR_SYMBYTES] = nonce3;

  shake256x4_absorb_once(&state, buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, OSKR_SYMBYTES+1);
  shake256x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, NOISE_NBLOCKS, &state);

  poly_cbd_eta1(r0, buf[0].vec);
  poly_cbd_eta1(r1, buf[1].vec);
  poly_cbd_eta1(r2, buf[2].vec);
  poly_cbd_eta1(r3, buf[3].vec);
}

void poly_getnoise_eta1122_4x(poly *r0,
                              poly *r1,
                              poly *r2,
                              poly *r3,
                              const uint8_t seed[32],
                              uint8_t nonce0,
                              uint8_t nonce1,
                              uint8_t nonce2,
                              uint8_t nonce3)
{
  ALIGNED_UINT8(NOISE_NBLOCKS*SHAKE256_RATE) buf[4];
  __m256i f;
  keccakx4_state state;

  f = _mm256_loadu_si256((__m256i *)seed);
  _mm256_store_si256(buf[0].vec, f);
  _mm256_store_si256(buf[1].vec, f);
  _mm256_store_si256(buf[2].vec, f);
  _mm256_store_si256(buf[3].vec, f);

  buf[0].coeffs[32] = nonce0;
  buf[1].coeffs[32] = nonce1;
  buf[2].coeffs[32] = nonce2;
  buf[3].coeffs[32] = nonce3;

  shake256x4_absorb_once(&state, buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 33);
  shake256x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, NOISE_NBLOCKS, &state);

  poly_cbd_eta1(r0, buf[0].vec);
  poly_cbd_eta1(r1, buf[1].vec);
  poly_cbd_eta2(r2, buf[2].vec);
  poly_cbd_eta2(r3, buf[3].vec);
}

/*************************************************
* Name:        poly_ntt
*
* Description: Computes negacyclic number-theoretic transform (NTT) of
*              a polynomial in place.
*              Input coefficients assumed to be in normal order,
*              output coefficients are in special order that is natural
*              for the vectorization. Input coefficients are assumed to be
*              bounded by q in absolute value, output coefficients are bounded
*              by 16118 in absolute value.
*
* Arguments:   - poly *r: pointer to in/output polynomial
**************************************************/
void poly_ntt(poly *r)
{
  ntt_avx(r->vec, qdata.vec);
#if (OSKR_N == 512)
  ntt_avx(&r->vec[16], qdata.vec);
#endif
}

/*************************************************
* Name:        poly_invntt_tomont
*
* Description: Computes inverse of negacyclic number-theoretic transform (NTT)
*              of a polynomial in place;
*              Input coefficients assumed to be in special order from vectorized
*              forward ntt, output in normal order. Input coefficients can be
*              arbitrary 16-bit integers, output coefficients are bounded by 14870
*              in absolute value.
*
* Arguments:   - poly *a: pointer to in/output polynomial
**************************************************/
void poly_invntt_tomont(poly *r)
{
  invntt_avx(r->vec, qdata.vec);
#if (OSKR_N == 512)
  poly t;
  invntt_avx(&r->vec[16], qdata.vec);
  combine_avx(r->vec,&t.vec);
  r=&t;
#endif
}

void poly_nttunpack(poly *r)
{
  nttunpack_avx(r->vec, qdata.vec);
#if (OSKR_N == 512)
  nttunpack_avx(&r->vec[16], qdata.vec);
#endif
}

/*************************************************
* Name:        poly_basemul_montgomery
*
* Description: Multiplication of two polynomials in NTT domain.
*              One of the input polynomials needs to have coefficients
*              bounded by q, the other polynomial can have arbitrary
*              coefficients. Output coefficients are bounded by 6656.
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_basemul_montgomery(poly *r, const poly *a, const poly *b)
{
  basemul_avx(r->vec, a->vec, b->vec, qdata.vec);
}

void poly_basemul_montgomery_multi(poly *r, const poly *a, const poly *b)
{
  basemul_multi_avx(a->vec,b->vec,r->vec,qdata.vec);
}
/*************************************************
* Name:        poly_tomont
*
* Description: Inplace conversion of all coefficients of a polynomial
*              from normal domain to Montgomery domain
*
* Arguments:   - poly *r: pointer to input/output polynomial
**************************************************/
void poly_tomont(poly *r)
{
  tomont_avx(r->vec, qdata.vec);
#if (OSKR_N == 512)
  tomont_avx(&r->vec[16], qdata.vec);
#endif
}

/*************************************************
* Name:        poly_reduce
*
* Description: Applies Barrett reduction to all coefficients of a polynomial
*              for details of the Barrett reduction see comments in reduce.c
*
* Arguments:   - poly *r: pointer to input/output polynomial
**************************************************/
void poly_reduce(poly *r)
{
  reduce_avx(r->vec, qdata.vec);
#if (OSKR_N == 512)
  reduce_avx(&r->vec[16], qdata.vec);
#endif
}

/*************************************************
* Name:        poly_add
*
* Description: Add two polynomials. No modular reduction
*              is performed.
*
* Arguments: - poly *r: pointer to output polynomial
*            - const poly *a: pointer to first input polynomial
*            - const poly *b: pointer to second input polynomial
**************************************************/
void poly_add(poly *r, const poly *a, const poly *b)
{
  unsigned int i;
  __m256i f0, f1;

  for(i=0;i<OSKR_N/16;i++) {
    f0 = _mm256_load_si256(&a->vec[i]);
    f1 = _mm256_load_si256(&b->vec[i]);
    f0 = _mm256_add_epi16(f0, f1);
    _mm256_store_si256(&r->vec[i], f0);
  }
}

/*************************************************
* Name:        poly_sub
*
* Description: Subtract two polynomials. No modular reduction
*              is performed.
*
* Arguments: - poly *r: pointer to output polynomial
*            - const poly *a: pointer to first input polynomial
*            - const poly *b: pointer to second input polynomial
**************************************************/
void poly_sub(poly *r, const poly *a, const poly *b)
{
  unsigned int i;
  __m256i f0, f1;

  for(i=0;i<OSKR_N/16;i++) {
    f0 = _mm256_load_si256(&a->vec[i]);
    f1 = _mm256_load_si256(&b->vec[i]);
    f0 = _mm256_sub_epi16(f0, f1);
    _mm256_store_si256(&r->vec[i], f0);
  }
}