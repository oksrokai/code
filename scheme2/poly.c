#include <stdint.h>
#include <stdio.h>
#include "poly.h"
#include "ntt.h"
#include "polyvec.h"
#include "reduce.h"
#include "cbd.h"
#include "symmetric.h"
#include "consts.h"

void poly_con(uint8_t c[OKAI_POLYCOMPRESSEDBYTES], const uint8_t m[OKAI_INDCPA_MSGBYTES], poly *s)
{
  con_avx_4(c, m, &s->coeffs[0]);
#if (OKAI_N == 512)
  con_avx_4(c + (OKAI_POLYCOMPRESSEDBYTES >> 1), m + (OKAI_INDCPA_MSGBYTES >> 1), &s->coeffs[256]);
#endif
}

void poly_rec(uint8_t m[OKAI_INDCPA_MSGBYTES], const uint8_t c[OKAI_POLYCOMPRESSEDBYTES], poly *s)
{
  rec_avx_4(m, c, &s->coeffs[0]);
#if (OKAI_N == 512)
  rec_avx_4(m + OKAI_INDCPA_MSGBYTES/2, c + OKAI_POLYCOMPRESSEDBYTES/2, &s->coeffs[256]);
#endif
}

/*************************************************
* Name:        poly_tobytes
*
* Description: Serialization of a polynomial
*
* Arguments:   - unsigned char *r: pointer to output byte array
*              - const poly *a:    pointer to input polynomial
**************************************************/
void poly_tobytes(unsigned char *r, const poly *a) {
  poly_tobytes_avx2(r, &a->coeffs[0]);
#if (OKAI_N == 512)
  poly_tobytes_avx2(r + OKAI_POLYBYTES/2, &a->coeffs[256]);
#endif
}

/*************************************************
* Name:        poly_frombytes
*
* Description: De-serialization of a polynomial;
*              inverse of poly_tobytes
*
* Arguments:   - poly *r:                pointer to output polynomial
*              - const unsigned char *a: pointer to input byte array
**************************************************/
void poly_frombytes(poly *r, const unsigned char *a) 
{
  poly_frombytes_avx2(&r->coeffs[0], a);
#if (OKAI_N == 512)
  poly_frombytes_avx2(&r->coeffs[256], a + OKAI_POLYBYTES/2);
#endif
}

/*************************************************
* Name:        poly_frommsg
*
* Description: Convert 32-byte message to polynomial
*
* Arguments:   - poly *r:            pointer to output polynomial
*              - const uint8_t *msg: pointer to input message
**************************************************/
void poly_frommsg(poly * restrict r,
                  const uint8_t msg[OKAI_INDCPA_MSGBYTES])
{
  __m256i f, g0, g1, g2, g3, h0, h1, h2, h3;
  const __m256i shift = _mm256_broadcastsi128_si256(_mm_set_epi32(0,1,2,3));
  const __m256i idx = _mm256_broadcastsi128_si256(_mm_set_epi8(15,14,11,10,7,6,3,2,13,12,9,8,5,4,1,0));
  const __m256i hqs = _mm256_set1_epi16((OKAI_Q+1)/2);

#define FROMMSG64(i)						\
  g3 = _mm256_shuffle_epi32(f,0x55*i);				\
  g3 = _mm256_sllv_epi32(g3,shift);				\
  g3 = _mm256_shuffle_epi8(g3,idx);				\
  g0 = _mm256_slli_epi16(g3,12);				\
  g1 = _mm256_slli_epi16(g3,8);					\
  g2 = _mm256_slli_epi16(g3,4);					\
  g0 = _mm256_srai_epi16(g0,15);				\
  g1 = _mm256_srai_epi16(g1,15);				\
  g2 = _mm256_srai_epi16(g2,15);				\
  g3 = _mm256_srai_epi16(g3,15);				\
  g0 = _mm256_and_si256(g0,hqs);  /* 19 18 17 16  3  2  1  0 */	\
  g1 = _mm256_and_si256(g1,hqs);  /* 23 22 21 20  7  6  5  4 */	\
  g2 = _mm256_and_si256(g2,hqs);  /* 27 26 25 24 11 10  9  8 */	\
  g3 = _mm256_and_si256(g3,hqs);  /* 31 30 29 28 15 14 13 12 */	\
  h0 = _mm256_unpacklo_epi64(g0,g1);				\
  h2 = _mm256_unpackhi_epi64(g0,g1);				\
  h1 = _mm256_unpacklo_epi64(g2,g3);				\
  h3 = _mm256_unpackhi_epi64(g2,g3);				\
  g0 = _mm256_permute2x128_si256(h0,h1,0x20);			\
  g2 = _mm256_permute2x128_si256(h0,h1,0x31);			\
  g1 = _mm256_permute2x128_si256(h2,h3,0x20);			\
  g3 = _mm256_permute2x128_si256(h2,h3,0x31);			\
  _mm256_store_si256((__m256i *)&r->coeffs[  0+32*i+ 0],g0);	\
  _mm256_store_si256((__m256i *)&r->coeffs[  0+32*i+16],g1);	\
  _mm256_store_si256((__m256i *)&r->coeffs[128+32*i+ 0],g2);	\
  _mm256_store_si256((__m256i *)&r->coeffs[128+32*i+16],g3)

  f = _mm256_load_si256((__m256i *)msg);
  FROMMSG64(0);
  FROMMSG64(1);
  FROMMSG64(2);
  FROMMSG64(3);
}

/*************************************************
* Name:        poly_tomsg
*
* Description: Convert polynomial to 32-byte message
*
* Arguments:   - uint8_t *msg: pointer to output message
*              - poly *a:      pointer to input polynomial
**************************************************/
void poly_tomsg(uint8_t msg[OKAI_INDCPA_MSGBYTES], poly * restrict a)
{
  unsigned int i;
  uint32_t small;
  __m256i f0, f1, g0, g1;
  const __m256i hqs = _mm256_set1_epi16((OKAI_Q - 1)/2);
  const __m256i hhqs = _mm256_set1_epi16((OKAI_Q - 5)/4);

  for(i=0;i<OKAI_N/32;i++) {
    f0 = _mm256_load_si256((__m256i *)&a->coeffs[32*i]);
    f1 = _mm256_load_si256((__m256i *)&a->coeffs[32*i+16]);
    f0 = _mm256_sub_epi16(hqs, f0);
    f1 = _mm256_sub_epi16(hqs, f1);
    g0 = _mm256_srai_epi16(f0, 15);
    g1 = _mm256_srai_epi16(f1, 15);
    f0 = _mm256_xor_si256(f0, g0);
    f1 = _mm256_xor_si256(f1, g1);
    f0 = _mm256_sub_epi16(hhqs, f0);
    f1 = _mm256_sub_epi16(hhqs, f1);
    f0 = _mm256_packs_epi16(f0, f1);
    small = _mm256_movemask_epi8(f0);
    small = ~small;
    msg[4*i+0] = small;
    msg[4*i+1] = small >> 16;
    msg[4*i+2] = small >>  8;
    msg[4*i+3] = small >> 24;
  }
}

/*************************************************
* Name:        poly_getnoise
*
* Description: Sample a polynomial deterministically from a seed and a nonce,
*              with output polynomial close to centered binomial distribution
*              with parameter OKAI_ETA
*
* Arguments:   - poly *r:                   pointer to output polynomial
*              - const unsigned char *seed: pointer to input seed
*              - unsigned char nonce:       one-byte input nonce
**************************************************/
void poly_getnoises(poly *r, const uint8_t seed[OKAI_SYMBYTES], uint8_t nonce)
{
  __attribute__((aligned(32)))
  uint8_t buf[OKAI_SETA*OKAI_N/4];
  prf(buf, sizeof(buf), seed, nonce);
  cbd1(r, buf);
}

void poly_getnoisee(poly *r, const uint8_t seed[OKAI_SYMBYTES], uint8_t nonce)
{
  __attribute__((aligned(32)))
  uint8_t buf[OKAI_EETA*OKAI_N/4];
  prf(buf, sizeof(buf), seed, nonce);
#if (OKAI_N == 256)
  cbd4(r, buf);
#elif (OKAI_N == 512)
  cbd6(r, buf);
#endif
}

/*************************************************
* Name:        poly_ntt
*
* Description: Computes negacyclic number-theoretic transform (NTT) of
*              a polynomial in place;
*              inputs assumed to be in normal order, output in bitreversed order
*
* Arguments:   - uint16_t *r: pointer to in/output polynomial
**************************************************/
void poly_ntt(poly *r)
{
  ntt_level0_avx(r->coeffs, zetas_exp);
  ntt_level0_avx(r->coeffs + 64, zetas_exp);
  ntt_levels1t6_avx(r->coeffs, zetas_exp + 4);
  ntt_levels1t6_avx(r->coeffs + 128, zetas_exp + 200);
#if (OKAI_N == 512)
  ntt_level0_avx(r->coeffs + 256, zetas_exp);
  ntt_level0_avx(r->coeffs + 320, zetas_exp);
  ntt_levels1t6_avx(r->coeffs + 256, zetas_exp + 4);
  ntt_levels1t6_avx(r->coeffs + 384, zetas_exp + 200);
#endif
}

/*************************************************
* Name:        poly_invntt
*
* Description: Computes inverse of negacyclic number-theoretic transform (NTT) of
*              a polynomial in place;
*              inputs assumed to be in bitreversed order, output in normal order
*
* Arguments:   - uint16_t *a: pointer to in/output polynomial
**************************************************/
void poly_invntt(poly *r)
{
  invntt_levels0t5_avx(r->coeffs, zetas_inv_exp);
  invntt_levels0t5_avx(r->coeffs + 128, zetas_inv_exp + 196);
  invntt_level6_avx(r->coeffs, zetas_inv_exp + 392);
#if (OKAI_N == 512)
  invntt_levels0t5_avx(r->coeffs + 256, zetas_inv_exp);
  invntt_levels0t5_avx(r->coeffs + 384, zetas_inv_exp + 196);
  invntt_level6_avx(r->coeffs + 256, zetas_inv_exp + 392);
#endif
}

void poly_nttpack(poly *r)
{
  nttpack_avx(r->coeffs);
  nttpack_avx(r->coeffs + 128);
#if (OKAI_N == 512)
  nttpack_avx(r->coeffs + 256);
  nttpack_avx(r->coeffs + 384);
#endif
}

void poly_nttunpack(poly *r)
{
  nttunpack_avx(r->coeffs);
  nttunpack_avx(r->coeffs + 128);
#if (OKAI_N == 512)
  nttunpack_avx(r->coeffs + 256);
  nttunpack_avx(r->coeffs + 384);
#endif
}

static void poly_add_256(int16_t * restrict r,
              const int16_t * restrict a,
              const int16_t * restrict b)
{
    unsigned int i;
    __m256i f0, f1;

    for(i=0;i<256;i+=16) {
        f0 = _mm256_load_si256((__m256i *)&a[i]);
        f1 = _mm256_load_si256((__m256i *)&b[i]);
        f0 = _mm256_add_epi16(f0, f1);
        _mm256_store_si256((__m256i *)&r[i], f0);
    }
}

static void poly_sub_256(int16_t * restrict r,
              const int16_t * restrict a,
              const int16_t * restrict b)
{
    unsigned int i;
    __m256i f0, f1;

    for(i=0;i<256;i+=16) {
        f0 = _mm256_load_si256((__m256i *)&a[i]);
        f1 = _mm256_load_si256((__m256i *)&b[i]);
        f0 = _mm256_sub_epi16(f0, f1);
        _mm256_store_si256((__m256i *)&r[i], f0);
    }
}

static void basemul_multi_avx(int16_t *a, int16_t *b, int16_t *r, __m256i *vec);
void poly_basemul_montgomery_multi(poly *r, const poly *a, const poly *b)
{
  basemul_multi_avx(a->coeffs,b->coeffs,r->coeffs,qdata.vec);
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
  reduce_avx(r->coeffs);
  reduce_avx(r->coeffs + 128);
#if (OKAI_N == 512)
  reduce_avx(r->coeffs + 256);
  reduce_avx(r->coeffs + 384);
#endif
}

/*************************************************
* Name:        poly_add
*
* Description: Add two polynomials
*
* Arguments: - poly *r:       pointer to output polynomial
*            - const poly *a: pointer to first input polynomial
*            - const poly *b: pointer to second input polynomial
**************************************************/
void poly_add(poly * restrict r,
              const poly * restrict a,
              const poly * restrict b)
{
  unsigned int i;
  __m256i f0, f1;

  for(i=0;i<OKAI_N;i+=16) {
    f0 = _mm256_load_si256((__m256i *)&a->coeffs[i]);
    f1 = _mm256_load_si256((__m256i *)&b->coeffs[i]);
    f0 = _mm256_add_epi16(f0, f1);
    _mm256_store_si256((__m256i *)&r->coeffs[i], f0);
  }
}

/*************************************************
* Name:        poly_sub
*
* Description: Subtract two polynomials
*
* Arguments: - poly *r:       pointer to output polynomial
*            - const poly *a: pointer to first input polynomial
*            - const poly *b: pointer to second input polynomial
**************************************************/
void poly_sub(poly * restrict r,
              const poly * restrict a,
              const poly * restrict b)
{
  unsigned int i;
  __m256i f0, f1;

  for(i=0;i<OKAI_N;i+=16) {
    f0 = _mm256_load_si256((__m256i *)&a->coeffs[i]);
    f1 = _mm256_load_si256((__m256i *)&b->coeffs[i]);
    f0 = _mm256_sub_epi16(f0, f1);
    _mm256_store_si256((__m256i *)&r->coeffs[i], f0);
  }
}
