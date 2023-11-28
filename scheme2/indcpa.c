#include <string.h>
#include <stdio.h>
#include "indcpa.h"
#include "poly.h"
#include "polyvec.h"
#include "randombytes.h"
#include "ntt.h"
#include "symmetric.h"
#include "rejsample.h"
#include "cbd.h"
#include "consts.h"

#define gen_a(A,B)  gen_matrix(A,B,0)
#define gen_at(A,B) gen_matrix(A,B,1)

void pack_sk(unsigned char *r, polyvec *sk) {
  int i;
  for(i=0; i<OKAI_K; i++) {
    poly_nttpack(&sk->vec[i]);
  }
  polyvec_tobytes(r, sk);
}

void unpack_sk(polyvec *sk, const unsigned char *packedsk) {
  int i;
  polyvec_frombytes(sk, packedsk);
  for(i=0; i<OKAI_K; i++)
    poly_nttunpack(&sk->vec[i]);
}

void pack_pk(unsigned char *r, const polyvec *pk, const unsigned char *seed) {
  int i;
  polyvec_compress(r, pk);
  for(i=0;i<OKAI_SYMBYTES;i++)
    r[i+OKAI_POLYVECCOMPRESSEDPKBYTES] = seed[i];
}

void unpack_pk(polyvec *pk, unsigned char *seed, const unsigned char *packedpk) {
  int i;
  polyvec_decompress(pk, packedpk);
  memcpy(seed, packedpk+OKAI_POLYVECCOMPRESSEDPKBYTES, OKAI_SYMBYTES);
}

static unsigned int rej_uniform(int16_t *r, unsigned int len, const unsigned char *buf, unsigned int buflen) {
  unsigned int ctr, pos;
  uint16_t val;

  ctr = pos = 0;
  while(ctr < len && pos + 2 <= buflen)
  {
    val = buf[pos] | ((uint16_t)buf[pos+1] << 8);
    pos += 2;

    if(val < 8*OKAI_Q)
    {
      val -= ((int32_t)(val*8737) >> 26) * OKAI_Q;
      r[ctr++] = (int16_t)val;
    }
  }

  return ctr;
}

/*************************************************
* Name:        gen_matrix
*
* Description: Deterministically generate matrix A (or the transpose of A)
*              from a seed. Entries of the matrix are polynomials that look
*              uniformly random. Performs rejection sampling on output of
*              a XOF
*
* Arguments:   - polyvec *a:                pointer to ouptput matrix A
*              - const unsigned char *seed: pointer to input seed
*              - int transposed:            boolean deciding whether A or A^T is generated
**************************************************/
#define GEN_MATRIX_NBLOCKS ((2*OKAI_N*(1U << 16)/(8*OKAI_Q) \
                             + XOF_BLOCKBYTES)/XOF_BLOCKBYTES)
#if OKAI_K == 2
void gen_matrix(polyvec *a, const uint8_t seed[32], int transposed)
{
  unsigned int ctr0, ctr1, ctr2, ctr3;
  __attribute__((aligned(32)))
  uint8_t buf[4][(GEN_MATRIX_NBLOCKS*XOF_BLOCKBYTES+31)/32*32];
  __m256i f;
  keccakx4_state state;

  f = _mm256_load_si256((__m256i *)seed);
  _mm256_store_si256((__m256i *)buf[0], f);
  _mm256_store_si256((__m256i *)buf[1], f);
  _mm256_store_si256((__m256i *)buf[2], f);
  _mm256_store_si256((__m256i *)buf[3], f);
#if (OKAI_N == 512)
  f = _mm256_load_si256((__m256i *)(&seed[32]));
  _mm256_store_si256((__m256i *)(&buf[0][32]), f);
  _mm256_store_si256((__m256i *)(&buf[1][32]), f);
  _mm256_store_si256((__m256i *)(&buf[2][32]), f);
  _mm256_store_si256((__m256i *)(&buf[3][32]), f);
#endif
  if(transposed) {
    buf[0][OKAI_SYMBYTES+0] = 0;
    buf[0][OKAI_SYMBYTES+1] = 0;
    buf[1][OKAI_SYMBYTES+0] = 0;
    buf[1][OKAI_SYMBYTES+1] = 1;
    buf[2][OKAI_SYMBYTES+0] = 1;
    buf[2][OKAI_SYMBYTES+1] = 0;
    buf[3][OKAI_SYMBYTES+0] = 1;
    buf[3][OKAI_SYMBYTES+1] = 1;
  }
  else {
    buf[0][OKAI_SYMBYTES+0] = 0;
    buf[0][OKAI_SYMBYTES+1] = 0;
    buf[1][OKAI_SYMBYTES+0] = 1;
    buf[1][OKAI_SYMBYTES+1] = 0;
    buf[2][OKAI_SYMBYTES+0] = 0;
    buf[2][OKAI_SYMBYTES+1] = 1;
    buf[3][OKAI_SYMBYTES+0] = 1;
    buf[3][OKAI_SYMBYTES+1] = 1;
  }

  shake128x4_absorb(&state, buf[0], buf[1], buf[2], buf[3], OKAI_SYMBYTES+2);
  shake128x4_squeezeblocks(buf[0], buf[1], buf[2], buf[3], GEN_MATRIX_NBLOCKS,
                           &state);

  ctr0 = rej_uniform_avx(a[0].vec[0].coeffs, buf[0]);
  ctr1 = rej_uniform_avx(a[0].vec[1].coeffs, buf[1]);
  ctr2 = rej_uniform_avx(a[1].vec[0].coeffs, buf[2]);
  ctr3 = rej_uniform_avx(a[1].vec[1].coeffs, buf[3]);

  while(ctr0 < OKAI_N || ctr1 < OKAI_N || ctr2 < OKAI_N || ctr3 < OKAI_N) {
    shake128x4_squeezeblocks(buf[0], buf[1], buf[2], buf[3], 1, &state);

    ctr0 += rej_uniform(a[0].vec[0].coeffs + ctr0, OKAI_N - ctr0, buf[0],
                        XOF_BLOCKBYTES);
    ctr1 += rej_uniform(a[0].vec[1].coeffs + ctr1, OKAI_N - ctr1, buf[1],
                        XOF_BLOCKBYTES);
    ctr2 += rej_uniform(a[1].vec[0].coeffs + ctr2, OKAI_N - ctr2, buf[2],
                        XOF_BLOCKBYTES);
    ctr3 += rej_uniform(a[1].vec[1].coeffs + ctr3, OKAI_N - ctr3, buf[3],
                        XOF_BLOCKBYTES);
  }


  poly_nttunpack(&a[0].vec[0]);
  poly_nttunpack(&a[0].vec[1]);
  poly_nttunpack(&a[1].vec[0]);
  poly_nttunpack(&a[1].vec[1]);
}
#elif OKAI_K == 3
void gen_matrix(polyvec *a, const uint8_t seed[32], int transposed)
{
  unsigned int ctr0, ctr1, ctr2, ctr3;
  __attribute__((aligned(32)))
  uint8_t buf[4][(GEN_MATRIX_NBLOCKS*XOF_BLOCKBYTES+31)/32*32];
  __m256i f;
  keccakx4_state state;
  keccak_state state1x;

  f = _mm256_load_si256((__m256i *)seed);
  _mm256_store_si256((__m256i *)buf[0], f);
  _mm256_store_si256((__m256i *)buf[1], f);
  _mm256_store_si256((__m256i *)buf[2], f);
  _mm256_store_si256((__m256i *)buf[3], f);

  if(transposed) {
    buf[0][OKAI_SYMBYTES+0] = 0;
    buf[0][OKAI_SYMBYTES+1] = 0;
    buf[1][OKAI_SYMBYTES+0] = 0;
    buf[1][OKAI_SYMBYTES+1] = 1;
    buf[2][OKAI_SYMBYTES+0] = 0;
    buf[2][OKAI_SYMBYTES+1] = 2;
    buf[3][OKAI_SYMBYTES+0] = 1;
    buf[3][OKAI_SYMBYTES+1] = 0;
  }
  else {
    buf[0][OKAI_SYMBYTES+0] = 0;
    buf[0][OKAI_SYMBYTES+1] = 0;
    buf[1][OKAI_SYMBYTES+0] = 1;
    buf[1][OKAI_SYMBYTES+1] = 0;
    buf[2][OKAI_SYMBYTES+0] = 2;
    buf[2][OKAI_SYMBYTES+1] = 0;
    buf[3][OKAI_SYMBYTES+0] = 0;
    buf[3][OKAI_SYMBYTES+1] = 1;
  }

  shake128x4_absorb(&state, buf[0], buf[1], buf[2], buf[3], OKAI_SYMBYTES+2);
  shake128x4_squeezeblocks(buf[0], buf[1], buf[2], buf[3], GEN_MATRIX_NBLOCKS,
                           &state);

  ctr0 = rej_uniform_avx(a[0].vec[0].coeffs, buf[0]);
  ctr1 = rej_uniform_avx(a[0].vec[1].coeffs, buf[1]);
  ctr2 = rej_uniform_avx(a[0].vec[2].coeffs, buf[2]);
  ctr3 = rej_uniform_avx(a[1].vec[0].coeffs, buf[3]);

  while(ctr0 < OKAI_N || ctr1 < OKAI_N || ctr2 < OKAI_N || ctr3 < OKAI_N) {
    shake128x4_squeezeblocks(buf[0], buf[1], buf[2], buf[3], 1, &state);

    ctr0 += rej_uniform(a[0].vec[0].coeffs + ctr0, OKAI_N - ctr0, buf[0],
                        XOF_BLOCKBYTES);
    ctr1 += rej_uniform(a[0].vec[1].coeffs + ctr1, OKAI_N - ctr1, buf[1],
                        XOF_BLOCKBYTES);
    ctr2 += rej_uniform(a[0].vec[2].coeffs + ctr2, OKAI_N - ctr2, buf[2],
                        XOF_BLOCKBYTES);
    ctr3 += rej_uniform(a[1].vec[0].coeffs + ctr3, OKAI_N - ctr3, buf[3],
                        XOF_BLOCKBYTES);
  }

  poly_nttunpack(&a[0].vec[0]);
  poly_nttunpack(&a[0].vec[1]);
  poly_nttunpack(&a[0].vec[2]);
  poly_nttunpack(&a[1].vec[0]);

  f = _mm256_load_si256((__m256i *)seed);
  _mm256_store_si256((__m256i *)buf[0], f);
  _mm256_store_si256((__m256i *)buf[1], f);
  _mm256_store_si256((__m256i *)buf[2], f);
  _mm256_store_si256((__m256i *)buf[3], f);

  if(transposed) {
    buf[0][OKAI_SYMBYTES+0] = 1;
    buf[0][OKAI_SYMBYTES+1] = 1;
    buf[1][OKAI_SYMBYTES+0] = 1;
    buf[1][OKAI_SYMBYTES+1] = 2;
    buf[2][OKAI_SYMBYTES+0] = 2;
    buf[2][OKAI_SYMBYTES+1] = 0;
    buf[3][OKAI_SYMBYTES+0] = 2;
    buf[3][OKAI_SYMBYTES+1] = 1;
  }
  else {
    buf[0][OKAI_SYMBYTES+0] = 1;
    buf[0][OKAI_SYMBYTES+1] = 1;
    buf[1][OKAI_SYMBYTES+0] = 2;
    buf[1][OKAI_SYMBYTES+1] = 1;
    buf[2][OKAI_SYMBYTES+0] = 0;
    buf[2][OKAI_SYMBYTES+1] = 2;
    buf[3][OKAI_SYMBYTES+0] = 1;
    buf[3][OKAI_SYMBYTES+1] = 2;
  }

  shake128x4_absorb(&state, buf[0], buf[1], buf[2], buf[3], OKAI_SYMBYTES+2);
  shake128x4_squeezeblocks(buf[0], buf[1], buf[2], buf[3], GEN_MATRIX_NBLOCKS,
                           &state);

  ctr0 = rej_uniform_avx(a[1].vec[1].coeffs, buf[0]);
  ctr1 = rej_uniform_avx(a[1].vec[2].coeffs, buf[1]);
  ctr2 = rej_uniform_avx(a[2].vec[0].coeffs, buf[2]);
  ctr3 = rej_uniform_avx(a[2].vec[1].coeffs, buf[3]);

  while(ctr0 < OKAI_N || ctr1 < OKAI_N || ctr2 < OKAI_N || ctr3 < OKAI_N) {
    shake128x4_squeezeblocks(buf[0], buf[1], buf[2], buf[3], 1, &state);

    ctr0 += rej_uniform(a[1].vec[1].coeffs + ctr0, OKAI_N - ctr0, buf[0],
                        XOF_BLOCKBYTES);
    ctr1 += rej_uniform(a[1].vec[2].coeffs + ctr1, OKAI_N - ctr1, buf[1],
                        XOF_BLOCKBYTES);
    ctr2 += rej_uniform(a[2].vec[0].coeffs + ctr2, OKAI_N - ctr2, buf[2],
                        XOF_BLOCKBYTES);
    ctr3 += rej_uniform(a[2].vec[1].coeffs + ctr3, OKAI_N - ctr3, buf[3],
                        XOF_BLOCKBYTES);
  }

  poly_nttunpack(&a[1].vec[1]);
  poly_nttunpack(&a[1].vec[2]);
  poly_nttunpack(&a[2].vec[0]);
  poly_nttunpack(&a[2].vec[1]);

  f = _mm256_load_si256((__m256i *)seed);
  _mm256_store_si256((__m256i *)buf[0], f);
  buf[0][OKAI_SYMBYTES+0] = 2;
  buf[0][OKAI_SYMBYTES+1] = 2;
  shake128_absorb(&state1x, buf[0], OKAI_SYMBYTES+2);
  shake128_squeezeblocks(buf[0], GEN_MATRIX_NBLOCKS, &state1x);
  ctr0 = rej_uniform_avx(a[2].vec[2].coeffs, buf[0]);
  while(ctr0 < OKAI_N)
  {
    shake128_squeezeblocks(buf[0], 1, &state1x);
    ctr0 += rej_uniform(a[2].vec[2].coeffs + ctr0, OKAI_N - ctr0, buf[0],
                        XOF_BLOCKBYTES);
  }

  poly_nttunpack(&a[2].vec[2]);
}
#endif

void indcpa_keypair(unsigned char *pk,
                   unsigned char *sk)
{
  unsigned int i;
  __attribute__((aligned(32)))
  uint8_t buf[OKAI_SYMBYTES+OKAI_SYMBYTES];
  const uint8_t *publicseed = buf;
  const uint8_t *noiseseed = buf+OKAI_SYMBYTES;
  polyvec a[OKAI_K], e, pkpv, skpv;

  randombytes(buf, OKAI_SYMBYTES);
  hash_g(buf, buf, OKAI_SYMBYTES);

  gen_matrix(a, publicseed, 0);

#if OKAI_K == 2
  poly_getnoises(skpv.vec, noiseseed, 0);
  poly_getnoises(skpv.vec + 1, noiseseed, 1);
  poly_getnoisee(e.vec, noiseseed, 2);
  poly_getnoisee(e.vec + 1, noiseseed, 3);
#elif OKAI_K == 3
  poly_getnoises(skpv.vec, noiseseed, 0);
  poly_getnoises(skpv.vec + 1, noiseseed, 1);
  poly_getnoises(skpv.vec + 2, noiseseed, 2);
  poly_getnoisee(e.vec, noiseseed, 3);
  poly_getnoisee(e.vec + 1, noiseseed, 4);
  poly_getnoisee(e.vec + 2, noiseseed, 5);
#endif

  polyvec_ntt(&skpv);

  // matrix-vector multiplication
  for(i=0;i<OKAI_K;i++)
    polyvec_pointwise_acc(&pkpv.vec[i],&skpv,a+i);
  
  polyvec_invntt(&pkpv);
  polyvec_add(&pkpv,&pkpv,&e);
  polyvec_reduce(&pkpv);

  pack_sk(sk, &skpv);
  pack_pk(pk, &pkpv, publicseed);
}


void indcpa_enc(unsigned char *c,
               const unsigned char *m,
               const unsigned char *pk,
               const unsigned char *coins)
{
  unsigned int i;
  __attribute__((aligned(32)))
  uint8_t seed[OKAI_SYMBYTES];
  polyvec sp, pkpv, ep, at[OKAI_K], bp;
  poly v, epp;

  unpack_pk(&pkpv, seed, pk);

  polyvec_ntt(&pkpv);
  gen_matrix(at, seed, 1);

#if OKAI_K == 2
  poly_getnoises(sp.vec, coins, 0);
  poly_getnoises(sp.vec + 1, coins, 1);
  poly_getnoisee(ep.vec, coins, 2);
  poly_getnoisee(ep.vec + 1, coins, 3);
  poly_getnoisee(&epp, coins, 4);
#elif OKAI_K == 3
  poly_getnoises(sp.vec, coins, 0);
  poly_getnoises(sp.vec + 1, coins, 1);
  poly_getnoises(sp.vec + 2, coins, 2);
  poly_getnoisee(ep.vec, coins, 3);
  poly_getnoisee(ep.vec + 1, coins, 4);
  poly_getnoisee(ep.vec + 2, coins, 5);
  poly_getnoisee(&epp, coins, 6);
#endif

  polyvec_ntt(&sp);

  // matrix-vector multiplication
  for(i=0;i<OKAI_K;i++)
    polyvec_pointwise_acc(&bp.vec[i],&sp,at+i);
  polyvec_invntt(&bp);


  polyvec_pointwise_acc(&v, &pkpv, &sp);
  poly_invntt(&v);


  polyvec_add(&bp, &bp, &ep);
  polyvec_reduce(&bp);

  poly_add(&v, &v, &epp);
  
  poly_reduce(&v);
  polyvec_compressc(c, &bp);

  poly_con(c+OKAI_POLYVECCOMPRESSEDCBYTES, m, &v);
}


void indcpa_dec(unsigned char *m,
               const unsigned char *c,
               const unsigned char *sk)
{
  polyvec bp, skpv;
  poly mp;

  polyvec_decompressc(&bp, c);
  
  unpack_sk(&skpv, sk);

  polyvec_ntt(&bp);

  polyvec_pointwise_acc(&mp,&skpv,&bp);
  poly_invntt(&mp);
  poly_reduce(&mp);

  poly_rec(m, c+OKAI_POLYVECCOMPRESSEDCBYTES, &mp);
}