#include <stddef.h>
#include <stdint.h>
#include <immintrin.h>
#include <string.h>
#include "align.h"
#include "params.h"
#include "indcpa.h"
#include "polyvec.h"
#include "poly.h"
#include "ntt.h"
#include "cbd.h"
#include "rejsample.h"
#include "symmetric.h"
#include "randombytes.h"

/*************************************************
* Name:        pack_pk
*
* Description: Serialize the public key as concatenation of the
*              serialized vector of polynomials pk and the
*              public seed used to generate the matrix A.
*              The polynomial coefficients in pk are assumed to
*              lie in the invertal [0,q], i.e. pk must be reduced
*              by polyvec_reduce().
*
* Arguments:   uint8_t *r: pointer to the output serialized public key
*              polyvec *pk: pointer to the input public-key polyvec
*              const uint8_t *seed: pointer to the input public seed
**************************************************/
static void pack_pk(uint8_t r[OSKR_INDCPA_PUBLICKEYBYTES],
                    polyvec *pk,
                    const uint8_t seed[OSKR_SYMBYTES])
{
  polyvec_tobytes(r, pk);
  memcpy(r+OSKR_POLYVECBYTES, seed, OSKR_SYMBYTES);
}

/*************************************************
* Name:        unpack_pk
*
* Description: De-serialize public key from a byte array;
*              approximate inverse of pack_pk
*
* Arguments:   - polyvec *pk: pointer to output public-key polynomial vector
*              - uint8_t *seed: pointer to output seed to generate matrix A
*              - const uint8_t *packedpk: pointer to input serialized public key
**************************************************/
void unpack_pk(polyvec *pk,
                      uint8_t seed[OSKR_SYMBYTES],
                      const uint8_t packedpk[OSKR_INDCPA_PUBLICKEYBYTES])
{
  polyvec_frombytes(pk, packedpk);
  memcpy(seed, packedpk+OSKR_POLYVECBYTES, OSKR_SYMBYTES);
}

/*************************************************
* Name:        pack_sk
*
* Description: Serialize the secret key.
*              The polynomial coefficients in sk are assumed to
*              lie in the invertal [0,q], i.e. sk must be reduced
*              by polyvec_reduce().
*
* Arguments:   - uint8_t *r: pointer to output serialized secret key
*              - polyvec *sk: pointer to input vector of polynomials (secret key)
**************************************************/
static void pack_sk(uint8_t r[OSKR_INDCPA_SECRETKEYBYTES], polyvec *sk)
{
  polyvec_tobytes(r, sk);
}

/*************************************************
* Name:        unpack_sk
*
* Description: De-serialize the secret key; inverse of pack_sk
*
* Arguments:   - polyvec *sk: pointer to output vector of polynomials (secret key)
*              - const uint8_t *packedsk: pointer to input serialized secret key
**************************************************/
void unpack_sk(polyvec *sk, const uint8_t packedsk[OSKR_INDCPA_SECRETKEYBYTES])
{
  polyvec_frombytes(sk, packedsk);
}

/*************************************************
* Name:        rej_uniform
*
* Description: Run rejection sampling on uniform random bytes to generate
*              uniform random integers mod q
*
* Arguments:   - int16_t *r: pointer to output array
*              - unsigned int len: requested number of 16-bit integers (uniform mod q)
*              - const uint8_t *buf: pointer to input buffer (assumed to be uniformly random bytes)
*              - unsigned int buflen: length of input buffer in bytes
*
* Returns number of sampled 16-bit integers (at most len)
**************************************************/
static unsigned int rej_uniform(int16_t *r,
                                unsigned int len,
                                const uint8_t *buf,
                                unsigned int buflen)
{
  unsigned int ctr, pos;
  uint16_t val0, val1;

  ctr = pos = 0;
  while(ctr < len && pos <= buflen - 3) {
    val0 = ((buf[pos+0] >> 0) | ((uint16_t)buf[pos+1] << 8)) & 0xFFF;
    val1 = ((buf[pos+1] >> 4) | ((uint16_t)buf[pos+2] << 4)) & 0xFFF;
    pos += 3;

    if(val0 < OSKR_Q)
      r[ctr++] = val0;
    if(ctr < len && val1 < OSKR_Q)
      r[ctr++] = val1;
  }

  return ctr;
}

#define gen_a(A,B)  gen_matrix(A,B,0)
#define gen_at(A,B) gen_matrix(A,B,1)

/*************************************************
* Name:        gen_matrix
*
* Description: Deterministically generate matrix A (or the transpose of A)
*              from a seed. Entries of the matrix are polynomials that look
*              uniformly random. Performs rejection sampling on output of
*              a XOF
*
* Arguments:   - polyvec *a: pointer to ouptput matrix A
*              - const uint8_t *seed: pointer to input seed
*              - int transposed: boolean deciding whether A or A^T is generated
**************************************************/
#define GEN_MATRIX_NBLOCKS (REJ_UNIFORM_AVX_BUFLEN/XOF_BLOCKBYTES)

#if OSKR_K == 2
void gen_matrix(polyvec *a, const uint8_t seed[OSKR_SYMBYTES], int transposed)
{
  unsigned int ctr0, ctr1, ctr2, ctr3;
  ALIGNED_UINT8((GEN_MATRIX_NBLOCKS*XOF_BLOCKBYTES+31)/32*32) buf[4];
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
  
  if(transposed) {
    buf[0].coeffs[OSKR_SYMBYTES+0] = 0;
    buf[0].coeffs[OSKR_SYMBYTES+1] = 0;
    buf[1].coeffs[OSKR_SYMBYTES+0] = 0;
    buf[1].coeffs[OSKR_SYMBYTES+1] = 1;
    buf[2].coeffs[OSKR_SYMBYTES+0] = 1;
    buf[2].coeffs[OSKR_SYMBYTES+1] = 0;
    buf[3].coeffs[OSKR_SYMBYTES+0] = 1;
    buf[3].coeffs[OSKR_SYMBYTES+1] = 1;
  }
  else {
    buf[0].coeffs[OSKR_SYMBYTES+0] = 0;
    buf[0].coeffs[OSKR_SYMBYTES+1] = 0;
    buf[1].coeffs[OSKR_SYMBYTES+0] = 1;
    buf[1].coeffs[OSKR_SYMBYTES+1] = 0;
    buf[2].coeffs[OSKR_SYMBYTES+0] = 0;
    buf[2].coeffs[OSKR_SYMBYTES+1] = 1;
    buf[3].coeffs[OSKR_SYMBYTES+0] = 1;
    buf[3].coeffs[OSKR_SYMBYTES+1] = 1;
  }
  
  shake128x4_absorb_once(&state, buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, OSKR_SYMBYTES+2);
  shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, GEN_MATRIX_NBLOCKS, &state);

  ctr0 = rej_uniform_avx(a[0].vec[0].coeffs, buf[0].coeffs);
  ctr1 = rej_uniform_avx(a[0].vec[1].coeffs, buf[1].coeffs);
  ctr2 = rej_uniform_avx(a[1].vec[0].coeffs, buf[2].coeffs);
  ctr3 = rej_uniform_avx(a[1].vec[1].coeffs, buf[3].coeffs);

  while(ctr0 < OSKR_N || ctr1 < OSKR_N || ctr2 < OSKR_N || ctr3 < OSKR_N) {
    shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 1, &state);

    ctr0 += rej_uniform(a[0].vec[0].coeffs + ctr0, OSKR_N - ctr0, buf[0].coeffs, SHAKE128_RATE);
    ctr1 += rej_uniform(a[0].vec[1].coeffs + ctr1, OSKR_N - ctr1, buf[1].coeffs, SHAKE128_RATE);
    ctr2 += rej_uniform(a[1].vec[0].coeffs + ctr2, OSKR_N - ctr2, buf[2].coeffs, SHAKE128_RATE);
    ctr3 += rej_uniform(a[1].vec[1].coeffs + ctr3, OSKR_N - ctr3, buf[3].coeffs, SHAKE128_RATE);
  }

  poly_nttunpack(&a[0].vec[0]);
  poly_nttunpack(&a[0].vec[1]);
  poly_nttunpack(&a[1].vec[0]);
  poly_nttunpack(&a[1].vec[1]);
}
#elif OSKR_K == 3
void gen_matrix(polyvec *a, const uint8_t seed[32], int transposed)
{
  unsigned int ctr0, ctr1, ctr2, ctr3;
  ALIGNED_UINT8(REJ_UNIFORM_AVX_NBLOCKS*SHAKE128_RATE) buf[4];
  __m256i f;
  keccakx4_state state;
  keccak_state state1x;

  f = _mm256_loadu_si256((__m256i *)seed);
  _mm256_store_si256(buf[0].vec, f);
  _mm256_store_si256(buf[1].vec, f);
  _mm256_store_si256(buf[2].vec, f);
  _mm256_store_si256(buf[3].vec, f);

  if(transposed) {
    buf[0].coeffs[OSKR_SYMBYTES+0] = 0;
    buf[0].coeffs[OSKR_SYMBYTES+1] = 0;
    buf[1].coeffs[OSKR_SYMBYTES+0] = 0;
    buf[1].coeffs[OSKR_SYMBYTES+1] = 1;
    buf[2].coeffs[OSKR_SYMBYTES+0] = 0;
    buf[2].coeffs[OSKR_SYMBYTES+1] = 2;
    buf[3].coeffs[OSKR_SYMBYTES+0] = 1;
    buf[3].coeffs[OSKR_SYMBYTES+1] = 0;
  }
  else {
    buf[0].coeffs[OSKR_SYMBYTES+0] = 0;
    buf[0].coeffs[OSKR_SYMBYTES+1] = 0;
    buf[1].coeffs[OSKR_SYMBYTES+0] = 1;
    buf[1].coeffs[OSKR_SYMBYTES+1] = 0;
    buf[2].coeffs[OSKR_SYMBYTES+0] = 2;
    buf[2].coeffs[OSKR_SYMBYTES+1] = 0;
    buf[3].coeffs[OSKR_SYMBYTES+0] = 0;
    buf[3].coeffs[OSKR_SYMBYTES+1] = 1;
  }

  shake128x4_absorb_once(&state, buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, OSKR_SYMBYTES+2);
  shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, REJ_UNIFORM_AVX_NBLOCKS, &state);

  ctr0 = rej_uniform_avx(a[0].vec[0].coeffs, buf[0].coeffs);
  ctr1 = rej_uniform_avx(a[0].vec[1].coeffs, buf[1].coeffs);
  ctr2 = rej_uniform_avx(a[0].vec[2].coeffs, buf[2].coeffs);
  ctr3 = rej_uniform_avx(a[1].vec[0].coeffs, buf[3].coeffs);

  while(ctr0 < OSKR_N || ctr1 < OSKR_N || ctr2 < OSKR_N || ctr3 < OSKR_N) {
    shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 1, &state);

    ctr0 += rej_uniform(a[0].vec[0].coeffs + ctr0, OSKR_N - ctr0, buf[0].coeffs, SHAKE128_RATE);
    ctr1 += rej_uniform(a[0].vec[1].coeffs + ctr1, OSKR_N - ctr1, buf[1].coeffs, SHAKE128_RATE);
    ctr2 += rej_uniform(a[0].vec[2].coeffs + ctr2, OSKR_N - ctr2, buf[2].coeffs, SHAKE128_RATE);
    ctr3 += rej_uniform(a[1].vec[0].coeffs + ctr3, OSKR_N - ctr3, buf[3].coeffs, SHAKE128_RATE);
  }

  poly_nttunpack(&a[0].vec[0]);
  poly_nttunpack(&a[0].vec[1]);
  poly_nttunpack(&a[0].vec[2]);
  poly_nttunpack(&a[1].vec[0]);

  f = _mm256_loadu_si256((__m256i *)seed);
  _mm256_store_si256(buf[0].vec, f);
  _mm256_store_si256(buf[1].vec, f);
  _mm256_store_si256(buf[2].vec, f);
  _mm256_store_si256(buf[3].vec, f);

  if(transposed) {
    buf[0].coeffs[OSKR_SYMBYTES+0] = 1;
    buf[0].coeffs[OSKR_SYMBYTES+1] = 1;
    buf[1].coeffs[OSKR_SYMBYTES+0] = 1;
    buf[1].coeffs[OSKR_SYMBYTES+1] = 2;
    buf[2].coeffs[OSKR_SYMBYTES+0] = 2;
    buf[2].coeffs[OSKR_SYMBYTES+1] = 0;
    buf[3].coeffs[OSKR_SYMBYTES+0] = 2;
    buf[3].coeffs[OSKR_SYMBYTES+1] = 1;
  }
  else {
    buf[0].coeffs[OSKR_SYMBYTES+0] = 1;
    buf[0].coeffs[OSKR_SYMBYTES+1] = 1;
    buf[1].coeffs[OSKR_SYMBYTES+0] = 2;
    buf[1].coeffs[OSKR_SYMBYTES+1] = 1;
    buf[2].coeffs[OSKR_SYMBYTES+0] = 0;
    buf[2].coeffs[OSKR_SYMBYTES+1] = 2;
    buf[3].coeffs[OSKR_SYMBYTES+0] = 1;
    buf[3].coeffs[OSKR_SYMBYTES+1] = 2;
  }

  shake128x4_absorb_once(&state, buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, OSKR_SYMBYTES+2);
  shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, REJ_UNIFORM_AVX_NBLOCKS, &state);

  ctr0 = rej_uniform_avx(a[1].vec[1].coeffs, buf[0].coeffs);
  ctr1 = rej_uniform_avx(a[1].vec[2].coeffs, buf[1].coeffs);
  ctr2 = rej_uniform_avx(a[2].vec[0].coeffs, buf[2].coeffs);
  ctr3 = rej_uniform_avx(a[2].vec[1].coeffs, buf[3].coeffs);

  while(ctr0 < OSKR_N || ctr1 < OSKR_N || ctr2 < OSKR_N || ctr3 < OSKR_N) {
    shake128x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 1, &state);

    ctr0 += rej_uniform(a[1].vec[1].coeffs + ctr0, OSKR_N - ctr0, buf[0].coeffs, SHAKE128_RATE);
    ctr1 += rej_uniform(a[1].vec[2].coeffs + ctr1, OSKR_N - ctr1, buf[1].coeffs, SHAKE128_RATE);
    ctr2 += rej_uniform(a[2].vec[0].coeffs + ctr2, OSKR_N - ctr2, buf[2].coeffs, SHAKE128_RATE);
    ctr3 += rej_uniform(a[2].vec[1].coeffs + ctr3, OSKR_N - ctr3, buf[3].coeffs, SHAKE128_RATE);
  }

  poly_nttunpack(&a[1].vec[1]);
  poly_nttunpack(&a[1].vec[2]);
  poly_nttunpack(&a[2].vec[0]);
  poly_nttunpack(&a[2].vec[1]);

  f = _mm256_loadu_si256((__m256i *)seed);
  _mm256_store_si256(buf[0].vec, f);
  buf[0].coeffs[OSKR_SYMBYTES+0] = 2;
  buf[0].coeffs[OSKR_SYMBYTES+1] = 2;
  shake128_absorb_once(&state1x, buf[0].coeffs, OSKR_SYMBYTES+2);
  shake128_squeezeblocks(buf[0].coeffs, REJ_UNIFORM_AVX_NBLOCKS, &state1x);
  ctr0 = rej_uniform_avx(a[2].vec[2].coeffs, buf[0].coeffs);
  while(ctr0 < OSKR_N) {
    shake128_squeezeblocks(buf[0].coeffs, 1, &state1x);
    ctr0 += rej_uniform(a[2].vec[2].coeffs + ctr0, OSKR_N - ctr0, buf[0].coeffs, SHAKE128_RATE);
  }

  poly_nttunpack(&a[2].vec[2]);
}
#endif

/*************************************************
* Name:        indcpa_keypair
*
* Description: Generates public and private key for the CPA-secure
*              public-key encryption scheme underlying OSKR
*
* Arguments:   - uint8_t *pk: pointer to output public key
*                             (of length OSKR_INDCPA_PUBLICKEYBYTES bytes)
*              - uint8_t *sk: pointer to output private key
                              (of length OSKR_INDCPA_SECRETKEYBYTES bytes)
(pk,sk)
**************************************************/
void indcpa_keypair(uint8_t pk[OSKR_INDCPA_PUBLICKEYBYTES],
                    uint8_t sk[OSKR_INDCPA_SECRETKEYBYTES])
{
  unsigned int i;
  uint8_t buf[2*OSKR_SYMBYTES];
  const uint8_t *publicseed = buf;
  const uint8_t *noiseseed = buf + OSKR_SYMBYTES;
  polyvec a[OSKR_K], e, pkpv, skpv;

  randombytes(buf, OSKR_SYMBYTES);
  hash_g(buf, buf, OSKR_SYMBYTES);

  gen_a(a, publicseed);

#if OSKR_K == 2
  poly_getnoise_eta1_4x(skpv.vec+0, skpv.vec+1, e.vec+0, e.vec+1, noiseseed, 0, 1, 2, 3);
#elif OSKR_K == 3
  poly_getnoise_eta1_4x(skpv.vec+0, skpv.vec+1, skpv.vec+2, e.vec+0, noiseseed, 0, 1, 2, 3);
  poly_getnoise_eta1_4x(e.vec+1, e.vec+2, pkpv.vec+0, pkpv.vec+1, noiseseed, 4, 5, 6, 7);

#endif

  polyvec_ntt(&skpv);
  polyvec_reduce(&skpv);
  polyvec_ntt(&e);

  for(i=0;i<OSKR_K;i++) {
    polyvec_basemul_acc_montgomery(&pkpv.vec[i], &a[i], &skpv);
    poly_tomont(&pkpv.vec[i]);
  }

  polyvec_add(&pkpv, &pkpv, &e);
  polyvec_reduce(&pkpv);

  pack_sk(sk, &skpv);
  pack_pk(pk, &pkpv, publicseed);
}

/*************************************************
* Name:        indcpa_enc
*
* Description: Encryption function of the CPA-secure
*              public-key encryption scheme underlying OSKR.
*
* Arguments:   - uint8_t *c: pointer to output ciphertext
*                            (of length OSKR_INDCPA_BYTES bytes)
*              - const uint8_t *m: pointer to input message
*                                  (of length OSKR_INDCPA_MSGBYTES bytes)
*              - const uint8_t *pk: pointer to input public key
*                                   (of length OSKR_INDCPA_PUBLICKEYBYTES)
*              - const uint8_t *coins: pointer to input random coins used as seed
*                                      (of length OSKR_SYMBYTES) to deterministically
*                                      generate all randomness
**************************************************/
void indcpa_enc(uint8_t c[OSKR_INDCPA_BYTES],
                const uint8_t m[OSKR_INDCPA_MSGBYTES],
                const uint8_t pk[OSKR_INDCPA_PUBLICKEYBYTES],
                const uint8_t coins[OSKR_SYMBYTES])
{
  unsigned int i;
  uint8_t seed[OSKR_SYMBYTES];
  polyvec sp, pkpv, ep, at[OSKR_K], b;
  poly v, epp;

  unpack_pk(&pkpv, seed, pk);

  gen_at(at, seed);
  
#if (OSKR_K == 2 && OSKR_N == 256)
  poly_getnoise_eta1122_4x(sp.vec+0, sp.vec+1, ep.vec+0, ep.vec+1, coins, 0, 1, 2, 3);
  poly_getnoise_eta2(&epp, coins, 4);
#elif OSKR_K == 3
  poly_getnoise_eta1_4x(sp.vec+0, sp.vec+1, sp.vec+2, ep.vec+0, coins, 0, 1, 2 ,3);
  poly_getnoise_eta1_4x(ep.vec+1, ep.vec+2, &epp, b.vec+0, coins,  4, 5, 6, 7);
#elif (OSKR_K == 2 && OSKR_N == 512)
  poly_getnoise_eta1_4x(sp.vec+0, sp.vec+1, ep.vec+0, ep.vec+1, coins, 0, 1, 2, 3);
  poly_getnoise_eta2(&epp, coins, 4);
#endif
  polyvec_ntt(&sp);

  // matrix-vector multiplication
  for(i=0;i<OSKR_K;i++)
    polyvec_basemul_acc_montgomery(&b.vec[i], &at[i], &sp);
  polyvec_basemul_acc_montgomery(&v, &pkpv, &sp);

  polyvec_invntt_tomont(&b);
  poly_invntt_tomont(&v);

  polyvec_add(&b, &b, &ep);
  poly_add(&v, &v, &epp);
  poly_reduce(&v);
  polyvec_reduce(&b);
  polyvec_compress(c, &b);

  poly_con(c+OSKR_POLYVECCOMPRESSEDBYTES, m, &v);
}

/*************************************************
* Name:        indcpa_dec
*
* Description: Decryption function of the CPA-secure
*              public-key encryption scheme underlying OSKR.
*
* Arguments:   - uint8_t *m: pointer to output decrypted message
*                            (of length OSKR_INDCPA_MSGBYTES)
*              - const uint8_t *c: pointer to input ciphertext
*                                  (of length OSKR_INDCPA_BYTES)
*              - const uint8_t *sk: pointer to input secret key
*                                   (of length OSKR_INDCPA_SECRETKEYBYTES)
**************************************************/
void indcpa_dec(uint8_t m[OSKR_INDCPA_MSGBYTES],
                const uint8_t c[OSKR_INDCPA_BYTES],
                const uint8_t sk[OSKR_INDCPA_SECRETKEYBYTES])
{
  polyvec b, skpv;
  poly mp;

  polyvec_decompress(&b, c);
  unpack_sk(&skpv, sk);

  polyvec_ntt(&b);
  polyvec_basemul_acc_montgomery(&mp, &skpv, &b);
  poly_invntt_tomont(&mp);
  poly_reduce(&mp);

  poly_rec(m, c+OSKR_POLYVECCOMPRESSEDBYTES, &mp);
}
