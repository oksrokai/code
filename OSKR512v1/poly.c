#include "poly.h"

#include "cbd.h"
#include "ntt.h"
#include "params.h"
#include "symmetric.h"

#include <stdint.h>


/*************************************************
* Name:        poly_compress
*
* Description: Serialization of a polynomial and subsequent compression of a polynomial;
*
* Arguments:   - unsigned char *r: pointer to output byte array (of length OSKR_POLYCOMPRESSEDBYTES)
*              - const poly *a:    pointer to input polynomial to be serialized
*************************************************/
void poly_compress(unsigned char *r, poly *a)
{
  uint8_t t[8];
  int i,j,k=0;

#if (OSKR_POLYCOMPRESSEDBYTES == 128)
  for(i=0;i<OSKR_N;i+=8)
  {
    for(j=0;j<8;j++)
      t[j] = ((((uint32_t)a->coeffs[i+j] << 4) + OSKR_Q/2) / OSKR_Q) & 15;

    r[k]   = t[0] | (t[1] << 4);
    r[k+1] = t[2] | (t[3] << 4);
    r[k+2] = t[4] | (t[5] << 4);
    r[k+3] = t[6] | (t[7] << 4);
    k += 4;
  }
#elif (OSKR_POLYCOMPRESSEDBYTES == 160)
  for(i=0;i<OSKR_N;i+=8)
  {
    for(j=0;j<8;j++)
      t[j] = ((((uint32_t)a->coeffs[i+j] << 5) + OSKR_Q/2) / OSKR_Q) & 31;

    r[k]   =  t[0]       | (t[1] << 5);
    r[k+1] = (t[1] >> 3) | (t[2] << 2) | (t[3] << 7);
    r[k+2] = (t[3] >> 1) | (t[4] << 4);
    r[k+3] = (t[4] >> 4) | (t[5] << 1) | (t[6] << 6);
    r[k+4] = (t[6] >> 2) | (t[7] << 3);
    k += 5;
  }
#else
#error "OSKR_POLYCOMPRESSEDBYTES needs to be in {96, 128, 160}"
#endif
}

/*************************************************
* Name:        poly_decompress
*
* Description: De-serialization and subsequent decompression of a polynomial;
*              approximate inverse of poly_compress
*
* Arguments:   - poly *r:                pointer to output polynomial
*              - const unsigned char *a: pointer to input byte array (of length OSKR_POLYCOMPRESSEDBYTES bytes)
**************************************************/
void poly_decompress(poly *r, const unsigned char *a)
{
  int i;
#if (OSKR_POLYCOMPRESSEDBYTES == 128)
  for(i=0;i<OSKR_N;i+=8)
  {
    r->coeffs[i+0] = (((a[0] & 15) * OSKR_Q) + 8) >> 4;
    r->coeffs[i+1] = (((a[0] >> 4) * OSKR_Q) + 8) >> 4;
    r->coeffs[i+2] = (((a[1] & 15) * OSKR_Q) + 8) >> 4;
    r->coeffs[i+3] = (((a[1] >> 4) * OSKR_Q) + 8) >> 4;
    r->coeffs[i+4] = (((a[2] & 15) * OSKR_Q) + 8) >> 4;
    r->coeffs[i+5] = (((a[2] >> 4) * OSKR_Q) + 8) >> 4;
    r->coeffs[i+6] = (((a[3] & 15) * OSKR_Q) + 8) >> 4;
    r->coeffs[i+7] = (((a[3] >> 4) * OSKR_Q) + 8) >> 4;
    a += 4;
  }
#elif (OSKR_POLYCOMPRESSEDBYTES == 160)
  for(i=0;i<OSKR_N;i+=8)
  {
    r->coeffs[i+0] =  (((a[0] & 31) * OSKR_Q) + 16) >> 5;
    r->coeffs[i+1] = ((((a[0] >> 5) | ((a[1] & 3) << 3)) * OSKR_Q) + 16) >> 5;
    r->coeffs[i+2] = ((((a[1] >> 2) & 31) * OSKR_Q) + 16) >> 5;
    r->coeffs[i+3] = ((((a[1] >> 7) | ((a[2] & 15) << 1)) * OSKR_Q) + 16) >> 5;
    r->coeffs[i+4] = ((((a[2] >> 4) | ((a[3] &  1) << 4)) * OSKR_Q) + 16) >> 5;
    r->coeffs[i+5] = ((((a[3] >> 1) & 31) * OSKR_Q) + 16) >> 5;
    r->coeffs[i+6] = ((((a[3] >> 6) | ((a[4] &  7) << 2)) * OSKR_Q) + 16) >> 5;
    r->coeffs[i+7] =  (((a[4] >> 3) * OSKR_Q) + 16) >> 5;
    a += 5;
  }
#else
#error "OSKR_POLYCOMPRESSEDBYTES needs to be in {96, 128, 160}"
#endif
}

/*************************************************
* Name:        poly_packcompress
*
* Description: Serialization and subsequent compression of a polynomial of a polyvec,
*              writes to a byte string representation of the whole polyvec.
*              Used to compress a polyvec one poly at a time in a loop.
*
* Arguments:   - unsigned char *r:  pointer to output byte string representation of a polyvec (of length OSKR_POLYVECCOMPRESSEDBYTES)
*              - const poly *a:     pointer to input polynomial
*              - int i:             index of to be serialized polynomial in serialized polyec
**************************************************/
void poly_packcompress(unsigned char *r, poly *a, int i) {
    int j, k;

#if (OSKR_POLYVECCOMPRESSEDBYTES == (OSKR_K * 352))
  uint16_t t[8];

  for(j=0;j<OSKR_N/8;j++) {
    for(k=0;k<8;k++)
      t[k] = ((((uint32_t)a->coeffs[8*j+k] << 11) + OSKR_Q/2) / OSKR_Q) & 0x7ff;

    r[352*i+11*j+ 0] =  t[0] & 0xff;
    r[352*i+11*j+ 1] = (t[0] >>  8) | ((t[1] & 0x1f) << 3);
    r[352*i+11*j+ 2] = (t[1] >>  5) | ((t[2] & 0x03) << 6);
    r[352*i+11*j+ 3] = (t[2] >>  2) & 0xff;
    r[352*i+11*j+ 4] = (t[2] >> 10) | ((t[3] & 0x7f) << 1);
    r[352*i+11*j+ 5] = (t[3] >>  7) | ((t[4] & 0x0f) << 4);
    r[352*i+11*j+ 6] = (t[4] >>  4) | ((t[5] & 0x01) << 7);
    r[352*i+11*j+ 7] = (t[5] >>  1) & 0xff;
    r[352*i+11*j+ 8] = (t[5] >>  9) | ((t[6] & 0x3f) << 2);
    r[352*i+11*j+ 9] = (t[6] >>  6) | ((t[7] & 0x07) << 5);
    r[352*i+11*j+10] = (t[7] >>  3);
  }
#elif (OSKR_POLYVECCOMPRESSEDBYTES == (OSKR_K * 320))
    uint16_t t[4];

    for (j = 0; j < OSKR_N / 4; j++) {
        for (k = 0; k < 4; k++)
            t[k] = ((((uint32_t)a->coeffs[4 * j + k] << 10) + OSKR_Q / 2) / OSKR_Q) & 0x3ff;

        r[320*i+5*j+0] =   t[0] & 0xff;
        r[320*i+5*j+1] =  (t[0] >>  8) | ((t[1] & 0x3f) << 2);
        r[320*i+5*j+2] = ((t[1] >>  6) | ((t[2] & 0x0f) << 4)) & 0xff;
        r[320*i+5*j+3] = ((t[2] >>  4) | ((t[3] & 0x03) << 6)) & 0xff;
        r[320*i+5*j+4] =  (t[3] >>  2) & 0xff;
    }
#else
#error "OSKR_POLYVECCOMPRESSEDBYTES needs to in (OSKR_K * {352, 320})"
#endif
}

/*************************************************
* Name:        poly_unpackdecompress
*
* Description: Deserialization and subsequent compression of a polynomial of a polyvec,
*              Used to uncompress a polyvec one poly at a time in a loop.
*
* Arguments:   - const poly *r:     pointer to output polynomial
*              - unsigned char *a:  pointer to input byte string representation of a polyvec (of length OSKR_POLYVECCOMPRESSEDBYTES)
*              - int i:             index of poly in polyvec to decompress
**************************************************/
void poly_unpackdecompress(poly *r, const unsigned char *a, int i) {
  int j;
#if (OSKR_POLYVECCOMPRESSEDBYTES == (OSKR_K * 352))
    for(j=0;j<OSKR_N/8;j++)
    {
      r->coeffs[8*j+0] =  (((a[352*i+11*j+ 0]       | (((uint32_t)a[352*i+11*j+ 1] & 0x07) << 8)) * OSKR_Q) + 1024) >> 11;
      r->coeffs[8*j+1] = ((((a[352*i+11*j+ 1] >> 3) | (((uint32_t)a[352*i+11*j+ 2] & 0x3f) << 5)) * OSKR_Q) + 1024) >> 11;
      r->coeffs[8*j+2] = ((((a[352*i+11*j+ 2] >> 6) | (((uint32_t)a[352*i+11*j+ 3] & 0xff) << 2) | (((uint32_t)a[352*i+11*j+4] & 0x01) << 10)) * OSKR_Q) + 1024) >> 11;
      r->coeffs[8*j+3] = ((((a[352*i+11*j+ 4] >> 1) | (((uint32_t)a[352*i+11*j+ 5] & 0x0f) << 7)) * OSKR_Q) + 1024) >> 11;
      r->coeffs[8*j+4] = ((((a[352*i+11*j+ 5] >> 4) | (((uint32_t)a[352*i+11*j+ 6] & 0x7f) << 4)) * OSKR_Q) + 1024) >> 11;
      r->coeffs[8*j+5] = ((((a[352*i+11*j+ 6] >> 7) | (((uint32_t)a[352*i+11*j+ 7] & 0xff) << 1) | (((uint32_t)a[352*i+11*j+8] & 0x03) <<  9)) * OSKR_Q) + 1024) >> 11;
      r->coeffs[8*j+6] = ((((a[352*i+11*j+ 8] >> 2) | (((uint32_t)a[352*i+11*j+ 9] & 0x1f) << 6)) * OSKR_Q) + 1024) >> 11;
      r->coeffs[8*j+7] = ((((a[352*i+11*j+ 9] >> 5) | (((uint32_t)a[352*i+11*j+10] & 0xff) << 3)) * OSKR_Q) + 1024) >> 11;
    }
#elif (OSKR_POLYVECCOMPRESSEDBYTES == (OSKR_K * 320))
    for(j=0;j<OSKR_N/4;j++)
    {
      r->coeffs[4*j+0] =  (((a[320*i+5*j+ 0]       | (((uint32_t)a[320*i+5*j+ 1] & 0x03) << 8)) * OSKR_Q) + 512) >> 10;
      r->coeffs[4*j+1] = ((((a[320*i+5*j+ 1] >> 2) | (((uint32_t)a[320*i+5*j+ 2] & 0x0f) << 6)) * OSKR_Q) + 512) >> 10;
      r->coeffs[4*j+2] = ((((a[320*i+5*j+ 2] >> 4) | (((uint32_t)a[320*i+5*j+ 3] & 0x3f) << 4)) * OSKR_Q) + 512) >> 10;
      r->coeffs[4*j+3] = ((((a[320*i+5*j+ 3] >> 6) | (((uint32_t)a[320*i+5*j+ 4] & 0xff) << 2)) * OSKR_Q) + 512) >> 10;
    }
#else
#error "OSKR_POLYVECCOMPRESSEDBYTES needs to be in {320*OSKR_K, 352*OSKR_K}"
#endif
}


/*************************************************
* Name:        cmp_poly_compress
*
* Description: Serializes and consequently compares polynomial to a serialized polynomial
*
* Arguments:   - const unsigned char *r:    pointer to serialized polynomial to compare with
*              - poly *a:                   pointer to input polynomial to serialize and compare
* Returns:                                  boolean indicating whether the polynomials are equal
**************************************************/
int cmp_poly_compress(const unsigned char *r, poly *a) {
    unsigned char rc = 0;
    uint8_t t[8];
    int i, j, k = 0;

#if (OSKR_POLYCOMPRESSEDBYTES == 128)
    for (i = 0; i < OSKR_N; i += 8) {
        for (j = 0; j < 8; j++)
            t[j] = ((((uint32_t)a->coeffs[i + j] << 4) + OSKR_Q / 2) / OSKR_Q) & 15;

        rc |= r[k]      ^ (t[0] | (t[1] << 4));
        rc |= r[k + 1]  ^ (t[2] | (t[3] << 4));
        rc |= r[k + 2]  ^ (t[4] | (t[5] << 4));
        rc |= r[k + 3]  ^ (t[6] | (t[7] << 4));
        k += 4;
    }
#elif (OSKR_POLYCOMPRESSEDBYTES == 160)
    for(i=0;i<OSKR_N;i+=8)
    {
      for(j=0;j<8;j++)
        t[j] = ((((uint32_t)a->coeffs[i+j] << 5) + OSKR_Q/2) / OSKR_Q) & 31;

      rc |= r[k]   ^ (t[0]       | (t[1] << 5));
      rc |= r[k+1] ^ ((t[1] >> 3) | (t[2] << 2) | (t[3] << 7));
      rc |= r[k+2] ^ ((t[3] >> 1) | (t[4] << 4));
      rc |= r[k+3] ^ ((t[4] >> 4) | (t[5] << 1) | (t[6] << 6));
      rc |= r[k+4] ^ ((t[6] >> 2) | (t[7] << 3));
      k += 5;
    }
#else
#error "OSKR_POLYCOMPRESSEDBYTES needs to be in {96, 128, 160}"
#endif
    return rc;
}

/*************************************************
* Name:        cmp_poly_packcompress
*
* Description: Serializes and consequently compares poly of polyvec to a serialized polyvec
*              Should be called in a loop over all poly's of a polyvec.
*
* Arguments:   - const unsigned char *r:    pointer to serialized polyvec to compare with
*              - poly *a:                   pointer to input polynomial of polyvec to serialize and compare
*              - int i:                     index of poly in polyvec to compare with
* Returns:                                  boolean indicating whether the polyvecs are equal
**************************************************/
int cmp_poly_packcompress(const unsigned char *r, poly *a, int i) {
    unsigned char rc = 0;
    int j, k;

#if (OSKR_POLYVECCOMPRESSEDBYTES == (OSKR_K * 352))
  uint16_t t[8];
    for(j=0;j<OSKR_N/8;j++)
    {
      for(k=0;k<8;k++)
        t[k] = ((((uint32_t)a->coeffs[8*j+k] << 11) + OSKR_Q/2) / OSKR_Q) & 0x7ff;

      rc |= r[352*i+11*j+ 0] ^ (t[0] & 0xff);
      rc |= r[352*i+11*j+ 1] ^ ((t[0] >>  8) | ((t[1] & 0x1f) << 3));
      rc |= r[352*i+11*j+ 2] ^ ((t[1] >>  5) | ((t[2] & 0x03) << 6));
      rc |= r[352*i+11*j+ 3] ^ ((t[2] >>  2) & 0xff);
      rc |= r[352*i+11*j+ 4] ^ ((t[2] >> 10) | ((t[3] & 0x7f) << 1));
      rc |= r[352*i+11*j+ 5] ^ ((t[3] >>  7) | ((t[4] & 0x0f) << 4));
      rc |= r[352*i+11*j+ 6] ^ ((t[4] >>  4) | ((t[5] & 0x01) << 7));
      rc |= r[352*i+11*j+ 7] ^ ((t[5] >>  1) & 0xff);
      rc |= r[352*i+11*j+ 8] ^ ((t[5] >>  9) | ((t[6] & 0x3f) << 2));
      rc |= r[352*i+11*j+ 9] ^ ((t[6] >>  6) | ((t[7] & 0x07) << 5));
      rc |= r[352*i+11*j+10] ^ ((t[7] >>  3));
    }
#elif (OSKR_POLYVECCOMPRESSEDBYTES == (OSKR_K * 320))
    uint16_t t[4];
        for (j = 0; j < OSKR_N / 4; j++) {
            for (k = 0; k < 4; k++)
                t[k] = ((((uint32_t)a->coeffs[4 * j + k] << 10) + OSKR_Q / 2) / OSKR_Q) & 0x3ff;

            rc |= r[320*i+5*j+0] ^ (t[0] & 0xff);
            rc |= r[320*i+5*j+1] ^ ((t[0] >>  8) | ((t[1] & 0x3f) << 2));
            rc |= r[320*i+5*j+2] ^ (((t[1] >>  6) | ((t[2] & 0x0f) << 4)) & 0xff);
            rc |= r[320*i+5*j+3] ^ (((t[2] >>  4) | ((t[3] & 0x03) << 6)) & 0xff);
            rc |= r[320*i+5*j+4] ^ ((t[3] >>  2) & 0xff);
        }
#else
#error "OSKR_POLYVECCOMPRESSEDBYTES needs to be in {320*OSKR_K, 352*OSKR_K}"
#endif
    return rc;
}

/*************************************************
* Name:        poly_tobytes
*
* Description: Serialization of a polynomial
*
* Arguments:   - unsigned char *r: pointer to output byte array (needs space for OSKR_POLYBYTES bytes)
*              - const poly *a:    pointer to input polynomial
**************************************************/
void poly_tobytes(unsigned char *r, poly *a) {
    int i;
    uint16_t t0, t1;

    poly_reduce(a);

    for (i = 0; i < OSKR_N / 2; i++) {
        t0 = a->coeffs[2 * i];
        t1 = a->coeffs[2 * i + 1];
        r[3 * i] = t0 & 0xff;
        r[3 * i + 1] = (t0 >> 8) | ((t1 & 0xf) << 4);
        r[3 * i + 2] = (t1 >> 4) & 0xff;
    }
}

/*************************************************
* Name:        poly_frombytes
*
* Description: De-serialization of a polynomial;
*              inverse of poly_tobytes
*
* Arguments:   - poly *r:                pointer to output polynomial
*              - const unsigned char *a: pointer to input byte array (of OSKR_POLYBYTES bytes)
**************************************************/
void poly_frombytes(poly *r, const unsigned char *a) {
    int i;

    for (i = 0; i < OSKR_N / 2; i++) {
        r->coeffs[2 * i]     = a[3 * i]          | ((uint16_t)a[3 * i + 1] & 0x0f) << 8;
        r->coeffs[2 * i + 1] = a[3 * i + 1] >> 4 | ((uint16_t)a[3 * i + 2] & 0xff) << 4;
    }
}

/*************************************************
* Name:        poly_frombytes_mul_16_32
*
* Description: Multiplication of a polynomial with a de-serialization of another polynomial
*              Using strategy of better accumulation.
* Arguments:   - const poly *b:          pointer to input polynomial
*              - int32_t *r_tmp:         array for accumulating unreduced results
*              - const unsigned char *a: pointer to input byte array (of OSKR_POLYBYTES bytes)
**************************************************/
extern void frombytes_mul_asm_16_32(int32_t *r_tmp, const int16_t *b, const unsigned char *c, const int32_t zetas[64]);
void poly_frombytes_mul_16_32(int32_t *r_tmp, const poly *b, const unsigned char *a) {
    frombytes_mul_asm_16_32(r_tmp, b->coeffs, a, zetas);
}

/*************************************************
* Name:        poly_frombytes_mul_32_32
*
* Description: Multiplication of a polynomial with a de-serialization of another polynomial
*              Using strategy of better accumulation.
* Arguments:   - const poly *b:          pointer to input polynomial
*              - int32_t *r_tmp:         array for accumulating unreduced results
*              - const unsigned char *a: pointer to input byte array (of OSKR_POLYBYTES bytes)
**************************************************/
extern void frombytes_mul_asm_acc_32_32(int32_t *r_tmp, const int16_t *b, const unsigned char *c, const int32_t zetas[64]);
void poly_frombytes_mul_32_32(int32_t *r_tmp, const poly *b, const unsigned char *a) {
    frombytes_mul_asm_acc_32_32(r_tmp, b->coeffs, a, zetas);
}

/*************************************************
* Name:        poly_frombytes_mul_32_16
*
* Description: Multiplication of a polynomial with a de-serialization of another polynomial
*              Using strategy of better accumulation.
* Arguments:   - poly *r:                pointer to output polynomial
*              - const poly *b:          pointer to input polynomial
*              - const int32_t *r_tmp:   array containing unreduced results
*              - const unsigned char *a: pointer to input byte array (of OSKR_POLYBYTES bytes)
**************************************************/
extern void frombytes_mul_asm_acc_32_16(int16_t *r, const int16_t *b, const unsigned char *c, const int32_t zetas[64], const int32_t *r_tmp);
void poly_frombytes_mul_32_16(poly *r, const poly* b, const unsigned char *a, const int32_t *r_tmp) {
    frombytes_mul_asm_acc_32_16(r->coeffs, b->coeffs, a, zetas, r_tmp);
}

/*************************************************
* Name:        poly_getnoise_eta1
*
* Description: Sample a polynomial deterministically from a seed and a nonce,
*              with output polynomial close to centered binomial distribution
*              with parameter OSKR_ETA1
*
* Arguments:   - poly *r:                   pointer to output polynomial
*              - const unsigned char *seed: pointer to input seed (pointing to array of length OSKR_SYMBYTES bytes)
*              - unsigned char nonce:       one-byte input nonce
*              - int add:                   boolean to indicate to accumulate into r
**************************************************/
void poly_noise_eta1(poly *r, const unsigned char *seed, unsigned char nonce, int add) {
    unsigned char buf[OSKR_ETA1 * OSKR_N / 4];

    prf(buf, OSKR_ETA1 * OSKR_N / 4, seed, nonce);
    cbd_eta1(r, buf, add);
}

/*************************************************
* Name:        poly_getnoise_eta2
*
* Description: Sample a polynomial deterministically from a seed and a nonce,
*              with output polynomial close to centered binomial distribution
*              with parameter OSKR_ETA2
*
* Arguments:   - poly *r:                   pointer to output polynomial
*              - const unsigned char *seed: pointer to input seed (pointing to array of length OSKR_SYMBYTES bytes)
*              - unsigned char nonce:       one-byte input nonce
*              - int add:                   boolean to indicate to accumulate into r
**************************************************/
void poly_noise_eta2(poly *r, const unsigned char *seed, unsigned char nonce, int add) {
    unsigned char buf[OSKR_ETA2 * OSKR_N / 4];

    prf(buf, OSKR_ETA2 * OSKR_N / 4, seed, nonce);
    cbd_eta2(r, buf, add);
}

/*************************************************
* Name:        poly_basemul_opt_16_32
*
* Description: Multiplication of two polynomials using asymmetric multiplication.
*              Cached values are generated during matrix-vector product.
*              Using strategy of better accumulation (initial step).
* Arguments:   - const poly *a:       pointer to input polynomial
*              - const poly *b:       pointer to input polynomial
*              - const poly *a_prime: pointer to a pre-multiplied by zetas 
*              - int32_t *r_tmp:      array for accumulating unreduced results
**************************************************/
extern void basemul_asm_opt_16_32(int32_t *, const int16_t *, const int16_t *, const int16_t *);
void poly_basemul_opt_16_32(int32_t *r_tmp, const poly *a, const poly *b, const poly *a_prime) {
    basemul_asm_opt_16_32(r_tmp, a->coeffs, b->coeffs, a_prime->coeffs);
}

/*************************************************
* Name:        poly_basemul_acc_opt_32_32
*
* Description: Multiplication of two polynomials using asymmetric multiplication.
*              Cached values are generated during matrix-vector product.
*              Using strategy of better accumulation.
* Arguments:   - const poly *a:       pointer to input polynomial
*              - const poly *b:       pointer to input polynomial
*              - const poly *a_prime: pointer to a pre-multiplied by zetas 
*              - int32_t *r_tmp:      array for accumulating unreduced results
**************************************************/
extern void basemul_asm_acc_opt_32_32(int32_t *, const int16_t *, const int16_t *, const int16_t *);
void poly_basemul_acc_opt_32_32(int32_t *r, const poly *a, const poly *b, const poly *a_prime) {
    basemul_asm_acc_opt_32_32(r, a->coeffs, b->coeffs, a_prime->coeffs);
}

/*************************************************
* Name:        poly_basemul_acc_opt_32_16
*
* Description: Multiplication of two polynomials using asymmetric multiplication.
*              Cached values are generated during matrix-vector product.
*              Using strategy of better accumulation (final step).
* Arguments:   - const poly *a:        pointer to input polynomial
*              - const poly *b:        pointer to input polynomial
*              - const poly *a_prime:  pointer to a pre-multiplied by zetas 
*              - poly *r:              pointer to output polynomial
*              - const int32_t *r_tmp: array containing unreduced results
**************************************************/
extern void basemul_asm_acc_opt_32_16(int16_t *, const int16_t *, const int16_t *, const int16_t *, const int32_t *);
void poly_basemul_acc_opt_32_16(poly *r, const poly *a, const poly *b, const poly *a_prime, const int32_t * r_tmp) {
    basemul_asm_acc_opt_32_16(r->coeffs, a->coeffs, b->coeffs, a_prime->coeffs, r_tmp);
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
void poly_ntt(poly *r) {
    ntt(r->coeffs);
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
void poly_invntt(poly *r) {
    invntt(r->coeffs);
}

extern void asm_fromplant(int16_t *r);
/*************************************************
* Name:        poly_fromplantt
*
* Description: Inplace conversion of all coefficients of a polynomial
*              from Montgomery domain to normal domain
*
* Arguments:   - poly *r:       pointer to input/output polynomial
**************************************************/
void poly_fromplant(poly *r) {
  asm_fromplant(r->coeffs);
}

extern void asm_barrett_reduce(int16_t *r);
/*************************************************
* Name:        poly_reduce
*
* Description: Applies Barrett reduction to all coefficients of a polynomial
*              for details of the Barrett reduction see comments in reduce.c
*
* Arguments:   - poly *r:       pointer to input/output polynomial
**************************************************/
void poly_reduce(poly *r) {
  asm_barrett_reduce(r->coeffs);
}

extern void pointwise_add(int16_t *, const int16_t *, const int16_t *);
/*************************************************
* Name:        poly_add
*
* Description: Add two polynomials
*
* Arguments: - poly *r:       pointer to output polynomial
*            - const poly *a: pointer to first input polynomial
*            - const poly *b: pointer to second input polynomial
**************************************************/
void poly_add(poly *r, const poly *a, const poly *b) {
    pointwise_add(r->coeffs,a->coeffs,b->coeffs);
}


extern void pointwise_sub(int16_t *, const int16_t *, const int16_t *);
/*************************************************
* Name:        poly_sub
*
* Description: Subtract two polynomials
*
* Arguments: - poly *r:       pointer to output polynomial
*            - const poly *a: pointer to first input polynomial
*            - const poly *b: pointer to second input polynomial
**************************************************/
void poly_sub(poly *r, const poly *a, const poly *b) {
    pointwise_sub(r->coeffs,a->coeffs,b->coeffs);
}

/*************************************************
* Name:        poly_frommsg
*
* Description: Convert 32-byte message to polynomial
*
* Arguments:   - poly *r:                  pointer to output polynomial
*              - const unsigned char *msg: pointer to input message
**************************************************/
void poly_frommsg(poly *r, const unsigned char msg[OSKR_SYMBYTES]) {
    int i, j;
    uint16_t mask;

    for (i = 0; i < OSKR_SYMBYTES; i++) {
        for (j = 0; j < 8; j++) {
            mask = -((msg[i] >> j) & 1);
            r->coeffs[8 * i + j] = mask & ((OSKR_Q + 1) / 2);
        }
    }
}

/*************************************************
* Name:        poly_tomsg
*
* Description: Convert polynomial to 32-byte message
*
* Arguments:   - unsigned char *msg: pointer to output message
*              - const poly *a:      pointer to input polynomial
**************************************************/
void poly_tomsg(unsigned char msg[OSKR_SYMBYTES], poly *a) {
    uint16_t t;
    int i, j;

    for (i = 0; i < OSKR_SYMBYTES; i++) {
        msg[i] = 0;
        for (j = 0; j < 8; j++) {
            t = (((a->coeffs[8 * i + j] << 1) + OSKR_Q / 2) / OSKR_Q) & 1;
            msg[i] |= t << j;
        }
    }
}

/*************************************************
* Name:        poly_zeroize
*
* Description: Zeros a polynomial
*
* Arguments:   - poly *p: pointer to polynomial
**************************************************/
void poly_zeroize(poly *p) {
  int i;
  for(i = 0; i < OSKR_N; i++)
   p->coeffs[i] = 0;
}

extern void con(uint8_t *, uint32_t *, int16_t *);
void poly_con(uint8_t c[OSKR_POLYCOMPRESSEDBYTES],
              const uint8_t m[OSKR_INDCPA_MSGBYTES],
              poly *s)
{
  size_t i,j;
  int16_t sp[32];
  uint32_t *pm = (uint32_t *)m;

  for(i=0;i<OSKR_N/32;i++) {
    for(j=0;j<32;j++) {
      sp[j] = s->coeffs[8*j+i];
    }
    con(c, &pm[i], sp);
    c += 4;
  }
}

int cmp_poly_con(const uint8_t c[OSKR_POLYCOMPRESSEDBYTES],
              const uint8_t m[OSKR_INDCPA_MSGBYTES],
              poly *s)
{
  size_t i,j;
  int k=0;
  uint8_t ctmp[128];
  int16_t sp[32];
  uint32_t *pm = (uint32_t *)m;
  unsigned char rc=0;

  for(i=0;i<OSKR_N/32;i++) {
    for(j=0;j<32;j++) {
      sp[j] = s->coeffs[8*j+i];
    }
    con(ctmp+4*i, &pm[i], sp);
    rc |= c[k] ^ ctmp[k];
    rc |= c[k+1] ^ ctmp[k+1];
    rc |= c[k+2] ^ ctmp[k+2];
    rc |= c[k+3] ^ ctmp[k+3];

    rc |= c[k+32] ^ ctmp[k+32];
    rc |= c[k+33] ^ ctmp[k+33];
    rc |= c[k+34] ^ ctmp[k+34];
    rc |= c[k+35] ^ ctmp[k+35];

    rc |= c[k+64] ^ ctmp[k+64];
    rc |= c[k+65] ^ ctmp[k+65];
    rc |= c[k+66] ^ ctmp[k+66];
    rc |= c[k+67] ^ ctmp[k+67];

    rc |= c[k+96] ^ ctmp[k+96];
    rc |= c[k+97] ^ ctmp[k+97];
    rc |= c[k+98] ^ ctmp[k+98];
    rc |= c[k+99] ^ ctmp[k+99];
    k += 4;
    //ctmp+=4;
  }
  return rc;
}

extern void rec(uint32_t *, const uint8_t *, int16_t *);
void poly_rec(uint8_t m[OSKR_INDCPA_MSGBYTES],
              const uint8_t c[OSKR_POLYCOMPRESSEDBYTES],
              poly *s)
{
  
  int16_t sp[32];
  unsigned int i,j;
  uint32_t *pm = (uint32_t *)m;

  for(i=0;i<OSKR_N/32;i++) {
    for(j=0;j<32;j++)
    {
      sp[j] = s->coeffs[8*j+i];
    }
    rec(&pm[i], c, sp);
    c += 4;
  }
}
