#include <stdint.h>
#include <stddef.h>
#include "poly.h"
#include "ntt.h"
#include "polyvec.h"
#include "reduce.h"
#include "cbd.h"
#include "fips202.h"

/*************************************************
* Name:        poly_compress
*
* Description: Compression and subsequent serialization of a polynomial
*
* Arguments:   - unsigned char *r: pointer to output byte array
*              - poly *a:    pointer to input polynomial
**************************************************/
void poly_compress(unsigned char *r, poly *a)
{
  uint32_t t[8];
  unsigned int i,j,k=0;

  for(i=0;i<OKAI_N;i+=8)
  {
    for(j=0;j<8;j++)
      t[j] = (((freeze(a->coeffs[i+j]) << 4) + OKAI_Q/2)/OKAI_Q) & 15;

    r[k]   = t[0] | (t[1] << 4);
    r[k+1] = t[2] | (t[3] << 4);
    r[k+2] = t[4] | (t[5] << 4);
    r[k+3] = t[6] | (t[7] << 4);
    k += 4;
  }
}

/*************************************************
* Name:        poly_decompress
*
* Description: De-serialization and subsequent decompression of a polynomial;
*              approximate inverse of poly_compress
*
* Arguments:   - poly *r:                pointer to output polynomial
*              - const unsigned char *a: pointer to input byte array
**************************************************/
void poly_decompress(poly *r, const unsigned char *a)
{
  unsigned int i;
  for(i=0;i<OKAI_N;i+=8)
  {
    r->coeffs[i+0] = (((a[0] & 15) * OKAI_Q) + 8)>> 4;
    r->coeffs[i+1] = (((a[0] >> 4) * OKAI_Q)+ 8) >> 4;
    r->coeffs[i+2] = (((a[1] & 15) * OKAI_Q) + 8)>> 4;
    r->coeffs[i+3] = (((a[1] >> 4) * OKAI_Q) + 8)>> 4;
    r->coeffs[i+4] = (((a[2] & 15) * OKAI_Q) + 8)>> 4;
    r->coeffs[i+5] = (((a[2] >> 4) * OKAI_Q) + 8)>> 4;
    r->coeffs[i+6] = (((a[3] & 15) * OKAI_Q) + 8)>> 4;
    r->coeffs[i+7] = (((a[3] >> 4) * OKAI_Q) + 8)>> 4;
    a += 4;
  }
}

/*************************************************
* Name:        poly_tobytes
*
* Description: Serialization of a polynomial
*
* Arguments:   - unsigned char *r: pointer to output byte array
*              - const poly *a:    pointer to input polynomial
**************************************************/
void poly_tobytes(unsigned char *r, const poly *a)
{
  int16_t t[16];

  for(size_t i = 0; i < 16; i++) {
    for(size_t j = 0; j < 16; j++)
      t[j] = freeze(a->coeffs[16*j + i]);
    r[2*i +   0] = (t[0] >>  0);
    r[2*i +   1] = (t[0] >>  8) + (t[1] << 5);
    r[2*i +  32] = (t[1] >>  3);
    r[2*i +  33] = (t[1] >> 11) + (t[2]  << 2);
    r[2*i +  64] = (t[2] >>  6) + (t[3]  << 7);
    r[2*i +  65] = (t[3] >>  1);
    r[2*i +  96] = (t[3] >>  9) + (t[4]  << 4);
    r[2*i +  97] = (t[4] >>  4);
    r[2*i + 128] = (t[4] >> 12) + (t[5]  << 1);
    r[2*i + 129] = (t[5] >>  7) + (t[6]  << 6);
    r[2*i + 160] = (t[6] >>  2);
    r[2*i + 161] = (t[6] >> 10) + (t[7]  << 3);
    r[2*i + 192] = (t[7] >>  5);
    r[2*i + 193] = (t[8] >>  0);
    r[2*i + 224] = (t[8] >>  8) + (t[9] << 5);
    r[2*i + 225] = (t[9] >>  3);
    r[2*i + 256] = (t[9] >> 11) + (t[10]  << 2);
    r[2*i + 257] = (t[10] >>  6) + (t[11]  << 7);
    r[2*i + 288] = (t[11] >>  1);
    r[2*i + 289] = (t[11] >>  9) + (t[12]  << 4);
    r[2*i + 320] = (t[12] >>  4);
    r[2*i + 321] = (t[12] >> 12) + (t[13]  << 1);
    r[2*i + 352] = (t[13] >>  7) + (t[14]  << 6);
    r[2*i + 353] = (t[14] >>  2);
    r[2*i + 384] = (t[14] >> 10) + (t[15]  << 3);
    r[2*i + 385] = (t[15] >>  5);
  }
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
  unsigned char t[26];

  for(size_t i = 0; i < 16; i++) {
    for(size_t k = 0; k < 13; k++) {
      t[2*k]   = a[2*i + 32*k];
      t[2*k+1] = a[2*i + 32*k + 1];
    }
    r->coeffs[i +   0]  = t[0];
    r->coeffs[i +   0] += ((int16_t)t[1] & 0x1f) << 8;
    r->coeffs[i +  16]  = t[1] >> 5;
    r->coeffs[i +  16] += (int16_t)t[2] << 3;
    r->coeffs[i +  16] += ((int16_t)t[3] & 0x03) << 11;
    r->coeffs[i +  32]  = t[3] >> 2;
    r->coeffs[i +  32] += ((int16_t)t[4] & 0x7f) << 6;
    r->coeffs[i +  48]  = t[4] >> 7;
    r->coeffs[i +  48] += (int16_t)t[5] << 1;
    r->coeffs[i +  48] += ((int16_t)t[6] & 0x0f) <<  9;
    r->coeffs[i +  64]  = t[6] >> 4;
    r->coeffs[i +  64] += (int16_t)t[7] << 4;
    r->coeffs[i +  64] += ((int16_t)t[8] & 0x01) << 12;
    r->coeffs[i +  80]  = t[8] >> 1;
    r->coeffs[i +  80] += ((int16_t)t[9] & 0x3f) << 7;
    r->coeffs[i +  96]  = t[9] >> 6;
    r->coeffs[i +  96] += (int16_t)t[10] << 2;
    r->coeffs[i +  96] += ((int16_t)t[11] & 0x07) << 10;
    r->coeffs[i + 112]  = t[11] >> 3;
    r->coeffs[i + 112] += (int16_t)t[12] << 5;
    r->coeffs[i + 128]  = t[13];
    r->coeffs[i + 128] += ((int16_t)t[14] & 0x1f) << 8;
    r->coeffs[i + 144]  = t[14] >> 5;
    r->coeffs[i + 144] += (int16_t)t[15] << 3;
    r->coeffs[i + 144] += ((int16_t)t[16] & 0x03) << 11;
    r->coeffs[i + 160]  = t[16] >> 2;
    r->coeffs[i + 160] += ((int16_t)t[17] & 0x7f) << 6;
    r->coeffs[i + 176]  = t[17] >> 7;
    r->coeffs[i + 176] += (int16_t)t[18] << 1;
    r->coeffs[i + 176] += ((int16_t)t[19] & 0x0f) <<  9;
    r->coeffs[i + 192]  = t[19] >> 4;
    r->coeffs[i + 192] += (int16_t)t[20] << 4;
    r->coeffs[i + 192] += ((int16_t)t[21] & 0x01) << 12;
    r->coeffs[i + 208]  = t[21] >> 1;
    r->coeffs[i + 208] += ((int16_t)t[22] & 0x3f) << 7;
    r->coeffs[i + 224]  = t[22] >> 6;
    r->coeffs[i + 224] += (int16_t)t[23] << 2;
    r->coeffs[i + 224] += ((int16_t)t[24] & 0x07) << 10;
    r->coeffs[i + 240]  = t[24] >> 3;
    r->coeffs[i + 240] += (int16_t)t[25] << 5;
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
void poly_getnoise(poly *r,const unsigned char *seed, unsigned char nonce)
{
  unsigned char buf[OKAI_NETA*OKAI_N/4];
  unsigned char extseed[OKAI_SYMBYTES+1];
  unsigned int i;

  for(i=0;i<OKAI_SYMBYTES;i++)
    extseed[i] = seed[i];
  extseed[OKAI_SYMBYTES] = nonce;

  shake256(buf,OKAI_NETA*OKAI_N/4,extseed,OKAI_SYMBYTES+1);
  cbd4(r, buf);
}

void poly_getsecret(poly *r,const unsigned char *seed, unsigned char nonce)
{
  unsigned char buf[OKAI_NETA*OKAI_N/4];
  unsigned char extseed[OKAI_SYMBYTES+1];
  unsigned int i;

  for(i=0;i<OKAI_SYMBYTES;i++)
    extseed[i] = seed[i];
  extseed[OKAI_SYMBYTES] = nonce;

  shake256(buf,OKAI_NETA*OKAI_N/4,extseed,OKAI_SYMBYTES+1);
  cbd(r, buf);
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
void poly_invntt(poly *r)
{
  invntt(r->coeffs);
}

/*************************************************
* Name:        poly_basemul
*
* Description: Multiplication of two polynomials in NTT domain
*
* Arguments:   - poly *r:       pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
extern void basemul_asm(int16_t *, const int16_t *, const int16_t *, const int16_t *);
void poly_basemul(poly *r, const poly *a, const poly *b)
{
  basemul_asm(r->coeffs, a->coeffs, b->coeffs, zetas);
}

/*************************************************
* Name:        poly_reduce
*
* Description: Applies Barrett reduction to all coefficients of a polynomial
*              for details of the Barrett reduction see comments in reduce.c
*
* Arguments:   - poly *r:       pointer to input/output polynomial
**************************************************/
extern void asm_barrett_reduce(int16_t *r);
void poly_reduce(poly *r)
{
  asm_barrett_reduce(r->coeffs);
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
extern void pointwise_add(int16_t *, const int16_t *, const int16_t *);
void poly_add(poly *r, const poly *a, const poly *b)
{
  pointwise_add(r->coeffs,a->coeffs,b->coeffs);
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
extern void pointwise_sub(int16_t *, const int16_t *, const int16_t *);
void poly_sub(poly *r, const poly *a, const poly *b)
{
  pointwise_sub(r->coeffs,a->coeffs,b->coeffs);
}

/*************************************************
* Name:        poly_csubq
*
* Description: Applies conditional subtraction of q to each coefficient of a polynomial
*              for details of conditional subtraction of q see comments in reduce.c
*
* Arguments:   - poly *r:       pointer to input/output polynomial
**************************************************/
void poly_csubq(poly *r)
{
  int i;

  for(i=0;i<OKAI_N;i++)
    r->coeffs[i] = csubq(r->coeffs[i]);
}

/*************************************************
* Name:        poly_frommsg
*
* Description: Convert 32-byte message to polynomial
*
* Arguments:   - poly *r:                  pointer to output polynomial
*              - const unsigned char *msg: pointer to input message
**************************************************/
void poly_frommsg(poly *r, const unsigned char msg[OKAI_SYMBYTES])
{
  unsigned int i, j;
  int16_t mask;

  for(i=0;i<OKAI_SYMBYTES;i++)
  {
    for(j=0;j<8;j++)
    {
      mask = -((msg[i] >> j)&1);
      r->coeffs[8*i+j] = mask & ((OKAI_Q+1)/2);
    }
  }
}

/*************************************************
* Name:        poly_tomsg
*
* Description: Convert polynomial to 32-byte message
*
* Arguments:   - unsigned char *msg: pointer to output message
*              - poly *a:      pointer to input polynomial
**************************************************/
void poly_tomsg(unsigned char msg[OKAI_SYMBYTES], poly *a)
{
  uint16_t t;
  unsigned int i,j;

  for(i=0;i<OKAI_SYMBYTES;i++)
  {
    msg[i] = 0;
    for(j=0;j<8;j++)
    {
      t = (((freeze(a->coeffs[8*i+j]) << 1) + OKAI_Q/2)/OKAI_Q) & 1;
      msg[i] |= t << j;
      // msg[i] |= a->coeffs[8*i+j] << j;
    }
  }
}

void poly_con_4(uint8_t *c,
                       const uint8_t *m,
                       int16_t scoeffs[256]) {
    unsigned int i, j;
    uint8_t t[32];
    int16_t u, mask;
    uint32_t *pm = (uint32_t *)m;
    poly k;

    for(i=0;i<8;i++) {
        for(j=0;j<32;j++) {
            mask = -(int16_t) ((pm[i] >> j) & 1);// 0 or -1
            k.coeffs[8*j+i] = mask & ((OKAI_Q+1)/2);// 0 or Q/2 | round(q/m)·ki
            scoeffs[8*j+i] += k.coeffs[8*j+i];// σ2+Decompress(k,dm)
            // scoeffs[8*j+i] = barrett_reduce(scoeffs[8*j+i]);// σ2+Decompress(k,dm)
        }
    }
    for (i = 0; i < 8; i++) {
        for (j = 0; j < 32; j++) {
            // map to positive standard representatives
            u = scoeffs[8 * j + i];
            u += (u >> 15) & OKAI_Q;// map to [0,q-1]
            t[j] = ((((uint16_t) u << 4) + OKAI_Q / 2) / OKAI_Q) & 15;
        }
        for (j = 0; j < 32; j+=8) {
            c[4*i + 4*j + 0] = t[j + 0] | (t[j + 1] << 4);
            c[4*i + 4*j + 1] = t[j + 2] | (t[j + 3] << 4);
            c[4*i + 4*j + 2] = t[j + 4] | (t[j + 5] << 4);
            c[4*i + 4*j + 3] = t[j + 6] | (t[j + 7] << 4);
        }
    }
}

void poly_rec_4(uint8_t *m,
                       const uint8_t *c,
                       int16_t scoeffs[256]) {
    poly v;
    unsigned int i, j, k, x;
    int32_t t1;
    uint32_t *pm = (uint32_t *)m;
    unsigned char tmp[16];

    for(i = 0; i < 8; i++) {
        for(k = 0; k < 4; k++) {// recover v
            tmp[4*k]   = c[4*i + 32*k];
            tmp[4*k+1] = c[4*i + 32*k + 1];
            tmp[4*k+2] = c[4*i + 32*k + 2];
            tmp[4*k+3] = c[4*i + 32*k + 3];
        }
        for(x = 0; x < 16; x++) {
            v.coeffs[i + 16 * x]  = (tmp[x] & 15);
            v.coeffs[i + 16 * x + 8]  = (tmp[x] >> 4) & 15;
        }
    }
    for(i=0;i<8;i++) {
        pm[i] = 0;
        for(j=0;j<32;j++) {
            t1 = (v.coeffs[8*j+i] + 4) * OKAI_Q - (scoeffs[8*j+i] << 4);// qv - σ*g
            t1 = ((t1 * 8737) >> 26) + (t1 >> 31);
            t1 = (t1 >> 3) & 1;
            pm[i] |= t1 << j;
        }
    }
}
