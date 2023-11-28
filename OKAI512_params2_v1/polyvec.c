#include <stdio.h>
#include "polyvec.h"
#include "fips202.h"
#include "cbd.h"
#include "reduce.h"
#include "poly.h"

/*************************************************
* Name:        polyvec_compress
*
* Description: Compress and serialize vector of polynomials
*
* Arguments:   - unsigned char *r: pointer to output byte array
*              - const polyvec *a: pointer to input vector of polynomials
**************************************************/
static void poly_compress2(uint8_t *r, const poly *a)
{
  unsigned int i,j;
  uint16_t t[16];
  for(i = 0; i < 16; i++) {
    for(j = 0; j < 16; j++)
      t[j] = ((((uint32_t)(a->coeffs[16*j + i]) << 9) + OKAI_Q/2)/ OKAI_Q) & 0x1ff;
      
    r[2*i +   0] =  t[0] & 0xff;
    r[2*i +   1] = (t[0] >> 8) | ((t[1] & 0x7f) << 1);
    r[2*i +  32] = (t[1] >> 7) | ((t[2] & 0x3f) << 2);
    r[2*i +  33] = (t[2] >> 6) | ((t[3] & 0x1f) << 3);
    r[2*i +  64] = (t[3] >> 5) | ((t[4] & 0x0f) << 4);
    r[2*i +  65] = (t[4] >> 4) | ((t[5] & 0x07) << 5);
    r[2*i +  96] = (t[5] >> 3) | ((t[6] & 0x03) << 6);
    r[2*i +  97] = (t[6] >> 2) | ((t[7] & 0x01) << 7);
    r[2*i + 128] =  t[7] >> 1;

    r[2*i + 129] =  t[8] & 0xff;
    r[2*i + 160] = (t[8] >> 8) | ((t[9] & 0x7f) << 1);
    r[2*i + 161] = (t[9] >> 7) | ((t[10] & 0x3f) << 2);
    r[2*i + 192] = (t[10] >> 6) | ((t[11] & 0x1f) << 3);
    r[2*i + 193] = (t[11] >> 5) | ((t[12] & 0x0f) << 4);
    r[2*i + 224] = (t[12] >> 4) | ((t[13] & 0x07) << 5);
    r[2*i + 225] = (t[13] >> 3) | ((t[14] & 0x03) << 6);
    r[2*i + 256] = (t[14] >> 2) | ((t[15] & 0x01) << 7);
    r[2*i + 257] =  t[15] >> 1;
  }
}

static void poly_decompress2(poly *r, const uint8_t *a)
{
  unsigned char t[18];
  unsigned int i, k;

  for(i = 0; i < 16; i++) {
    for(k = 0; k < 9; k++) {
      t[2*k]   = a[2*i + 32*k];
      t[2*k+1] = a[2*i + 32*k + 1];
    }

    r->coeffs[i +   0]  = t[0];
    r->coeffs[i +   0] += ((int16_t)t[1] & 0x01) << 8;
    r->coeffs[i +  16]  = t[1] >> 1;
    r->coeffs[i +  16] += ((int16_t)t[2] & 0x03) << 7;
    r->coeffs[i +  32]  = t[2] >> 2;
    r->coeffs[i +  32] += ((int16_t)t[3] & 0x07) << 6;
    r->coeffs[i +  48]  = t[3] >> 3;
    r->coeffs[i +  48] += ((int16_t)t[4] & 0x0f) << 5;
    r->coeffs[i +  64]  = t[4] >> 4;
    r->coeffs[i +  64] += ((int16_t)t[5] & 0x1f) << 4;
    r->coeffs[i +  80]  = t[5] >> 5;
    r->coeffs[i +  80] += ((int16_t)t[6] & 0x3f) << 3;
    r->coeffs[i +  96]  = t[6] >> 6;
    r->coeffs[i +  96] += ((int16_t)t[7] & 0x7f) << 2;
    r->coeffs[i + 112]  = t[7] >> 7;
    r->coeffs[i + 112] += ((int16_t)t[8] & 0xff) << 1;

    r->coeffs[i + 128]  = t[9];
    r->coeffs[i + 128] += ((int16_t)t[10] & 0x01) << 8;
    r->coeffs[i + 144]  = t[10] >> 1;
    r->coeffs[i + 144] += ((int16_t)t[11] & 0x03) << 7;
    r->coeffs[i + 160]  = t[11] >> 2;
    r->coeffs[i + 160] += ((int16_t)t[12] & 0x07) << 6;
    r->coeffs[i + 176]  = t[12] >> 3;
    r->coeffs[i + 176] += ((int16_t)t[13] & 0x0f) << 5;
    r->coeffs[i + 192]  = t[13] >> 4;
    r->coeffs[i + 192] += ((int16_t)t[14] & 0x1f) << 4;
    r->coeffs[i + 208]  = t[14] >> 5;
    r->coeffs[i + 208] += ((int16_t)t[15] & 0x3f) << 3;
    r->coeffs[i + 224]  = t[15] >> 6;
    r->coeffs[i + 224] += ((int16_t)t[16] & 0x7f) << 2;
    r->coeffs[i + 240]  = t[16] >> 7;
    r->coeffs[i + 240] += ((int16_t)t[17] & 0xff) << 1;
  }

  for(i = 0; i < OKAI_N; i++) {
    r->coeffs[i] = ((r->coeffs[i] * OKAI_Q) + 256) >> 9;
  }
}

void polyvec_compresspk(unsigned char *r, const polyvec *a)
{
  int i;
  for(i=0;i<OKAI_K;i++)
    poly_compress2(r+288*i, &a->vec[i]);
}

void polyvec_compressc(unsigned char *r, const polyvec *a)
{
  int i, j, k;
  uint8_t t[16];
  for(i = 0; i < OKAI_K; i++)
  {
    for(j = 0; j < 16; j++) {
      for(k = 0; k < 16; k++)
        t[k] = ((((uint32_t)(a->vec[i].coeffs[16*k + j]) << 8) + OKAI_Q/2)/ OKAI_Q) & 0xff;
      r[2*j +   0] = t[0];
      r[2*j +   1] = t[1];
      r[2*j +  32] = t[2];
      r[2*j +  33] = t[3];
      r[2*j +  64] = t[4];
      r[2*j +  65] = t[5];
      r[2*j +  96] = t[6];
      r[2*j +  97] = t[7];
      r[2*j + 128] = t[8];
      r[2*j + 129] = t[9];
      r[2*j + 160] = t[10];
      r[2*j + 161] = t[11];
      r[2*j + 192] = t[12];
      r[2*j + 193] = t[13];
      r[2*j + 224] = t[14];
      r[2*j + 225] = t[15];
    }
    r += 256;
  }
}

/*************************************************
* Name:        polyvec_decompress
*
* Description: De-serialize and decompress vector of polynomials;
*              approximate inverse of polyvec_compress
*
* Arguments:   - polyvec *r:       pointer to output vector of polynomials
*              - unsigned char *a: pointer to input byte array
**************************************************/
void polyvec_decompressc(polyvec *r, const unsigned char *a)
{
  int i, j, k;
  unsigned char t[16];

  for(i = 0; i < OKAI_K; i++)
  {
    for(j = 0; j < 16; j++) {
      for(k = 0; k < 8; k++) {
        t[2*k]   = a[2*j + 32*k];
        t[2*k+1] = a[2*j + 32*k + 1];
      }
      r->vec[i].coeffs[j +   0]  = t[0];
      r->vec[i].coeffs[j +  16]  = t[1];
      r->vec[i].coeffs[j +  32]  = t[2];
      r->vec[i].coeffs[j +  48]  = t[3];
      r->vec[i].coeffs[j +  64]  = t[4];
      r->vec[i].coeffs[j +  80]  = t[5];
      r->vec[i].coeffs[j +  96]  = t[6];
      r->vec[i].coeffs[j + 112]  = t[7];
      r->vec[i].coeffs[j + 128]  = t[8];
      r->vec[i].coeffs[j + 144]  = t[9];
      r->vec[i].coeffs[j + 160]  = t[10];
      r->vec[i].coeffs[j + 176]  = t[11];
      r->vec[i].coeffs[j + 192]  = t[12];
      r->vec[i].coeffs[j + 208]  = t[13];
      r->vec[i].coeffs[j + 224]  = t[14];
      r->vec[i].coeffs[j + 240]  = t[15];
    }

    for(j = 0; j < OKAI_N; j++) {
      r->vec[i].coeffs[j] = ((r->vec[i].coeffs[j] * OKAI_Q) + 128) >> 8;
    }
    a += 256;
  }
}

void polyvec_decompresspk(polyvec *r, const unsigned char *a)
{
  int i;
  for(i = 0; i < OKAI_K; i++)
    poly_decompress2(&r->vec[i], a+288*i);
}

/*************************************************
* Name:        polyvec_tobytes
*
* Description: Serialize vector of polynomials
*
* Arguments:   - unsigned char *r: pointer to output byte array
*              - const polyvec *a: pointer to input vector of polynomials
**************************************************/
void polyvec_tobytes(unsigned char *r, const polyvec *a)
{
  unsigned int i;
  for(i=0;i<OKAI_K;i++)
    poly_tobytes(r+i*OKAI_POLYBYTES, &a->vec[i]);
}

/*************************************************
* Name:        polyvec_frombytes
*
* Description: De-serialize vector of polynomials;
*              inverse of polyvec_tobytes
*
* Arguments:   - unsigned char *r: pointer to output byte array
*              - const polyvec *a: pointer to input vector of polynomials
**************************************************/
void polyvec_frombytes(polyvec *r, const unsigned char *a)
{
  unsigned int i;
  for(i=0;i<OKAI_K;i++)
    poly_frombytes(&r->vec[i], a+i*OKAI_POLYBYTES);
}

/*************************************************
* Name:        polyvec_ntt
*
* Description: Apply forward NTT to all elements of a vector of polynomials
*
* Arguments:   - polyvec *r: pointer to in/output vector of polynomials
**************************************************/
void polyvec_ntt(polyvec *r)
{
  unsigned int i;
  for(i=0;i<OKAI_K;i++)
    poly_ntt(&r->vec[i]);
}

/*************************************************
* Name:        polyvec_invntt
*
* Description: Apply inverse NTT to all elements of a vector of polynomials
*
* Arguments:   - polyvec *r: pointer to in/output vector of polynomials
**************************************************/
void polyvec_invntt(polyvec *r)
{
  unsigned int i;
  for(i=0;i<OKAI_K;i++)
    poly_invntt(&r->vec[i]);
}

/*************************************************
* Name:        polyvec_pointwise_acc
*
* Description: Pointwise multiply elements of a and b and accumulate into r
*
* Arguments: - poly *r:          pointer to output polynomial
*            - const polyvec *a: pointer to first input vector of polynomials
*            - const polyvec *b: pointer to second input vector of polynomials
**************************************************/
void polyvec_pointwise_acc(poly *r, const polyvec *a, const polyvec *b)
{
  unsigned int i;
  poly t;

  poly_basemul(r, &a->vec[0], &b->vec[0]);
  for(i=1;i<OKAI_K;i++) {
    poly_basemul(&t, &a->vec[i], &b->vec[i]);
    poly_add(r, r, &t);
    poly_reduce(r);
  }
}

/*************************************************
* Name:        polyvec_reduce
*
* Description: Applies Barrett reduction to each coefficient
*              of each element of a vector of polynomials
*              for details of the Barrett reduction see comments in reduce.c
*
* Arguments:   - poly *r: pointer to input/output polynomial
**************************************************/
void polyvec_reduce(polyvec *r)
{
  unsigned int i;
  for(i=0;i<OKAI_K;i++)
    poly_reduce(&r->vec[i]);
}

/*************************************************
* Name:        polyvec_add
*
* Description: Add vectors of polynomials
*
* Arguments: - polyvec *r:       pointer to output vector of polynomials
*            - const polyvec *a: pointer to first input vector of polynomials
*            - const polyvec *b: pointer to second input vector of polynomials
**************************************************/
void polyvec_add(polyvec *r, const polyvec *a, const polyvec *b)
{
  unsigned int i;
  for(i=0;i<OKAI_K;i++)
    poly_add(&r->vec[i], &a->vec[i], &b->vec[i]);
}
