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
void polyvec_compress(unsigned char *r, const polyvec *a)
{
  int i,j,k;
  uint16_t t[4];
  for(i=0;i<OKAI_K;i++)
  {
    for(j=0;j<OKAI_N/4;j++)
    {
      for(k=0;k<4;k++)
        t[k] = ((((uint32_t)freeze(a->vec[i].coeffs[4*j+k]) << 10) + OKAI_Q/2) / OKAI_Q) & 0x3ff;

      r[5*j+ 0] =  t[0] & 0xff;
      r[5*j+ 1] = (t[0] >> 8) | ((t[1] & 0x3f) << 2);
      r[5*j+ 2] = (t[1] >> 6) | ((t[2] & 0x0f) << 4);
      r[5*j+ 3] = (t[2] >> 4) | ((t[3] & 0x03) << 6);
      r[5*j+ 4] = (t[3] >> 2);
    }
    r += 320*2;
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
void polyvec_decompress(polyvec *r, const unsigned char *a)
{
  int i,j;

  for(i=0;i<OKAI_K;i++)
  {
    for(j=0;j<OKAI_N/4;j++)
    {
      r->vec[i].coeffs[4*j+0] =  (((a[5*j+ 0]       | (((uint32_t)a[5*j+ 1] & 0x03) << 8)) * OKAI_Q) + 512) >> 10;
      r->vec[i].coeffs[4*j+1] = ((((a[5*j+ 1] >> 2) | (((uint32_t)a[5*j+ 2] & 0x0f) << 6)) * OKAI_Q) + 512) >> 10;
      r->vec[i].coeffs[4*j+2] = ((((a[5*j+ 2] >> 4) | (((uint32_t)a[5*j+ 3] & 0x3f) << 4)) * OKAI_Q) + 512) >> 10;
      r->vec[i].coeffs[4*j+3] = ((((a[5*j+ 3] >> 6) | (((uint32_t)a[5*j+ 4] & 0xff) << 2)) * OKAI_Q) + 512) >> 10;
    }
    a += 320*2;
  }
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
  poly t;

  poly_basemul(r, &a->vec[0], &b->vec[0]);
  poly_basemul((poly *)((int16_t *)(r)+256), (poly *)((int16_t *)(a->vec) + 256), (poly *)((int16_t *)(b->vec) + 256));
  poly_basemul(&t, &a->vec[1], &b->vec[1]);
  poly_basemul((poly *)((int16_t *)(&t)+256), (poly *)((int16_t *)(a->vec + 1) + 256), (poly *)((int16_t *)(b->vec + 1) + 256));
  poly_add(r, r, &t);
  poly_reduce(r);
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
