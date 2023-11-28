#include <stdint.h>
#include "polyvec.h"
#include "poly.h"
#include "ntt.h"
#include "reduce.h"

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
#if (OKAI_N == 512)
  poly_compress_to_10_avx2(r, &a->vec[0].coeffs[0]);
  poly_compress_to_10_avx2(r+320*1, &a->vec[0].coeffs[256]);
  poly_compress_to_10_avx2(r+320*2, &a->vec[1].coeffs[0]);
  poly_compress_to_10_avx2(r+320*3, &a->vec[1].coeffs[256]);
#else
  unsigned int i;
  for(i = 0; i < OKAI_K; i++)
    poly_compress_to_9_avx2(r+288*i, &a->vec[i]);
#endif
}

void polyvec_compressc(unsigned char *r, const polyvec *a)
{
#if (OKAI_N == 512)
  poly_compress_to_10_avx2(r, &a->vec[0].coeffs[0]);
  poly_compress_to_10_avx2(r+320*1, &a->vec[0].coeffs[256]);
  poly_compress_to_10_avx2(r+320*2, &a->vec[1].coeffs[0]);
  poly_compress_to_10_avx2(r+320*3, &a->vec[1].coeffs[256]);
#elif (OKAI_K == 3)
  unsigned int i;
  for(i = 0; i < OKAI_K; i++)
    poly_compress_to_9_avx2(r+288*i, &a->vec[i]);
#elif (OKAI_K == 2)
  poly_compress_to_8_avx2(r, &a->vec[0]);
  poly_compress_to_8_avx2(r+256, &a->vec[1]);
#endif
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
#if (OKAI_N == 512)
  poly_decompress_from_10_avx2(&r->vec[0].coeffs[0], a);
  poly_decompress_from_10_avx2(&r->vec[0].coeffs[256], a+320*1);
  poly_decompress_from_10_avx2(&r->vec[1].coeffs[0], a+320*2);
  poly_decompress_from_10_avx2(&r->vec[1].coeffs[256], a+320*3);
#else
  unsigned int i;
  for(i = 0; i < OKAI_K; i++)
    poly_decompress_from_9_avx2(&r->vec[i], a+288*i);
#endif
}

void polyvec_decompressc(polyvec *r, const unsigned char *a)
{
#if (OKAI_N == 512)
  poly_decompress_from_10_avx2(&r->vec[0].coeffs[0], a);
  poly_decompress_from_10_avx2(&r->vec[0].coeffs[256], a+320*1);
  poly_decompress_from_10_avx2(&r->vec[1].coeffs[0], a+320*2);
  poly_decompress_from_10_avx2(&r->vec[1].coeffs[256], a+320*3);
#elif (OKAI_K == 3)
  unsigned int i;
  for(i = 0; i < OKAI_K; i++)
    poly_decompress_from_9_avx2(&r->vec[i], a+288*i);
#elif (OKAI_K == 2)
  poly_decompress_from_8_avx2(&r->vec[0], a);
  poly_decompress_from_8_avx2(&r->vec[1], a+256);
#endif
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
#if (OKAI_N == 256)
  basemul_acc_avx(r->coeffs, a->vec->coeffs, b->vec->coeffs, zetas_exp + 152);
  basemul_acc_avx(r->coeffs + 64, a->vec->coeffs + 64, b->vec->coeffs + 64, zetas_exp + 184);
  basemul_acc_avx(r->coeffs + 128, a->vec->coeffs + 128, b->vec->coeffs + 128, zetas_exp + 348);
  basemul_acc_avx(r->coeffs + 192, a->vec->coeffs + 192, b->vec->coeffs + 192, zetas_exp + 380);
  poly_reduce(r);
#elif (OKAI_N == 512)
  poly t;
  poly_basemul_montgomery_multi(r, &a->vec[0], &b->vec[0]);
  poly_basemul_montgomery_multi(&t, &a->vec[1], &b->vec[1]);
  poly_add(r, r, &t);
  poly_reduce(r);
#endif
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
  int i;
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
  int i;
  for(i=0;i<OKAI_K;i++)
    poly_add(&r->vec[i], &a->vec[i], &b->vec[i]);
}
