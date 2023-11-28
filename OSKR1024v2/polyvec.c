#include <stdint.h>
#include "polyvec.h"
#include "poly.h"

/*************************************************
* Name:        polyvec_tobytes
*
* Description: Serialize vector of polynomials
*
* Arguments:   - unsigned char *r: pointer to output byte array (needs space for OSKR_POLYVECBYTES)
*              - const polyvec *a: pointer to input vector of polynomials
**************************************************/
void polyvec_tobytes(unsigned char *r, polyvec *a)
{
  int i;
  for(i=0;i<OSKR_K;i++)
    poly_tobytes(r+i*OSKR_POLYBYTES, &a->vec[i]);
}

/*************************************************
* Name:        polyvec_frombytes
*
* Description: De-serialize vector of polynomials;
*              inverse of polyvec_tobytes
*
* Arguments:   - unsigned char *r: pointer to output byte array
*              - const polyvec *a: pointer to input vector of polynomials (of length OSKR_POLYVECBYTES)
**************************************************/
void polyvec_frombytes(polyvec *r, const unsigned char *a)
{
  int i;
  for(i=0;i<OSKR_K;i++)
    poly_frombytes(&r->vec[i], a+i*OSKR_POLYBYTES);
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
  int i;
  for(i=0;i<OSKR_K;i++)
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
  int i;
  for(i=0;i<OSKR_K;i++)
    poly_invntt(&r->vec[i]);
}

/*************************************************
* Name:        polyvec_reduce
*
* Description: Applies Barrett reduction to each coefficient
*              of each element of a vector of polynomials
*              for details of the Barrett reduction see comments in reduce.c
*
* Arguments:   - poly *r:       pointer to input/output polynomial
**************************************************/
void polyvec_reduce(polyvec *r)
{
  int i;
  for(i=0;i<OSKR_K;i++)
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
  for(i=0;i<OSKR_K;i++)
    poly_add(&r->vec[i], &a->vec[i], &b->vec[i]);
}
