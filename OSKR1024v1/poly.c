#include "poly.h"

#include "cbd.h"
#include "ntt.h"
#include "params.h"
#include "symmetric.h"

#include <stdint.h>

extern void con5(uint8_t *, uint32_t *, int32_t *);
static void con5S(uint8_t c[160], const uint8_t m[32], int16_t coeffs[256]){
  size_t i,j;
  int32_t sp[32];
  uint32_t *pm = (uint32_t *)m;
  for(i=0;i<8;i++) {
    for(j=0;j<32;j++) {
      sp[j] = *(coeffs + (8*j+i));
    }
    con5(c, pm[i], sp);
    c += 4;
  }
}

void poly_con(uint8_t c[OSKR_POLYCOMPRESSEDBYTES],
              const uint8_t m[OSKR_INDCPA_MSGBYTES],
              poly *s)
{
  con5S(&c[0], &m[0], &s->coeffs[0]);
  con5S(&c[160], &m[32], &s->coeffs[256]);
}

int cmp_poly_con_part(uint8_t *c, const uint8_t *m, int16_t *coeffs)
{
  size_t i,j;
  int k=0;
  uint8_t ctmp[160];
  int16_t sp[32];
  uint32_t *pm = (uint32_t *)m;
  unsigned char rc=0;

  for(i=0;i<8;i++) {
    for(j=0;j<32;j++) {
      sp[j] = coeffs[8*j+i];
    }
    con5(ctmp+4*i, &pm[i], sp);
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

    rc |= c[k+128] ^ ctmp[k+128];
    rc |= c[k+129] ^ ctmp[k+129];
    rc |= c[k+130] ^ ctmp[k+130];
    rc |= c[k+131] ^ ctmp[k+131];
    k += 4;
    //ctmp+=4;
  }
  return rc;
}

int cmp_poly_con(const uint8_t c[OSKR_POLYCOMPRESSEDBYTES],
              const uint8_t m[OSKR_INDCPA_MSGBYTES],
              poly *s)
{
  unsigned char rc=0;
  rc |= cmp_poly_con_part(&c[0], &m[0], &s->coeffs[0]);
  rc |= cmp_poly_con_part(&c[160], &m[32], &s->coeffs[256]);
  return rc;
}

extern void rec5(uint32_t *,const uint8_t *, int16_t *);
static void rec5S(uint8_t *m,const uint8_t *c,int16_t *s)
{
  int16_t sp[32];
  unsigned int i,j;
  uint32_t *pm = (uint32_t *)m;
  for(i=0;i<256/32;i++) 
  {
    for(j=0;j<32;j++)
    {
      sp[j] = s[8*j+i];
    }
    rec5(&pm[i], c, sp);
    c += 4;
  }
}
void poly_rec(uint8_t m[OSKR_INDCPA_MSGBYTES],
              const uint8_t c[OSKR_POLYCOMPRESSEDBYTES],
              poly *s)
{
  rec5S(&m[0], &c[0], &s->coeffs[0]);
  rec5S(&m[32], &c[160], &s->coeffs[256]);
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
  uint16_t t[8];
  for(j=0;j<OSKR_N/8;j++) {
    for(k=0;k<8;k++)
      t[k] = ((((uint32_t)a->coeffs[8*j+k] << 11) + OSKR_Q/2) / OSKR_Q) & 0x7ff;

    r[352*2*i+11*j+ 0] =  t[0] & 0xff;
    r[352*2*i+11*j+ 1] = (t[0] >>  8) | ((t[1] & 0x1f) << 3);
    r[352*2*i+11*j+ 2] = (t[1] >>  5) | ((t[2] & 0x03) << 6);
    r[352*2*i+11*j+ 3] = (t[2] >>  2) & 0xff;
    r[352*2*i+11*j+ 4] = (t[2] >> 10) | ((t[3] & 0x7f) << 1);
    r[352*2*i+11*j+ 5] = (t[3] >>  7) | ((t[4] & 0x0f) << 4);
    r[352*2*i+11*j+ 6] = (t[4] >>  4) | ((t[5] & 0x01) << 7);
    r[352*2*i+11*j+ 7] = (t[5] >>  1) & 0xff;
    r[352*2*i+11*j+ 8] = (t[5] >>  9) | ((t[6] & 0x3f) << 2);
    r[352*2*i+11*j+ 9] = (t[6] >>  6) | ((t[7] & 0x07) << 5);
    r[352*2*i+11*j+10] = (t[7] >>  3);
  }
}

int cmp_poly_packcompress(const unsigned char *r, poly *a, int i) {
    unsigned char rc = 0;
    int j, k;
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
    return rc;
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
  for(j=0;j<OSKR_N/8;j++)
  {
    r->coeffs[8*j+0] =  (((a[352*2*i+11*j+ 0]       | (((uint32_t)a[352*2*i+11*j+ 1] & 0x07) << 8)) * OSKR_Q) + 1024) >> 11;
    r->coeffs[8*j+1] = ((((a[352*2*i+11*j+ 1] >> 3) | (((uint32_t)a[352*2*i+11*j+ 2] & 0x3f) << 5)) * OSKR_Q) + 1024) >> 11;
    r->coeffs[8*j+2] = ((((a[352*2*i+11*j+ 2] >> 6) | (((uint32_t)a[352*2*i+11*j+ 3] & 0xff) << 2) | (((uint32_t)a[352*2*i+11*j+4] & 0x01) << 10)) * OSKR_Q) + 1024) >> 11;
    r->coeffs[8*j+3] = ((((a[352*2*i+11*j+ 4] >> 1) | (((uint32_t)a[352*2*i+11*j+ 5] & 0x0f) << 7)) * OSKR_Q) + 1024) >> 11;
    r->coeffs[8*j+4] = ((((a[352*2*i+11*j+ 5] >> 4) | (((uint32_t)a[352*2*i+11*j+ 6] & 0x7f) << 4)) * OSKR_Q) + 1024) >> 11;
    r->coeffs[8*j+5] = ((((a[352*2*i+11*j+ 6] >> 7) | (((uint32_t)a[352*2*i+11*j+ 7] & 0xff) << 1) | (((uint32_t)a[352*2*i+11*j+8] & 0x03) <<  9)) * OSKR_Q) + 1024) >> 11;
    r->coeffs[8*j+6] = ((((a[352*2*i+11*j+ 8] >> 2) | (((uint32_t)a[352*2*i+11*j+ 9] & 0x1f) << 6)) * OSKR_Q) + 1024) >> 11;
    r->coeffs[8*j+7] = ((((a[352*2*i+11*j+ 9] >> 5) | (((uint32_t)a[352*2*i+11*j+10] & 0xff) << 3)) * OSKR_Q) + 1024) >> 11;
  }
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
* Name:        poly_frombytes_mul
*
* Description: Multiplication of a polynomial with a de-serialization of another polynomial
*
* Arguments:   - poly *r:                pointer to output polynomial
*              - const poly *b:          pointer to input polynomial
*              - const unsigned char *a: pointer to input byte array (of OSKR_POLYBYTES bytes)
**************************************************/
extern void frombytes_mul_asm(int16_t *r, const int16_t *b, const unsigned char *c, const int32_t zetas[64]);
void poly_frombytes_mul(poly *r, const poly *b, const unsigned char *a) {
    frombytes_mul_asm(r->coeffs, b->coeffs, a, zetas);
    frombytes_mul_asm(r->coeffs+256, b->coeffs+256, a+384, zetas);
}

/*************************************************
* Name:        poly_frombytes_mul_acc
*
* Description: Multiplication of a polynomial with a de-serialization of another polynomial
*              Accumulation in r.
*
* Arguments:   - poly *r:                pointer to output polynomial
*              - const poly *b:          pointer to input polynomial
*              - const unsigned char *a: pointer to input byte array (of OSKR_POLYBYTES bytes)
**************************************************/
extern void frombytes_mul_asm_acc(int16_t *r, const int16_t *b, const unsigned char *c, const int32_t zetas[64]);
void poly_frombytes_mul_acc(poly *r, const poly *b, const unsigned char *a) {
    frombytes_mul_asm_acc(r->coeffs, b->coeffs, a, zetas);
    frombytes_mul_asm_acc(r->coeffs+256, b->coeffs+256, a+384, zetas);
}

/*************************************************
* Name:        poly_getnoise
*
* Description: Sample a polynomial deterministically from a seed and a nonce,
*              with output polynomial close to centered binomial distribution
*              with parameter OSKR_ETA
*
* Arguments:   - poly *r:                   pointer to output polynomial
*              - const unsigned char *seed: pointer to input seed (pointing to array of length OSKR_SYMBYTES bytes)
*              - unsigned char nonce:       one-byte input nonce
*              - int add:                   boolean to indicate to accumulate into r
**************************************************/
void poly_noise(poly *r, const unsigned char *seed, unsigned char nonce, int add) {
    unsigned char buf[OSKR_ETA * OSKR_N / 4];

    prf(buf, OSKR_ETA * OSKR_N / 4, seed, nonce);
    cbd(r, buf, add);
}

extern void basemul_asm(int16_t *, const int16_t *, const int16_t *, const int32_t *);
extern void multi_basemul_asm(int16_t *, const int16_t *, const int16_t *, const int32_t *);
extern void combine_asm(int16_t *, int16_t *, int16_t *);
/*************************************************
* Name:        poly_basemul
*
* Description: Multiplication of two polynomials in NTT domain
*
* Arguments:   - poly *r:       pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_basemul(poly *r, const poly *a, const poly *b) {
    basemul_asm(r->coeffs, a->coeffs, b->coeffs, zetas);
    basemul_asm(r->coeffs+256, a->coeffs+256, b->coeffs+256, zetas);
}

void poly_multi_basemul(poly *r, const poly *a, const poly *b) {
    multi_basemul_asm(r->coeffs,a->coeffs,b->coeffs,zetas);
}

extern void basemul_asm_acc(int16_t *, const int16_t *, const int16_t *, const int32_t *);
/*************************************************
* Name:        poly_basemul_acc
*
* Description: Multiplication of two polynomials in NTT domain, accumulating
*
* Arguments:   - poly *r:       pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_basemul_acc(poly *r, const poly *a, const poly *b) {
    basemul_asm_acc(r->coeffs, a->coeffs, b->coeffs, zetas);
    basemul_asm_acc(r->coeffs+256, a->coeffs+256, b->coeffs+256, zetas);
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
* Name:        poly_frommont
*
* Description: Inplace conversion of all coefficients of a polynomial
*              from Montgomery domain to normal domain
*
* Arguments:   - poly *r:       pointer to input/output polynomial
**************************************************/
void poly_frommont(poly *r) {
  asm_fromplant(r->coeffs);
  asm_fromplant(r->coeffs+256);
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
  asm_barrett_reduce(r->coeffs+256);
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
void poly_add(poly *r, poly *a, poly *b) {
    pointwise_add(r->coeffs,a->coeffs,b->coeffs);
    pointwise_add(r->coeffs+256,a->coeffs+256,b->coeffs+256);
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
    pointwise_sub(r->coeffs+256,a->coeffs+256,b->coeffs+256);
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
