#include <stdint.h>
#include <stddef.h>
#include "poly.h"
#include "ntt.h"
#include "polyvec.h"
#include "reduce.h"
#include "cbd.h"
#include "fips202.h"
int16_t p1_y[256]={
  0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1
};

void poly_compress_10(uint8_t r[320], const int16_t *acoeffs) {
    unsigned int i, j;
    uint16_t t[16];
    for (i = 0; i < 16; i++) {
        for (j = 0; j < 16; j++) {
            t[j] = acoeffs[16 * j + i];
            t[j] += ((int16_t) t[j] >> 15) & OKAI_Q;
            t[j] = ((((uint32_t) t[j] << 10) + OKAI_Q / 2) / OKAI_Q) & 0x3ff;
        }
        r[2*i +   0] = (t[0] >> 0);
        r[2*i +   1] = (t[0] >> 8) | ((t[1] & 0x3f) << 2);
        r[2*i +  32] = (t[1] >> 6) | ((t[2] & 0x0f) << 4);
        r[2*i +  33] = (t[2] >> 4) | ((t[3] & 0x03) << 6);
        r[2*i +  64] = (t[3] >> 2);

        r[2*i +  65] = (t[4] >>  0);
        r[2*i +  96] = (t[4] >> 8) | ((t[5] & 0x3f) << 2);
        r[2*i +  97] = (t[5] >> 6) | ((t[6] & 0x0f) << 4);
        r[2*i + 128] = (t[6] >> 4) | ((t[7] & 0x03) << 6);
        r[2*i + 129] = (t[7] >> 2);

        r[2*i + 160] = (t[8]  >> 0);
        r[2*i + 161] = (t[8]  >> 8) | ((t[9]  & 0x3f) << 2);
        r[2*i + 192] = (t[9]  >> 6) | ((t[10] & 0x0f) << 4);
        r[2*i + 193] = (t[10] >> 4) | ((t[11] & 0x03) << 6);
        r[2*i + 224] = (t[11] >> 2);

        r[2*i + 225] = (t[12] >> 0);
        r[2*i + 256] = (t[12] >> 8) | ((t[13] & 0x3f) << 2);
        r[2*i + 257] = (t[13] >> 6) | ((t[14] & 0x0f) << 4);
        r[2*i + 288] = (t[14] >> 4) | ((t[15] & 0x03) << 6);
        r[2*i + 289] = (t[15] >> 2);
    }
}

void poly_decompress_10(int16_t *rcoeffs, const uint8_t a[320]) {
    unsigned int i, j;
    uint16_t t[20];

    for (i = 0; i < 16; i++) {
        for(j = 0; j < 10; j++) {
            t[2*j]   = a[2*i + 32*j];
            t[2*j+1] = a[2*i + 32*j + 1];
        }
        rcoeffs[i +   0]  = t[0];
        rcoeffs[i +   0] += ((int16_t)t[1] & 0x03) << 8;
        rcoeffs[i +  16]  = t[1] >> 2;
        rcoeffs[i +  16] += ((int16_t)t[2] & 0x0f) << 6;
        rcoeffs[i +  32]  = t[2] >> 4;
        rcoeffs[i +  32] += ((int16_t)t[3] & 0x3f) << 4;
        rcoeffs[i +  48]  = t[3] >> 6;
        rcoeffs[i +  48] += ((int16_t)t[4] & 0xff) << 2;

        rcoeffs[i +  64]  = t[5];
        rcoeffs[i +  64] += ((int16_t)t[6] & 0x03) << 8;
        rcoeffs[i +  80]  = t[6] >> 2;
        rcoeffs[i +  80] += ((int16_t)t[7] & 0x0f) << 6;
        rcoeffs[i +  96]  = t[7] >> 4;
        rcoeffs[i +  96] += ((int16_t)t[8] & 0x3f) << 4;
        rcoeffs[i + 112]  = t[8] >> 6;
        rcoeffs[i + 112] += ((int16_t)t[9] & 0xff) << 2;

        rcoeffs[i + 128]  = t[10];
        rcoeffs[i + 128] += ((int16_t)t[11] & 0x03) << 8;
        rcoeffs[i + 144]  = t[11] >> 2;
        rcoeffs[i + 144] += ((int16_t)t[12] & 0x0f) << 6;
        rcoeffs[i + 160]  = t[12] >> 4;
        rcoeffs[i + 160] += ((int16_t)t[13] & 0x3f) << 4;
        rcoeffs[i + 176]  = t[13] >> 6;
        rcoeffs[i + 176] += ((int16_t)t[14] & 0xff) << 2;

        rcoeffs[i + 192]  = t[15];
        rcoeffs[i + 192] += ((int16_t)t[16] & 0x03) << 8;
        rcoeffs[i + 208]  = t[16] >> 2;
        rcoeffs[i + 208] += ((int16_t)t[17] & 0x0f) << 6;
        rcoeffs[i + 224]  = t[17] >> 4;
        rcoeffs[i + 224] += ((int16_t)t[18] & 0x3f) << 4;
        rcoeffs[i + 240]  = t[18] >> 6;
        rcoeffs[i + 240] += ((int16_t)t[19] & 0xff) << 2;
    }
    for (i = 0; i < 256; i++)
        rcoeffs[i] = ((uint32_t) (rcoeffs[i] & 0x3FF) * OKAI_Q + 512) >> 10;
}

void poly_compress512(uint8_t r[320*2], const poly *a) {
    poly_compress_10(r, &a->coeffs[0]);
    poly_compress_10(r+320, &a->coeffs[256]);
}

void poly_decompress512(poly *r, const uint8_t a[320*2]){
    poly_decompress_10(&r->coeffs[0], a);
    poly_decompress_10(&r->coeffs[256], a + 320);
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

void poly_con(uint8_t *c,
              const uint8_t *m,
              poly *s)
{
  poly_con_4(c, m, &s->coeffs[0]);
  poly_con_4(&c[128], &m[32], &s->coeffs[256]);
}

void poly_rec(uint8_t *m,
              const uint8_t *c,
              poly *s) {
    poly_rec_4(m, c, &s->coeffs[0]);
    poly_rec_4(&m[32], &c[128], &s->coeffs[256]);
}

void poly_compress(unsigned char *r, poly *a)
{
  uint32_t t[8];
  unsigned int i,j,k=0;

  for(i=0;i<OKAI_N;i+=8)
  {
    for(j=0;j<8;j++)
      t[j] = (((freeze(a->coeffs[i+j]) << 3) + OKAI_Q/2)/OKAI_Q) & 7; // mod 2^3

    r[k]   =  t[0]       | (t[1] << 3) | (t[2] << 6);
    r[k+1] = (t[2] >> 2) | (t[3] << 1) | (t[4] << 4) | (t[5] << 7);
    r[k+2] = (t[5] >> 1) | (t[6] << 2) | (t[7] << 5);
    k += 3;
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
    r->coeffs[i+0] = (((a[0] & 7) * OKAI_Q) + 4)>> 3;
    r->coeffs[i+1] = ((((a[0] >> 3) & 7) * OKAI_Q) + 4)>> 3;
    r->coeffs[i+2] = ((((a[0] >> 6) | ((a[1] << 2) & 4)) * OKAI_Q) + 4)>> 3;
    r->coeffs[i+3] = ((((a[1] >> 1) & 7) * OKAI_Q) + 4)>> 3;
    r->coeffs[i+4] = ((((a[1] >> 4) & 7) * OKAI_Q) + 4)>> 3;
    r->coeffs[i+5] = ((((a[1] >> 7) | ((a[2] << 1) & 6)) * OKAI_Q) + 4)>> 3;
    r->coeffs[i+6] = ((((a[2] >> 2) & 7) * OKAI_Q) + 4)>> 3;
    r->coeffs[i+7] = ((((a[2] >> 5)) * OKAI_Q) + 4)>> 3;
    a += 3;
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
  int i,j;
  uint16_t t[8];

  for(i=0;i<OKAI_N/8;i++)
  {
    for(j=0;j<8;j++)
      t[j] = freeze(a->coeffs[8*i+j]);

    r[13*i+ 0] =  t[0]        & 0xff;
    r[13*i+ 1] = (t[0] >>  8) | ((t[1] & 0x07) << 5);
    r[13*i+ 2] = (t[1] >>  3) & 0xff;
    r[13*i+ 3] = (t[1] >> 11) | ((t[2] & 0x3f) << 2);
    r[13*i+ 4] = (t[2] >>  6) | ((t[3] & 0x01) << 7);
    r[13*i+ 5] = (t[3] >>  1) & 0xff;
    r[13*i+ 6] = (t[3] >>  9) | ((t[4] & 0x0f) << 4);
    r[13*i+ 7] = (t[4] >>  4) & 0xff;
    r[13*i+ 8] = (t[4] >> 12) | ((t[5] & 0x7f) << 1);
    r[13*i+ 9] = (t[5] >>  7) | ((t[6] & 0x03) << 6);
    r[13*i+10] = (t[6] >>  2) & 0xff;
    r[13*i+11] = (t[6] >> 10) | ((t[7] & 0x1f) << 3);
    r[13*i+12] = (t[7] >>  5);
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
  int i;
  for(i=0;i<OKAI_N/8;i++)
  {
    r->coeffs[8*i+0] =  a[13*i+ 0]       | (((uint16_t)a[13*i+ 1] & 0x1f) << 8);
    r->coeffs[8*i+1] = (a[13*i+ 1] >> 5) | (((uint16_t)a[13*i+ 2]       ) << 3) | (((uint16_t)a[13*i+ 3] & 0x03) << 11);
    r->coeffs[8*i+2] = (a[13*i+ 3] >> 2) | (((uint16_t)a[13*i+ 4] & 0x7f) << 6);
    r->coeffs[8*i+3] = (a[13*i+ 4] >> 7) | (((uint16_t)a[13*i+ 5]       ) << 1) | (((uint16_t)a[13*i+ 6] & 0x0f) <<  9);
    r->coeffs[8*i+4] = (a[13*i+ 6] >> 4) | (((uint16_t)a[13*i+ 7]       ) << 4) | (((uint16_t)a[13*i+ 8] & 0x01) << 12);
    r->coeffs[8*i+5] = (a[13*i+ 8] >> 1) | (((uint16_t)a[13*i+ 9] & 0x3f) << 7);
    r->coeffs[8*i+6] = (a[13*i+ 9] >> 6) | (((uint16_t)a[13*i+10]       ) << 2) | (((uint16_t)a[13*i+11] & 0x07) << 10);
    r->coeffs[8*i+7] = (a[13*i+11] >> 3) | (((uint16_t)a[13*i+12]       ) << 5);
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
  cbd6(r, buf);
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
  hntt_ntt(r->coeffs);
  // poly_reduce(r);
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
  hntt_invntt(r->coeffs);
  // poly_reduce(r);
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
  basemul_asm(r->coeffs + 256, a->coeffs + 256, b->coeffs + 256, zetas);
}

extern void basemul_asm_acc(int16_t *, const int16_t *, const int16_t *, const int16_t *);
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
  basemul_asm_acc(r->coeffs + 256, a->coeffs + 256, b->coeffs + 256, zetas);
}

void basemul(int16_t r[2], const int16_t a[2], const int16_t b[2], int16_t zeta) {
    int16_t tmp1, tmp2;
    tmp1 = a[0] + a[1];
    tmp2 = b[0] + b[1];
    r[1] = montgomery_reduce((int32_t)tmp1 * tmp2); 
    tmp1 = montgomery_reduce((int32_t)a[1] * b[1]);
    tmp2 = montgomery_reduce((int32_t)a[0] * b[0]);
    r[1] = r[1] - tmp1 - tmp2;
    r[0]  = montgomery_reduce((int32_t)tmp1 * zeta);
    r[0] += tmp2;
}

/* yhat:even->0 odd->1 */
/* (0+x)*(a0+a1x) -> x*(a0+a1x) -> a0x+a1x^2 -> a1*(±zetas)+a0x*/
void basemulyhat(int16_t r[2], const int16_t a[2], int16_t zeta) {
    r[1] = a[0];
    r[0] = montgomery_reduce((int32_t)a[1] * zeta);//a1*(±zetas)
}

void poly_pointwise(int16_t r[256], const int16_t a[256], const int16_t b[256])
{
    unsigned int i;
    for(i=0;i<64;i++) {
        basemul(r+4*i, a+4*i, b+4*i, zetas[i]);
        basemul(r+4*i+2, a+4*i+2, b+4*i+2, -zetas[i]);
    }
}

void poly_pointwiseyhat(int16_t r[256], const int16_t a[256])
{
    unsigned int i;
    for(i=0;i<64;i++) {
        basemulyhat(r+4*i, a+4*i, zetas[i]);
        basemulyhat(r+4*i+2, a+4*i+2, -zetas[i]);
    }
}

extern void multi_basemul_asm(int16_t *, const int16_t *, const int16_t *, const int16_t *);
void poly_multi_basemul(poly *r, const poly *a, const poly *b) {
    int16_t p1_y[256];
    int16_t c0[256], c1[256], c3[256];
    int16_t a2[256], b2[256];
    basemul_asm(c0,a->coeffs,b->coeffs,zetas);
    basemul_asm(c1,a->coeffs+256,b->coeffs+256,zetas);
    pointwise_add(a2,a->coeffs,a->coeffs+256);
    pointwise_add(b2,b->coeffs,b->coeffs+256);
    basemul_asm(c3,a2,b2,zetas);
    poly_pointwiseyhat(p1_y,c1);
    pointwise_sub(r->coeffs+256,c3,c0);
    pointwise_sub(r->coeffs+256,r->coeffs+256,c1);
    pointwise_add(r->coeffs,c0,p1_y);
    poly_reduce(r);
}

void poly_single_basemul(int16_t *r,int16_t *a,int16_t *b,int16_t zeta){
  int i;
  int16_t c0[2],c1[2],c3[2],a2[2],b2[2];
  int16_t p1_y[2]={0,1};
  basemul(c0, a, b, zeta);
  basemul(c1, a+256, b+2, zeta);
  for(i = 0; i < 2; i ++)
  {
      a2[i] = (a[i] + a[256+i]);//a0+a1->a2
      b2[i] = (b[i] + b[256+i]);//b0+b1->b2
  }
  basemul(c3,a2,b2,zeta);//poly_pointwise(a0+a1,b0+b1)->c3
  basemulyhat(p1_y,c1,zeta);//yhat pointwise p1hat
  for(i = 0; i < 2; i ++)
  {
      r[256+i] += barrett_reduce(c3[i] - c0[i] - c1[i]);//phat-p0hat-p1hat h1hat
      r[i] += c0[i] + p1_y[i];//p0hat+yhat pointwise p1hat h0hat
  }
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
  asm_barrett_reduce(r->coeffs + 256);
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
  pointwise_add(r->coeffs, a->coeffs, b->coeffs);
  pointwise_add(r->coeffs + 256, a->coeffs + 256, b->coeffs + 256);
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
  pointwise_sub(r->coeffs, a->coeffs, b->coeffs);
  pointwise_sub(r->coeffs + 256, a->coeffs + 256, b->coeffs + 256);
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

/*************************************************
* Name:        poly_compress_5
*
* Description: Compress and serialize polynomial
*
* Arguments:   - unsigned char *r: pointer to output byte array
*              - const polyvec *a: pointer to input polynomial
**************************************************/
void poly_compress_5(unsigned char *r, poly *a)
{
  uint16_t t[4];
  int j, k;
  for (j = 0; j < OKAI_N / 4; j++)
  {
    for (k = 0; k < 4; k++)
      t[k] = ((((uint32_t)freeze(a->coeffs[4*j+k]) << 10) + OKAI_Q/2) / OKAI_Q) & 0x3ff;

    r[5*j+ 0] =  t[0] & 0xff;
    r[5*j+ 1] = (t[0] >> 8) | ((t[1] & 0x3f) << 2);
    r[5*j+ 2] = (t[1] >> 6) | ((t[2] & 0x0f) << 4);
    r[5*j+ 3] = (t[2] >> 4) | ((t[3] & 0x03) << 6);
    r[5*j+ 4] = (t[3] >> 2);
  }
}

void poly_decompress_5(poly *r, unsigned char *a)
{
  int j;
  for (j = 0; j < OKAI_N / 4; j++)
  {
    r->coeffs[4*j+0] =  (((a[5*j+ 0]       | (((uint32_t)a[5*j+ 1] & 0x03) << 8)) * OKAI_Q) + 512) >> 10;
    r->coeffs[4*j+1] = ((((a[5*j+ 1] >> 2) | (((uint32_t)a[5*j+ 2] & 0x0f) << 6)) * OKAI_Q) + 512) >> 10;
    r->coeffs[4*j+2] = ((((a[5*j+ 2] >> 4) | (((uint32_t)a[5*j+ 3] & 0x3f) << 4)) * OKAI_Q) + 512) >> 10;
    r->coeffs[4*j+3] = ((((a[5*j+ 3] >> 6) | (((uint32_t)a[5*j+ 4] & 0xff) << 2)) * OKAI_Q) + 512) >> 10;
  }
}

/*************************************************
* Name:        poly_zeroize
*
* Description: Zeros a polynomial
*
* Arguments:   - poly *p: pointer to polynomial
**************************************************/
void poly_zeroize(poly *p)
{
  int i;
  for (i = 0; i < OKAI_N; i++)
    p->coeffs[i] = 0;
}