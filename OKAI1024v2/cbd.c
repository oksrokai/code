#include "cbd.h"
#include "params.h"

#include <stddef.h>
#include <stdint.h>

/*************************************************
* Name:        load32_littleendian
*
* Description: load bytes into a 32-bit integer
*              in little-endian order
*
* Arguments:   - const uint8_t *x: pointer to input byte array
*
* Returns 32-bit unsigned integer loaded from x
**************************************************/
static uint32_t load32_littleendian(const uint8_t *x) {
    uint32_t r;
    r  = (uint32_t)x[0];
    r |= (uint32_t)x[1] << 8;
    r |= (uint32_t)x[2] << 16;
    r |= (uint32_t)x[3] << 24;
    return r;
}

static uint32_t load24_littleendian(const uint8_t x[3])
{
  uint32_t r;
  r  = (uint32_t)x[0];
  r |= (uint32_t)x[1] << 8;
  r |= (uint32_t)x[2] << 16;
  return r;
}

/*************************************************
* Name:        cbd
*
* Description: Given an array of uniformly random bytes, compute
*              polynomial with coefficients distributed according to
*              a centered binomial distribution with parameter OKAI_ETA
*              specialized for OKAI_ETA=2
*
* Arguments:   - poly *r:                  pointer to output polynomial
*              - const uint8_t *buf: pointer to input byte array
**************************************************/
void cbd(poly *r, const unsigned char *buf)
{
  uint32_t d,t;
  int16_t a,b,c,e;
  int i,j;

  for(i=0;i<OKAI_N/16;i++)
  {
    t = load32_littleendian(buf+4*i);
    d  = t & 0x55555555;
    d += (t>>1) & 0x55555555;

    for(j=0;j<16;j+=2)
    {
      a = (d >>  4*j)    & 0x1;
      b = (d >> (4*j+1)) & 0x1;
      c = (d >> (4*j+2)) & 0x1;
      e = (d >> (4*j+3)) & 0x1;
      r->coeffs[16*i+j] = a - b;
      r->coeffs[16*i+j+1] = c - e;
    }
  }
}

void cbd6(poly *r, const unsigned char *buf)
{
  unsigned int i,j;
  uint32_t t,d;
  int16_t a,b;

  for(i=0;i<OKAI_N/2;i++)
  {
    t  = load24_littleendian(buf+3*i);
    d  = t & 0x00041041;
    d += (t>>1) & 0x00041041;
    d += (t>>2) & 0x00041041;
    d += (t>>3) & 0x00041041;
    d += (t>>4) & 0x00041041;
    d += (t>>5) & 0x00041041;

    for(j=0;j<2;j++) {// 2 coeffs
      a = (d >> (12*j+0)) & 0x3f;// low 6b
      b = (d >> (12*j+6)) & 0x3f;
      r->coeffs[2*i+j] = a - b;
    }
  }
}
