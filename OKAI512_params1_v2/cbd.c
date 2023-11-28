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
void cbd2(poly *r, const uint8_t buf[2 * OKAI_N / 4]) {
    unsigned int i, j;
    uint32_t t, d;
    int16_t a, b;

    for (i = 0; i < OKAI_N / 8; i++) {
        t = load32_littleendian(buf + 4 * i);
        d = t & 0x55555555;
        d += (t >> 1) & 0x55555555;

        for (j = 0; j < 8; j++) {
            a = (d >> (4 * j + 0)) & 0x3;
            b = (d >> (4 * j + 2)) & 0x3;
            r->coeffs[8 * i + j] = a - b;
        }
    }
}

void cbd12(poly *r, const unsigned char *buf)
{
  unsigned int i;
  uint32_t t,d;
  int16_t a,b;

  for(i=0;i<OKAI_N;i++)
  {
    t  = load24_littleendian(buf+3*i);
    d  = t & 0x00001001;
    d += (t>>1) & 0x00001001;
    d += (t>>2) & 0x00001001;
    d += (t>>3) & 0x00001001;
    d += (t>>4) & 0x00001001;
    d += (t>>5) & 0x00001001;
    d += (t>>6) & 0x00001001;
    d += (t>>7) & 0x00001001;
    d += (t>>8) & 0x00001001;
    d += (t>>9) & 0x00001001;
    d += (t>>10) & 0x00001001;
    d += (t>>11) & 0x00001001;

    a = d & 0xfff;// low 12b
    b = (d >> 12) & 0xfff;
    r->coeffs[i] = a - b;
  }
}
