#include "cbd.h"
#include "params.h"

#include <stdint.h>

/*************************************************
* Name:        load32_littleendian
*
* Description: load bytes into a 32-bit integer
*              in little-endian order
*
* Arguments:   - const unsigned char *x: pointer to input byte array
*
* Returns 32-bit unsigned integer loaded from x
**************************************************/
static uint32_t load32_littleendian(const unsigned char *x) {
    uint32_t r;
    r  = (uint32_t)x[0];
    r |= (uint32_t)x[1] << 8;
    r |= (uint32_t)x[2] << 16;
    r |= (uint32_t)x[3] << 24;
    return r;
}

/*************************************************
* Name:        cbd
*
* Description: Given an array of uniformly random bytes, compute
*              polynomial with coefficients distributed according to
*              a centered binomial distribution with parameter OSKR_ETA
*              specialized for OSKR_ETA=2
*
* Arguments:   - poly *r:                  pointer to output polynomial
*              - const unsigned char *buf: pointer to input byte array
*              - int add:                  boolean to indicate to accumulate into r
**************************************************/
void cbd(poly *r, const unsigned char *buf, int add) {
    uint32_t d, t;
    int16_t a, b;
    int i, j;

    for (i = 0; i < OSKR_N / 8; i++) {
        t = load32_littleendian(buf + 4 * i);
        d  = t & 0x55555555;
        d += (t >> 1) & 0x55555555;

        for (j = 0; j < 8; j++) {
            a = (d >>  4 * j)    & 0x3;
            b = (d >> (4 * j + 2)) & 0x3;
            if (!add)
              r->coeffs[8 * i + j] = 0;
            r->coeffs[8 * i + j] = r->coeffs[8 * i + j] + (a - b);
        }
    }
}
