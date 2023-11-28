#ifndef CBD_H
#define CBD_H

#include <stdint.h>

#include "poly.h"

/* avx2 */
void cbd1(poly *r, const uint8_t *buf);
void cbd4(poly *r, const unsigned char *buf);
void cbd6(poly *r, const unsigned char *buf);
#endif
