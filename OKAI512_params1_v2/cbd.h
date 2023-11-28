#ifndef CBD_H
#define CBD_H

#include <stdint.h>
#include "poly.h"

// void cbd(poly *r, const unsigned char *buf);
// void cbd4(poly *r, const unsigned char *buf);

void cbd2(poly *r, const uint8_t buf[2 * OKAI_N / 4]);
void cbd12(poly *r, const unsigned char *buf);

#endif
