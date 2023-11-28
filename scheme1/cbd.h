#ifndef CBD_H
#define CBD_H

#include <stdint.h>
#include <immintrin.h>
#include "params.h"
#include "poly.h"

void poly_cbd_eta1(poly *r, const __m256i buf[OSKR_ETA1*OSKR_N/128+1]);
void poly_cbd_eta2(poly *r, const __m256i buf[OSKR_ETA2*OSKR_N/128]);

#endif
