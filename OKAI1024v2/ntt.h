#ifndef NTT_H
#define NTT_H

#include <stdint.h>

extern int16_t zetas[64];
extern int16_t zetasinv[128];
extern const int16_t tntt_zetas[128];
void hntt_ntt(int16_t r[512]);
void hntt_invntt(int16_t r[512]);
// void basemul(int16_t r[2], const int16_t a[2], const int16_t b[2], int16_t zeta);

#endif
