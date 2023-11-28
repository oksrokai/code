#ifndef REJSAMPLE_H
#define REJSAMPLE_H

#include <stdint.h>

unsigned int rej_uniform_avx(int16_t *r,
                             const unsigned char *buf);

int rej_sample(uint16_t *r, size_t rlen, const unsigned char *buf, size_t buflen);
#endif
