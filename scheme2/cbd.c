#include <immintrin.h>
#include <stdint.h>

#include "cbd.h"
#include "params.h"

void cbd1(poly *r, const uint8_t *buf)
{
  unsigned int i;
  __m256i vec0, vec1, vec2;
  const __m256i mask = _mm256_set1_epi16(1);

  vec0 = _mm256_load_si256((__m256i *)&buf[0]);
  for(i=0;i<16;i+=2)
  {
    vec1 = _mm256_srli_epi32(vec0, i);
    vec1 = _mm256_and_si256(vec1, mask);
    vec2 = _mm256_srli_epi32(vec0, i+1);
    vec2 = _mm256_and_si256(vec2, mask);
    vec1 = _mm256_sub_epi16(vec1, vec2);
    _mm256_store_si256((__m256i *)&r->coeffs[8*i], vec1);
  }
  vec0 = _mm256_load_si256((__m256i *)&buf[32]);
  for(i=0;i<16;i+=2)
  {
    vec1 = _mm256_srli_epi32(vec0, i);
    vec1 = _mm256_and_si256(vec1, mask);
    vec2 = _mm256_srli_epi32(vec0, i+1);
    vec2 = _mm256_and_si256(vec2, mask);
    vec1 = _mm256_sub_epi16(vec1, vec2);
    _mm256_store_si256((__m256i *)&r->coeffs[128+8*i], vec1);
  }
#if (OKAI_N == 512)
  vec0 = _mm256_load_si256((__m256i *)&buf[64]);
  for(i=0;i<16;i+=2)
  {
    vec1 = _mm256_srli_epi32(vec0, i);
    vec1 = _mm256_and_si256(vec1, mask);
    vec2 = _mm256_srli_epi32(vec0, i+1);
    vec2 = _mm256_and_si256(vec2, mask);
    vec1 = _mm256_sub_epi16(vec1, vec2);
    _mm256_store_si256((__m256i *)&r->coeffs[256+8*i], vec1);
  }
  vec0 = _mm256_load_si256((__m256i *)&buf[96]);
  for(i=0;i<16;i+=2)
  {
    vec1 = _mm256_srli_epi32(vec0, i);
    vec1 = _mm256_and_si256(vec1, mask);
    vec2 = _mm256_srli_epi32(vec0, i+1);
    vec2 = _mm256_and_si256(vec2, mask);
    vec1 = _mm256_sub_epi16(vec1, vec2);
    _mm256_store_si256((__m256i *)&r->coeffs[384+8*i], vec1);
  }
#endif
}

void cbd4(poly *r, const unsigned char *buf)
{
	int i;
	__m256i * pbuf = (__m256i *) buf;
  __m256i * pr = (__m256i *) r->coeffs;
     
	__m256i mask55 = _mm256_set1_epi8(0x55);
  __m256i mask33 = _mm256_set1_epi8(0x33);
  __m256i mask0f = _mm256_set1_epi8(0x0f);
  __m256i q16x = _mm256_set1_epi16(OKAI_Q);
     
  __m256i  t,d;
  __m128i *pd = (__m128i *) &d;
     
  for(i=0;i<OKAI_N/32;i++)
	{
    d = _mm256_and_si256(pbuf[i],mask55);
    t = _mm256_srli_epi16(pbuf[i],1);
    t = _mm256_and_si256(t,mask55);
    t = _mm256_add_epi8(d,t);
          
    d = _mm256_and_si256(t,mask33);
    t = _mm256_srli_epi16(t,2);
    t = _mm256_and_si256(t,mask33);
    d = _mm256_add_epi8(d,t);

    t = _mm256_and_si256(d,mask0f);
    d = _mm256_srli_epi16(d,4);
    d = _mm256_and_si256(d,mask0f);
          
    d = _mm256_sub_epi8(t,d);
    t = _mm256_cvtepi8_epi16(*pd);
    pr[2*i] = _mm256_add_epi16(t,q16x);
          
    d = _mm256_permute2x128_si256(d,d,0x21);
    t = _mm256_cvtepi8_epi16(*pd);
    pr[2*i+1] = _mm256_add_epi16(t,q16x);
	}
}

void cbd6(poly *r, const unsigned char *buf)
{
  int i;
  const __m256i maskand = _mm256_set1_epi32(0x00041041);
  const __m256i mask3f = _mm256_set1_epi32(0x3f);
  const __m256i idx8  = _mm256_set_epi8(-1,15,14,13,-1,12,11,10,
                                        -1, 9, 8, 7,-1, 6, 5, 4,
                                        -1,11,10, 9,-1, 8, 7, 6,
                                        -1, 5, 4, 3,-1, 2, 1, 0);
  __m256i  t,d,shiftres;
  for(i=0;i<OKAI_N/16;i++){
    d = _mm256_loadu_si256((__m256i *)&buf[24*i]);
    d = _mm256_permute4x64_epi64(d, 0x94);
    d = _mm256_shuffle_epi8(d, idx8);// load 24b

    t = _mm256_and_si256(d,maskand);// t & 0x00041041
    
    shiftres = _mm256_srli_epi32(d,1);// t>>1
    shiftres = _mm256_and_si256(shiftres,maskand);
    t = _mm256_add_epi32(shiftres,t);

    shiftres = _mm256_srli_epi32(d,2);// t>>2
    shiftres = _mm256_and_si256(shiftres,maskand);
    t = _mm256_add_epi32(shiftres,t);

    shiftres = _mm256_srli_epi32(d,3);// t>>3
    shiftres = _mm256_and_si256(shiftres,maskand);
    t = _mm256_add_epi32(shiftres,t);

    shiftres = _mm256_srli_epi32(d,4);// t>>4
    shiftres = _mm256_and_si256(shiftres,maskand);
    t = _mm256_add_epi32(shiftres,t);

    shiftres = _mm256_srli_epi32(d,5);// t>>5
    shiftres = _mm256_and_si256(shiftres,maskand);
    t = _mm256_add_epi32(shiftres,t);

    d = _mm256_and_si256(t,mask3f);// a
    shiftres = _mm256_and_si256(_mm256_srli_epi32(t,6),mask3f);// b
    shiftres = _mm256_sub_epi32(d,shiftres);// a-b

    d = _mm256_and_si256(_mm256_srli_epi32(t,12),mask3f);// a
    t = _mm256_and_si256(_mm256_srli_epi32(t,18),mask3f);// b
    t = _mm256_sub_epi32(d,t);// a-b

    t = _mm256_slli_epi32(t,16);
    shiftres = _mm256_or_si256(shiftres,t);

    _mm256_store_si256((__m256i *)&r->coeffs[16*i], shiftres);
  }
}