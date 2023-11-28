#include "ntt.h"
#include "poly.h"
#include "polyvec.h"
#include "symmetric.h"
#include "matacc.h"

extern void multi_basemul_asm_acc_neg(int16_t *, int16_t *, int16_t *, const int32_t *);
extern void multi_basemul_asm_acc_pos(int16_t *, int16_t *, int16_t *, const int32_t *);
/*************************************************
* Name:        matacc
*
* Description: Multiplies a row of A or A^T, generated on-the-fly,
*              with a vector of polynomials and accumulates into the result.
*
* Arguments:   - poly *r:                    pointer to output polynomial to accumulate in
*              - polyvec *b:                 pointer to input vector of polynomials to multiply with
*              - unsigned char i:            byte to indicate the index < OSKR_K of the row of A or A^T
*              - const unsigned char *seed:  pointer to the public seed used to generate A
*              - int transposed:             boolean indicatin whether A or A^T is generated
**************************************************/
// void matacc(poly* r, const polyvec *b, unsigned char i, const unsigned char *seed, int transposed) {
//   unsigned char buf[XOF_BLOCKBYTES+2];
//   xof_state state;
//   int16_t c[4];
//   int j = 0;
  
//   if (transposed)
//       xof_absorb(&state, seed, i, j);
//     else
//       xof_absorb(&state, seed, j, i);

//     xof_squeezeblocks(buf, 1, &state);
//     matacc_asm(r->coeffs, b->vec[j].coeffs, c, buf, zetas, &state);
//     matacc_asm(r->coeffs+256, b->vec[j].coeffs+256, c, buf, zetas, &state);
//   for(j=1;j<OSKR_K;j++) {

//     if (transposed)
//       xof_absorb(&state, seed, i, j);
//     else
//       xof_absorb(&state, seed, j, i);

//     xof_squeezeblocks(buf, 1, &state);

//     matacc_asm_acc(r->coeffs, b->vec[j].coeffs, c, buf, zetas, &state);
//     matacc_asm_acc(r->coeffs+256, b->vec[j].coeffs+256, c, buf, zetas, &state);
//   } 
// }

void matacc(poly* r, polyvec *b, unsigned char i, const unsigned char *seed, int transposed) {
  unsigned char buf[XOF_BLOCKBYTES+1];
  xof_state state;
  int ctr, pos, k;
  uint16_t val;
  int16_t c[4];

  poly_zeroize(r);

  for(int j=0;j<OSKR_K;j++) {
    ctr = pos = 0;
    if (transposed)
      xof_absorb(&state, seed, i, j);
    else
      xof_absorb(&state, seed, j, i);

    xof_squeezeblocks(buf, 1, &state);

    while (ctr < OSKR_N/4)
    {
      k = 0;
      while(k < 4) {
        val = buf[pos] | ((uint16_t)buf[pos + 1] << 8);
        if (val < 19 * OSKR_Q) {
          val -= (val>>12)*OSKR_Q; // Barrett reduction
          c[k++] = (int16_t) val;
        }

        pos += 2;
        if (pos + 2 > XOF_BLOCKBYTES) {
          xof_squeezeblocks(buf, 1, &state);
          pos = 0;
        }
      }
      if(ctr%2){
        multi_basemul_asm_acc_neg(&r->coeffs[2*ctr],&b->vec[j].coeffs[2*ctr],c, zetas+ctr/2);
      }
      else{
        multi_basemul_asm_acc_pos(&r->coeffs[2*ctr],&b->vec[j].coeffs[2*ctr],c, zetas+ctr/2);
      }
      ctr++;
    }
  }
}

void gena(poly* r, unsigned char i, const unsigned char *seed, int transposed) {
  unsigned char buf[XOF_BLOCKBYTES+1];
  xof_state state;
  int ctr, pos, k;
  uint16_t val;

  poly_zeroize(r);

  for(int j=0;j<OSKR_K;j++) {
    ctr = pos = 0;
    if (transposed)
      xof_absorb(&state, seed, i, j);
    else
      xof_absorb(&state, seed, j, i);

    xof_squeezeblocks(buf, 1, &state);

    while (ctr < OSKR_N/4)
    {
      k = 0;
      while(k < 4) {
        val = buf[pos] | ((uint16_t)buf[pos + 1] << 8);
        if (val < 19 * OSKR_Q) {
          val -= (val>>12)*OSKR_Q; // Barrett reduction
          r->coeffs[4*ctr+k] = (int16_t) val;
          k++;
        }
        pos += 2;
        if (pos + 2 > XOF_BLOCKBYTES) {
          xof_squeezeblocks(buf, 1, &state);
          pos = 0;
        }
      }
      ctr++;
    }
  }
}