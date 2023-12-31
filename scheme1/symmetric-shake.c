#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include "params.h"
#include "symmetric.h"
#include "fips202.h"

/*************************************************
* Name:        OSKR_shake128_absorb
*
* Description: Absorb step of the SHAKE128 specialized for the OSKR context.
*
* Arguments:   - keccak_state *state: pointer to (uninitialized) output Keccak state
*              - const uint8_t *seed: pointer to OSKR_SYMBYTES input to be absorbed into state
*              - uint8_t i: additional byte of input
*              - uint8_t j: additional byte of input
**************************************************/
void OSKR_shake128_absorb(keccak_state *state,
                           const uint8_t seed[OSKR_SYMBYTES],
                           uint8_t x,
                           uint8_t y)
{
  uint8_t extseed[OSKR_SYMBYTES+2];

  memcpy(extseed, seed, OSKR_SYMBYTES);
  extseed[OSKR_SYMBYTES+0] = x;
  extseed[OSKR_SYMBYTES+1] = y;

  shake128_absorb_once(state, extseed, sizeof(extseed));
}

/*************************************************
* Name:        OSKR_shake256_prf
*
* Description: Usage of SHAKE256 as a PRF, concatenates secret and public input
*              and then generates outlen bytes of SHAKE256 output
*
* Arguments:   - uint8_t *out: pointer to output
*              - size_t outlen: number of requested output bytes
*              - const uint8_t *key: pointer to the key (of length OSKR_SYMBYTES)
*              - uint8_t nonce: single-byte nonce (public PRF input)
**************************************************/
void OSKR_shake256_prf(uint8_t *out, size_t outlen, const uint8_t key[OSKR_SYMBYTES], uint8_t nonce)
{
  uint8_t extkey[OSKR_SYMBYTES+1];

  memcpy(extkey, key, OSKR_SYMBYTES);
  extkey[OSKR_SYMBYTES] = nonce;

  shake256(out, outlen, extkey, sizeof(extkey));
}
