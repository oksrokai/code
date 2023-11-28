#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include "params.h"
#include "kem.h"
#include "indcpa.h"
#include "verify.h"
#include "symmetric.h"
#include "randombytes.h"

/*************************************************
* Name:        crypto_kem_keypair
*
* Description: Generates public and private key
*              for CCA-secure OSKR key encapsulation mechanism
*
* Arguments:   - uint8_t *pk: pointer to output public key
*                (an already allocated array of OSKR_PUBLICKEYBYTES bytes)
*              - uint8_t *sk: pointer to output private key
*                (an already allocated array of OSKR_SECRETKEYBYTES bytes)
*(pk,sk)
* Returns 0 (success)
**************************************************/
int crypto_kem_keypair(uint8_t *pk,
                       uint8_t *sk)
{
  indcpa_keypair(pk, sk);
  memcpy(sk+OSKR_INDCPA_SECRETKEYBYTES, pk, OSKR_INDCPA_PUBLICKEYBYTES);
#ifdef CCS21
  randombytes(sk+OSKR_INDCPA_SECRETKEYBYTES+OSKR_INDCPA_PUBLICKEYBYTES+OSKR_SYMBYTES, OSKR_SYMBYTES);
#else
  hash_h(sk+OSKR_SECRETKEYBYTES-2*OSKR_SYMBYTES, pk, OSKR_PUBLICKEYBYTES);
  randombytes(sk+OSKR_SECRETKEYBYTES-OSKR_SYMBYTES, OSKR_SYMBYTES);
#endif
  return 0;
}

/*************************************************
* Name:        crypto_kem_enc
*
* Description: Generates cipher text and shared
*              secret for given public key
*
* Arguments:   - uint8_t *ct: pointer to output cipher text
*                (an already allocated array of OSKR_CIPHERTEXTBYTES bytes)
*              - uint8_t *ss: pointer to output shared secret
*                (an already allocated array of OSKR_SSBYTES bytes)
*              - const uint8_t *pk: pointer to input public key
*                (an already allocated array of OSKR_PUBLICKEYBYTES bytes)
*
* Returns 0 (success)
**************************************************/
int crypto_kem_enc(uint8_t *ct,
                   uint8_t *ss,
                   const uint8_t *pk)
{
  uint8_t kr[2*OSKR_SYMBYTES];
#ifdef CCS21
  uint8_t buf[PREFIXLEN+OSKR_SYMBYTES];
  randombytes(buf+PREFIXLEN, OSKR_SYMBYTES);
  memcpy(buf, pk, PREFIXLEN);
  kdf(kr, buf, PREFIXLEN+OSKR_SYMBYTES);
  indcpa_enc(ct, buf+PREFIXLEN, pk, kr+OSKR_SYMBYTES);
  memcpy(ss, kr, OSKR_SYMBYTES);
#else
  uint8_t buf[2*OSKR_SYMBYTES];
  /* Will contain key, coins */
  randombytes(buf, OSKR_SYMBYTES);
  /* Don't release system RNG output */
  hash_h(buf, buf, OSKR_SYMBYTES);
  /* Multitarget countermeasure for coins + contributory KEM */
  hash_h(buf+OSKR_SYMBYTES, pk, OSKR_PUBLICKEYBYTES);
  hash_g(kr, buf, 2*OSKR_SYMBYTES);
  /* coins are in kr+OSKR_SYMBYTES */
  indcpa_enc(ct, buf, pk, kr+OSKR_SYMBYTES);
  /* overwrite coins in kr with H(c) */
  sha3_256(kr+OSKR_SYMBYTES, ct, OSKR_CIPHERTEXTBYTES);
  /* hash concatenation of pre-k and H(c) to k */
  kdf(ss, kr, OSKR_SYMBYTES+32);
#endif
  return 0;
}

/*************************************************
* Name:        crypto_kem_dec
*
* Description: Generates shared secret for given
*              cipher text and private key
*
* Arguments:   - uint8_t *ss: pointer to output shared secret
*                (an already allocated array of OSKR_SSBYTES bytes)
*              - const uint8_t *ct: pointer to input cipher text
*                (an already allocated array of OSKR_CIPHERTEXTBYTES bytes)
*              - const uint8_t *sk: pointer to input private key
*                (an already allocated array of OSKR_SECRETKEYBYTES bytes)
*
* Returns 0.
*
* On failure, ss will contain a pseudo-random value.
**************************************************/
int crypto_kem_dec(uint8_t *ss,
                   const uint8_t *ct,
                   const uint8_t *sk)
{
  int fail;
  uint8_t kr[2*OSKR_SYMBYTES];
  ALIGNED_UINT8(OSKR_CIPHERTEXTBYTES) cmp;
  const uint8_t *pk = sk+OSKR_INDCPA_SECRETKEYBYTES;
#ifdef CCS21
  uint8_t buf[PREFIXLEN+OSKR_SYMBYTES];
  uint8_t buf2[PREFIXLEN+OSKR_SYMBYTES+OSKR_CIPHERTEXTBYTES];
  uint8_t buf3[2*OSKR_SYMBYTES];
  indcpa_dec(buf+PREFIXLEN, ct, sk);
  memcpy(buf, pk, PREFIXLEN);
  kdf(kr, buf, PREFIXLEN+OSKR_SYMBYTES);
  indcpa_enc(cmp.coeffs, buf+PREFIXLEN, pk, kr+OSKR_SYMBYTES);
#else
  uint8_t buf[2*OSKR_SYMBYTES];
  /* Will contain key, coins */
  indcpa_dec(buf, ct, sk);
  /* Multitarget countermeasure for coins + contributory KEM */
  memcpy(buf+OSKR_SYMBYTES, sk+OSKR_SECRETKEYBYTES-2*OSKR_SYMBYTES, OSKR_SYMBYTES);
  hash_g(kr, buf, 2*OSKR_SYMBYTES);
  /* coins are in kr+OSKR_SYMBYTES */
  indcpa_enc(cmp.coeffs, buf, pk, kr+OSKR_SYMBYTES);
#endif
  fail = verify(ct, cmp.coeffs, OSKR_CIPHERTEXTBYTES);
#ifdef CCS21
  memcpy(buf2, pk, PREFIXLEN);
  memcpy(buf2+PREFIXLEN, sk+OSKR_INDCPA_SECRETKEYBYTES+OSKR_INDCPA_PUBLICKEYBYTES+OSKR_SYMBYTES, OSKR_SYMBYTES);
  memcpy(buf2+PREFIXLEN+OSKR_SYMBYTES,ct,OSKR_CIPHERTEXTBYTES);
  kdf1(buf3, buf2, PREFIXLEN+OSKR_SYMBYTES+OSKR_CIPHERTEXTBYTES);
  for (int i = 0; i < OSKR_SYMBYTES; ++i)
    ss[i] = kr[i] ^ ((-fail) & (kr[i] ^ buf3[i]));
#else
  /* overwrite coins in kr with H(c) */
  sha3_256(kr+OSKR_SYMBYTES, ct, OSKR_CIPHERTEXTBYTES);
  /* Overwrite pre-k with z on re-encryption failure */
  cmov(kr, sk+OSKR_SECRETKEYBYTES-OSKR_SYMBYTES, OSKR_SYMBYTES, fail);
  /* hash concatenation of pre-k and H(c) to k */
  kdf(ss, kr, OSKR_SYMBYTES+32);
#endif
  return 0;
}