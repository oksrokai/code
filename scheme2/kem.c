#include "api.h"
#include "randombytes.h"
#include "fips202.h"
#include "params.h"
#include "verify.h"
#include "indcpa.h"
#include "symmetric.h"

/* Build a CCA-secure KEM from an IND-CPA-secure encryption scheme */

int crypto_kem_keypair(unsigned char *pk, unsigned char *sk)
{
  indcpa_keypair(pk, sk);
  memcpy(sk+OKAI_INDCPA_SECRETKEYBYTES, pk, OKAI_INDCPA_PUBLICKEYBYTES);
#ifdef CCS21
  randombytes(sk+OKAI_INDCPA_SECRETKEYBYTES+OKAI_INDCPA_PUBLICKEYBYTES+OKAI_SYMBYTES, OKAI_SYMBYTES);
#else
  hash_h(sk+OKAI_SECRETKEYBYTES-2*OKAI_SYMBYTES,pk,OKAI_PUBLICKEYBYTES);
  randombytes(sk+OKAI_SECRETKEYBYTES-OKAI_SYMBYTES,OKAI_SYMBYTES);
#endif
  return 0;
}

int crypto_kem_enc(unsigned char *ct, unsigned char *ss, const unsigned char *pk)
{
  __attribute__((aligned(32)))
  unsigned char kr[2*OKAI_SYMBYTES];
#ifdef CCS21
  __attribute__((aligned(32)))
  uint8_t buf[PREFIXLEN+OKAI_SYMBYTES];
  randombytes(buf+PREFIXLEN, OKAI_SYMBYTES);
  memcpy(buf, pk, PREFIXLEN);
  kdf(kr, buf, PREFIXLEN+OKAI_SYMBYTES);
  indcpa_enc(ct, buf+PREFIXLEN, pk, kr+OKAI_SYMBYTES);
  memcpy(ss, kr, OKAI_SYMBYTES);
#else
  __attribute__((aligned(32)))
  unsigned char buf[2*OKAI_SYMBYTES];
  /* Will contain key, coins */
  randombytes(buf, OKAI_SYMBYTES);
  hash_h(buf,buf,OKAI_SYMBYTES);  /* Don't release system RNG output */

  hash_h(buf+OKAI_SYMBYTES, pk, OKAI_PUBLICKEYBYTES);  /* Multitarget countermeasure for coins + contributory KEM */
  hash_g(kr, buf, 2*OKAI_SYMBYTES);

  indcpa_enc(ct, buf, pk, kr+OKAI_SYMBYTES);  /* coins are in kr+OKAI_SYMBYTES */

  hash_h(kr+OKAI_SYMBYTES, ct, OKAI_CIPHERTEXTBYTES);  /* overwrite coins in kr with H(c) */
  kdf(ss, kr, 2*OKAI_SYMBYTES);  /* hash concatenation of pre-k and H(c) to k */
#endif
  return 0;
}

int crypto_kem_dec(unsigned char *ss, const unsigned char *ct, const unsigned char *sk)
{
  int fail;
  unsigned char cmp[OKAI_CIPHERTEXTBYTES];
  __attribute__((aligned(32)))
  unsigned char kr[2*OKAI_SYMBYTES];
  const unsigned char *pk = sk+OKAI_INDCPA_SECRETKEYBYTES;
#ifdef CCS21
  uint8_t buf[PREFIXLEN+OKAI_SYMBYTES];
  uint8_t buf2[PREFIXLEN+OKAI_SYMBYTES+OKAI_CIPHERTEXTBYTES];
  uint8_t buf3[2*OKAI_SYMBYTES];
  indcpa_dec(buf+PREFIXLEN, ct, sk);
  memcpy(buf, pk, PREFIXLEN);
  kdf(kr, buf, PREFIXLEN+OKAI_SYMBYTES);
  indcpa_enc(cmp, buf+PREFIXLEN, pk, kr+OKAI_SYMBYTES);
#else
  __attribute__((aligned(32)))
  uint8_t buf[2*OKAI_SYMBYTES];
  indcpa_dec(buf, ct, sk);
  /* Multitarget countermeasure for coins + contributory KEM */
  memcpy(buf+OKAI_SYMBYTES,sk+OKAI_SECRETKEYBYTES-2*OKAI_SYMBYTES,OKAI_SYMBYTES);
  hash_g(kr, buf, 2*OKAI_SYMBYTES);
  indcpa_enc(cmp, buf, pk, kr+OKAI_SYMBYTES);
#endif

  fail = verify(ct, cmp, OKAI_CIPHERTEXTBYTES);

#ifdef CCS21
  memcpy(buf2, pk, PREFIXLEN);
  memcpy(buf2+PREFIXLEN, sk+OKAI_INDCPA_SECRETKEYBYTES+OKAI_INDCPA_PUBLICKEYBYTES+OKAI_SYMBYTES, OKAI_SYMBYTES);
  memcpy(buf2+PREFIXLEN+OKAI_SYMBYTES,ct,OKAI_CIPHERTEXTBYTES);
  kdf1(buf3, buf2, PREFIXLEN+OKAI_SYMBYTES+OKAI_CIPHERTEXTBYTES);
  for (int i = 0; i < OKAI_SYMBYTES; ++i)
    ss[i] = kr[i] ^ ((-fail) & (kr[i] ^ buf3[i]));
  return fail;
#else
  hash_h(kr+OKAI_SYMBYTES, ct, OKAI_CIPHERTEXTBYTES);                     /* overwrite coins in kr with H(c)  */
  cmov(kr, sk+OKAI_SECRETKEYBYTES-OKAI_SYMBYTES, OKAI_SYMBYTES, fail);     /* Overwrite pre-k with z on re-encryption failure */
  kdf(ss, kr, 2*OKAI_SYMBYTES);                                         /* hash concatenation of pre-k and H(c) to k */
  return -fail;
#endif
}