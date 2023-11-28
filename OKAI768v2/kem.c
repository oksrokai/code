#include "api.h"
// #include "rng.h"
// #include "fips202.h"
#include "params.h"
#include "verify.h"
#include "indcpa.h"
#include "symmetric.h"

/*************************************************
* Name:        crypto_kem_keypair
*
* Description: Generates public and private key
*              for CCA-secure OKAI key encapsulation mechanism
*
* Arguments:   - unsigned char *pk: pointer to output public key (an already allocated array of CRYPTO_PUBLICKEYBYTES bytes)
*              - unsigned char *sk: pointer to output private key (an already allocated array of CRYPTO_SECRETKEYBYTES bytes)
*
* Returns 0 (success)
**************************************************/
int crypto_kem_keypair(unsigned char *pk, unsigned char *sk)
{
  indcpa_keypair(pk, sk);//CPAPKE's pk sk
  memcpy(sk+OKAI_INDCPA_SECRETKEYBYTES, pk, OKAI_INDCPA_PUBLICKEYBYTES);//pos change (sk||pk)
  randombytes(sk+OKAI_INDCPA_SECRETKEYBYTES+OKAI_INDCPA_PUBLICKEYBYTES+OKAI_SYMBYTES, OKAI_SYMBYTES);//pos change (sk||pk|| ||z):=sk
  return 0;
}

/*************************************************
* Name:        crypto_kem_enc
*
* Description: Generates cipher text and shared
*              secret for given public key
*
* Arguments:   - unsigned char *ct:       pointer to output cipher text (an already allocated array of CRYPTO_CIPHERTEXTBYTES bytes)
*              - unsigned char *ss:       pointer to output shared secret (an already allocated array of CRYPTO_BYTES bytes)
*              - const unsigned char *pk: pointer to input public key (an already allocated array of CRYPTO_PUBLICKEYBYTES bytes)
*
* Returns 0 (success)
**************************************************/
int crypto_kem_enc(unsigned char *ct, unsigned char *K, const unsigned char *pk)
{
  uint8_t buf[PREFIXLEN+OKAI_SYMBYTES];
  /* Will contain key, coins */
  uint8_t kr[2*OKAI_SYMBYTES];

  randombytes(buf+PREFIXLEN, OKAI_SYMBYTES);//m

  memcpy(buf, pk, PREFIXLEN);

  /* Don't release system RNG output */
  kdf(kr, buf, PREFIXLEN+OKAI_SYMBYTES);//G(ID(pk),m)

  /* coins are in kr+OSKR_SYMBYTES */
  indcpa_enc(ct, buf+PREFIXLEN, pk, kr+OKAI_SYMBYTES);//CPAPKE.Enc(pk,m,r)->ct
  
  memcpy(K, kr, OKAI_SYMBYTES);//KDF(K`||H(c))->K(shared key)
  return 0;
}

/*************************************************
* Name:        crypto_kem_dec
*
* Description: Generates shared secret for given
*              cipher text and private key
*
* Arguments:   - unsigned char *ss:       pointer to output shared secret (an already allocated array of CRYPTO_BYTES bytes)
*              - const unsigned char *ct: pointer to input cipher text (an already allocated array of CRYPTO_CIPHERTEXTBYTES bytes)
*              - const unsigned char *sk: pointer to input private key (an already allocated array of CRYPTO_SECRETKEYBYTES bytes)
*
* Returns 0 for sucess or -1 for failure
*
* On failure, ss will contain a randomized value.
**************************************************/
int crypto_kem_dec(unsigned char *K, const unsigned char *ct, const unsigned char *sk)
{
  int fail;
  uint8_t buf[PREFIXLEN+OKAI_SYMBYTES];
  unsigned char cmp[OKAI_CIPHERTEXTBYTES];
  uint8_t buf2[PREFIXLEN+OKAI_SYMBYTES+OKAI_CIPHERTEXTBYTES];
  /* Will contain key, coins */
  uint8_t kr[2*OKAI_SYMBYTES];
  uint8_t buf3[2*OKAI_SYMBYTES];
  const uint8_t *pk = sk+OKAI_INDCPA_SECRETKEYBYTES;

  indcpa_dec(buf+PREFIXLEN, ct, sk);

  memcpy(buf, pk, PREFIXLEN);//64B hash(pk)

  kdf(kr, buf, PREFIXLEN+OKAI_SYMBYTES);//G(m`||h)->(K`,r)


  indcpa_enc(cmp, buf+PREFIXLEN, pk, kr+OKAI_SYMBYTES);                                /* coins are in kr+OKAI_SYMBYTES */

  fail = verify(ct, cmp, OKAI_CIPHERTEXTBYTES);

  memcpy(buf2, pk, PREFIXLEN);
  memcpy(buf2+PREFIXLEN, sk+OKAI_INDCPA_SECRETKEYBYTES+OKAI_INDCPA_PUBLICKEYBYTES+OKAI_SYMBYTES, OKAI_SYMBYTES);
  memcpy(buf2+PREFIXLEN+OKAI_SYMBYTES,ct,OKAI_CIPHERTEXTBYTES);

  /* overwrite coins in kr with H(c) */
  kdf1(buf3, buf2, PREFIXLEN+OKAI_SYMBYTES+OKAI_CIPHERTEXTBYTES);//(K`||H(c)) H(c):32B

  for (int i = 0; i < OKAI_SYMBYTES; ++i)
    K[i] = kr[i] ^ ((-fail) & (kr[i] ^ buf3[i]));
  return fail;
}
