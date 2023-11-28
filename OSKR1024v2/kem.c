#include "api.h"
#include "indcpa.h"
#include "params.h"
#include "randombytes.h"
#include "symmetric.h"
#include "verify.h"

#include <stdlib.h>
#include <string.h>
#include <stdlib.h>

/*************************************************
* Name:        crypto_kem_keypair
*
* Description: Generates public and private key
*              for CCA-secure OSKR key encapsulation mechanism
*
* Arguments:   - unsigned char *pk: pointer to output public key (an already allocated array of CRYPTO_PUBLICKEYBYTES bytes)
*              - unsigned char *sk: pointer to output private key (an already allocated array of CRYPTO_SECRETKEYBYTES bytes)
*
* Returns 0 (success)
**************************************************/
int crypto_kem_keypair(unsigned char *pk, unsigned char *sk) {
    indcpa_keypair(pk, sk);//CPAPKE's pk sk
    memcpy(sk+OSKR_INDCPA_SECRETKEYBYTES, pk, OSKR_INDCPA_PUBLICKEYBYTES);//pos change (sk||pk)
    randombytes(sk+OSKR_INDCPA_SECRETKEYBYTES+OSKR_INDCPA_PUBLICKEYBYTES+OSKR_SYMBYTES, OSKR_SYMBYTES);//pos change (sk||pk|| ||z):=sk
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
int crypto_kem_enc(unsigned char *ct, unsigned char *K, const unsigned char *pk) {
    uint8_t buf[PREFIXLEN+OSKR_SYMBYTES];
    /* Will contain key, coins */
    uint8_t kr[2*OSKR_SYMBYTES];

    randombytes(buf+PREFIXLEN, OSKR_SYMBYTES);//m

    memcpy(buf, pk, PREFIXLEN);

    /* Don't release system RNG output */
    kdf(kr, buf, PREFIXLEN+OSKR_SYMBYTES);//G(ID(pk),m)

    /* coins are in kr+OSKR_SYMBYTES */
    indcpa_enc(ct, buf+PREFIXLEN, pk, kr+OSKR_SYMBYTES);//CPAPKE.Enc(pk,m,r)->ct
    
    memcpy(K, kr, OSKR_SYMBYTES);//KDF(K`||H(c))->K(shared key)
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
* Returns 0.
*
* On failure, ss will contain a pseudo-random value.
**************************************************/
int crypto_kem_dec(unsigned char *K, const unsigned char *ct, const unsigned char *sk) {                                         /* Will contain key, coins */
    unsigned char cmp[OSKR_CIPHERTEXTBYTES];
    int fail;
    uint8_t buf[PREFIXLEN+OSKR_SYMBYTES];
    uint8_t buf2[PREFIXLEN+OSKR_SYMBYTES+OSKR_CIPHERTEXTBYTES];
    /* Will contain key, coins */
    uint8_t kr[2*OSKR_SYMBYTES];
    uint8_t buf3[2*OSKR_SYMBYTES];
    const uint8_t *pk = sk+OSKR_INDCPA_SECRETKEYBYTES;

    indcpa_dec(buf+PREFIXLEN, ct, sk);//CPAPKE.Dec(s,(u,v))->m`

    /* Multitarget countermeasure for coins + contributory KEM */
    memcpy(buf, pk, PREFIXLEN);//64B hash(pk)

    kdf(kr, buf, PREFIXLEN+OSKR_SYMBYTES);
    indcpa_enc(cmp, buf+PREFIXLEN, pk, kr + OSKR_SYMBYTES);                          /* coins are in kr+OSKR_SYMBYTES */
    fail = verify(ct, cmp, OSKR_CIPHERTEXTBYTES);
    
    memcpy(buf2, pk, PREFIXLEN);
    memcpy(buf2+PREFIXLEN, sk+OSKR_INDCPA_SECRETKEYBYTES+OSKR_INDCPA_PUBLICKEYBYTES+OSKR_SYMBYTES, OSKR_SYMBYTES);
    memcpy(buf2+PREFIXLEN+OSKR_SYMBYTES,ct,OSKR_CIPHERTEXTBYTES);

    /* overwrite coins in kr with H(c) */
    kdf1(buf3, buf2, PREFIXLEN+OSKR_SYMBYTES+OSKR_CIPHERTEXTBYTES);//(K`||H(c)) H(c):32B

    /* Overwrite pre-k with z on re-encryption failure */
    for (int i = 0; i < OSKR_SYMBYTES; ++i)
        K[i] = kr[i] ^ ((-fail) & (kr[i] ^ buf3[i]));
    return fail;
}
