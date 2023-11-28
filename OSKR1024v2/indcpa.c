#include "indcpa.h"
#include "ntt.h"
#include "poly.h"
#include "polyvec.h"
#include "randombytes.h"
#include "symmetric.h"
#include "matacc.h"

#include <string.h>
#include <stdint.h>

void unpack_sk(polyvec *sk, const uint8_t packedsk[OSKR_INDCPA_SECRETKEYBYTES])
{
  polyvec_frombytes(sk, packedsk);
}
/*************************************************
* Name:        indcpa_keypair
*
* Description: Generates public and private key for the CPA-secure
*              public-key encryption scheme underlying OSKR
*
* Arguments:   - unsigned char *pk: pointer to output public key (of length OSKR_INDCPA_PUBLICKEYBYTES bytes)
*              - unsigned char *sk: pointer to output private key (of length OSKR_INDCPA_SECRETKEYBYTES bytes)
**************************************************/
void indcpa_keypair(unsigned char *pk, unsigned char *sk) {
    polyvec skpv;
    poly pkp;
    unsigned char buf[2 * OSKR_SYMBYTES];
    unsigned char *publicseed = buf;
    unsigned char *noiseseed = buf + OSKR_SYMBYTES;
    int i;
    unsigned char nonce = 0;

    randombytes(buf, OSKR_SYMBYTES);
    hash_g(buf, buf, OSKR_SYMBYTES);

    for (i = 0; i < OSKR_K; i++)
        poly_getnoise(skpv.vec + i, noiseseed, nonce++);

    polyvec_ntt(&skpv);
    
    for (i = 0; i < OSKR_K; i++) {
        matacc(&pkp, &skpv, i, publicseed, 0);
        poly_invntt(&pkp);

        poly_addnoise(&pkp, noiseseed, nonce++);
        poly_ntt(&pkp);

        poly_tobytes(pk+i*OSKR_POLYBYTES, &pkp);
    }
    polyvec_tobytes(sk, &skpv);
    memcpy(pk + OSKR_POLYVECBYTES, publicseed, OSKR_SYMBYTES); // Pack the public seed in the public key
}

/*************************************************
* Name:        indcpa_enc
*
* Description: Encryption function of the CPA-secure
*              public-key encryption scheme underlying OSKR.
*
* Arguments:   - unsigned char *c:          pointer to output ciphertext (of length OSKR_INDCPA_BYTES bytes)
*              - const unsigned char *m:    pointer to input message (of length OSKR_INDCPA_MSGBYTES bytes)
*              - const unsigned char *pk:   pointer to input public key (of length OSKR_INDCPA_PUBLICKEYBYTES bytes)
*              - const unsigned char *coin: pointer to input random coins used as seed (of length OSKR_SYMBYTES bytes)
*                                           to deterministically generate all randomness
**************************************************/
void indcpa_enc(unsigned char *c,
               const unsigned char *m,
               const unsigned char *pk,
               const unsigned char *coins) {
    polyvec sp;
    poly bp,t;
    poly *pkp = &bp;
    //poly *k = &bp;
    poly *v = &sp.vec[0];
    const unsigned char *seed = pk+OSKR_POLYVECBYTES;
    int i;
    unsigned char nonce = 0;

    for (i = 0; i < OSKR_K; i++)
        poly_getnoise(sp.vec + i, coins, nonce++);

    polyvec_ntt(&sp);

    for (i = 0; i < OSKR_K; i++) {
        matacc(&bp, &sp, i, seed, 1);
        poly_invntt(&bp);

        poly_addnoise(&bp, coins, nonce++);
        poly_reduce(&bp);

        poly_packcompress(c, &bp, i);
    }

    poly_frombytes(pkp, pk);
    poly_multi_basemul(v, pkp, &sp.vec[0]);
    poly_frombytes(pkp, pk + OSKR_POLYBYTES);
    poly_multi_basemul(&t, pkp, &sp.vec[1]);
    poly_add(v,v,&t);
    poly_invntt(v);
    poly_addnoise(v, coins, nonce++);
    poly_reduce(v);
    poly_con(c+OSKR_POLYVECCOMPRESSEDBYTES, m, v);
}

/*************************************************
* Name:        indcpa_dec
*
* Description: Decryption function of the CPA-secure
*              public-key encryption scheme underlying OSKR.
*
* Arguments:   - unsigned char *m:        pointer to output decrypted message (of length OSKR_INDCPA_MSGBYTES)
*              - const unsigned char *c:  pointer to input ciphertext (of length OSKR_INDCPA_BYTES)
*              - const unsigned char *sk: pointer to input secret key (of length OSKR_INDCPA_SECRETKEYBYTES)
**************************************************/
void indcpa_dec(unsigned char *m,
                const unsigned char *c,
                const unsigned char *sk) {
    poly mp, bp,t;
    polyvec skpv;
    //poly *v = &bp;
    int i;

    poly_unpackdecompress(&mp, c, 0);
    poly_ntt(&mp);
    unpack_sk(&skpv, sk);//sT
    poly_multi_basemul(&mp, &(skpv.vec[0]), &mp);
    poly_unpackdecompress(&bp, c, 1);
    poly_ntt(&bp);
    poly_multi_basemul(&t, &(skpv.vec[1]), &bp);
    poly_add(&mp,&mp,&t);
    poly_invntt(&mp);
    poly_reduce(&mp);
    poly_rec(m, c+OSKR_POLYVECCOMPRESSEDBYTES, &mp);
}
