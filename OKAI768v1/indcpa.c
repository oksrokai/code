#include <stdint.h>
#include "indcpa.h"
#include "poly.h"
#include "polyvec.h"
// #include "rng.h"
#include "ntt.h"
#include "symmetric.h"

/*************************************************
* Name:        pack_pk
*
* Description: Serialize the public key as concatenation of the
*              compressed and serialized vector of polynomials pk
*              and the public seed used to generate the matrix A.
*
* Arguments:   unsigned char *r:          pointer to the output serialized public key
*              const poly *pk:            pointer to the input public-key polynomial
*              const unsigned char *seed: pointer to the input public seed
**************************************************/
static void pack_pk(unsigned char *r, polyvec *pk, const unsigned char *seed)
{
  int i;
  polyvec_compress(r, pk);
  for(i=0;i<OKAI_SYMBYTES;i++)
    r[i+OKAI_POLYVECCOMPRESSEDBYTES] = seed[i];
}

/*************************************************
* Name:        unpack_pk
*
* Description: De-serialize and decompress public key from a byte array;
*              approximate inverse of pack_pk
*
* Arguments:   - polyvec *pk:                   pointer to output public-key vector of polynomials
*              - unsigned char *seed:           pointer to output seed to generate matrix A
*              - const unsigned char *packedpk: pointer to input serialized public key
**************************************************/
static void unpack_pk(polyvec *pk, unsigned char *seed, const unsigned char *packedpk)
{
  int i;
  polyvec_decompress(pk, packedpk);

  for(i=0;i<OKAI_SYMBYTES;i++)
    seed[i] = packedpk[i+OKAI_POLYVECCOMPRESSEDBYTES];
}

/*************************************************
* Name:        pack_ciphertext
*
* Description: Serialize the ciphertext as concatenation of the
*              compressed and serialized vector of polynomials b
*              and the compressed and serialized polynomial v
*
* Arguments:   unsigned char *r:          pointer to the output serialized ciphertext
*              const poly *pk:            pointer to the input vector of polynomials b
*              const unsigned char *seed: pointer to the input polynomial v
**************************************************/
static void pack_ciphertext(unsigned char *r, polyvec *b, poly *v)
{
  polyvec_compress(r, b);
  // poly_compress(r+OKAI_POLYVECCOMPRESSEDBYTES, v);
}

/*************************************************
* Name:        unpack_ciphertext
*
* Description: De-serialize and decompress ciphertext from a byte array;
*              approximate inverse of pack_ciphertext
*
* Arguments:   - polyvec *b:             pointer to the output vector of polynomials b
*              - poly *v:                pointer to the output polynomial v
*              - const unsigned char *c: pointer to the input serialized ciphertext
**************************************************/
static void unpack_ciphertext(polyvec *b, poly *v, const unsigned char *c)
{
  polyvec_decompress(b, c);
  poly_decompress(v, c+OKAI_POLYVECCOMPRESSEDBYTES);
}

/*************************************************
* Name:        pack_sk
*
* Description: Serialize the secret key
*
* Arguments:   - unsigned char *r:  pointer to output serialized secret key
*              - const polyvec *sk: pointer to input vector of polynomials (secret key)
**************************************************/
static void pack_sk(unsigned char *r, const polyvec *sk)
{
  polyvec_tobytes(r, sk);
}

/*************************************************
* Name:        unpack_sk
*
* Description: De-serialize the secret key;
*              inverse of pack_sk
*
* Arguments:   - polyvec *sk:                   pointer to output vector of polynomials (secret key)
*              - const unsigned char *packedsk: pointer to input serialized secret key
**************************************************/
static void unpack_sk(polyvec *sk, const unsigned char *packedsk)
{
  polyvec_frombytes(sk, packedsk);
}

/*************************************************
* Name:        rej_uniform
*
* Description: Run rejection sampling on uniform random bytes to generate
*              uniform random integers mod q
*
* Arguments:   - int16_t *r:               pointer to output buffer
*              - unsigned int len:         requested number of 16-bit integers (uniform mod q)
*              - const unsigned char *buf: pointer to input buffer (assumed to be uniform random bytes)
*              - unsigned int buflen:      length of input buffer in bytes
*
* Returns number of sampled 16-bit integers (at most len)
**************************************************/
static unsigned int rej_uniform(int16_t *r, unsigned int len, const unsigned char *buf, unsigned int buflen)
{
  unsigned int ctr, pos;
  uint16_t val;
  const uint32_t v = (1U << 26)/OKAI_Q + 1;

  ctr = pos = 0;
  while(ctr < len && pos + 2 <= buflen)
  {
    val = buf[pos] | ((uint16_t)buf[pos+1] << 8);
    pos += 2;

    if(val < 8*OKAI_Q)
    {
      val -= ((int32_t)(v*val) >> 26) * OKAI_Q; // Barrett reduction
      r[ctr++] = (int16_t)val;
    }
  } 

  return ctr;
}

#define gen_a(A,B)  gen_matrix(A,B,0)
#define gen_at(A,B) gen_matrix(A,B,1)

/*************************************************
* Name:        gen_matrix
*
* Description: Deterministically generate matrix A (or the transpose of A)
*              from a seed. Entries of the matrix are polynomials that look
*              uniformly random. Performs rejection sampling on output of
*              a XOF
*
* Arguments:   - polyvec *a:                pointer to ouptput matrix A
*              - const unsigned char *seed: pointer to input seed
*              - int transposed:            boolean deciding whether A or A^T is generated
**************************************************/
void gen_matrix(polyvec *a, const unsigned char *seed, int transposed) // Not static for benchmarking
{
  unsigned int ctr, i, j;
  const unsigned int maxnblocks=(530+XOF_BLOCKBYTES)/XOF_BLOCKBYTES; /* 530 is expected number of required bytes */
  unsigned char buf[XOF_BLOCKBYTES*maxnblocks+1];
  xof_state state;

  for(i=0;i<OKAI_K;i++)
  {
    for(j=0;j<OKAI_K;j++)
    {
      if(transposed) {
        xof_absorb(&state, seed, i, j);
      }
      else {
        xof_absorb(&state, seed, j, i);
      }

      xof_squeezeblocks(buf, maxnblocks, &state);
      ctr = rej_uniform(a[i].vec[j].coeffs, OKAI_N, buf, maxnblocks*XOF_BLOCKBYTES);

      while(ctr < OKAI_N)
      {
        xof_squeezeblocks(buf, 1, &state);
        ctr += rej_uniform(a[i].vec[j].coeffs + ctr, OKAI_N - ctr, buf, XOF_BLOCKBYTES);
      }
    }
  }
}

/*************************************************
* Name:        indcpa_keypair
*
* Description: Generates public and private key for the CPA-secure
*              public-key encryption scheme underlying OKAI
*
* Arguments:   - unsigned char *pk: pointer to output public key
*              - unsigned char *sk: pointer to output private key
**************************************************/
void indcpa_keypair(unsigned char *pk,
                   unsigned char *sk)
{
  polyvec a[OKAI_K], e, pkpv, skpv;
  unsigned char buf[OKAI_SYMBYTES+OKAI_SYMBYTES];
  unsigned char *publicseed = buf;
  unsigned char *noiseseed = buf+OKAI_SYMBYTES;
  int i;
  unsigned char nonce=0;

  randombytes(buf, OKAI_SYMBYTES);
  hash_g(buf, buf, OKAI_SYMBYTES);

  gen_a(a, publicseed);

  for(i=0;i<OKAI_K;i++)
    poly_getsecret(skpv.vec+i,noiseseed,nonce++);

  polyvec_ntt(&skpv);

  for(i=0;i<OKAI_K;i++)
    poly_getnoise(e.vec+i,noiseseed,nonce++);

  // matrix-vector multiplication
  for(i=0;i<OKAI_K;i++)
    polyvec_pointwise_acc(&pkpv.vec[i],&skpv,a+i);

  polyvec_invntt(&pkpv);
  polyvec_add(&pkpv,&pkpv,&e);
  polyvec_reduce(&pkpv);

  pack_sk(sk, &skpv);
  pack_pk(pk, &pkpv, publicseed);
}


/*************************************************
* Name:        indcpa_enc
*
* Description: Encryption function of the CPA-secure
*              public-key encryption scheme underlying OKAI.
*
* Arguments:   - unsigned char *c:          pointer to output ciphertext
*              - const unsigned char *m:    pointer to input message (of length OKAI_SYMBYTES bytes)
*              - const unsigned char *pk:   pointer to input public key
*              - const unsigned char *coin: pointer to input random coins used as seed
*                                           to deterministically generate all randomness
**************************************************/
void indcpa_enc(unsigned char *c,
               const unsigned char *m,
               const unsigned char *pk,
               const unsigned char *coins)
{
  polyvec sp, pkpv, ep, at[OKAI_K], bp;
  poly v, k, epp;
  unsigned char seed[OKAI_SYMBYTES];
  int i;
  unsigned char nonce=0;

  unpack_pk(&pkpv, seed, pk);
  // poly_frommsg(&k, m);

  polyvec_ntt(&pkpv);
  gen_at(at, seed);

  for(i=0;i<OKAI_K;i++)
    poly_getsecret(sp.vec+i,coins,nonce++);
  for(i=0;i<OKAI_K;i++)
    poly_getnoise(ep.vec+i,coins,nonce++);
  poly_getnoise(&epp,coins,nonce++);

  polyvec_ntt(&sp);

  // matrix-vector multiplication
  for(i=0;i<OKAI_K;i++)
    polyvec_pointwise_acc(&bp.vec[i],&sp,at+i);
  polyvec_invntt(&bp);

  polyvec_pointwise_acc(&v, &pkpv, &sp);
  poly_invntt(&v);

  polyvec_add(&bp, &bp, &ep);
  polyvec_reduce(&bp);
  
  poly_add(&v, &v, &epp);
  // poly_add(&v, &v, &k);
  // poly_reduce(&v);

  // pack_ciphertext(c, &bp, &v);
  polyvec_compress(c, &bp);
  poly_con_4(c+OKAI_POLYVECCOMPRESSEDBYTES,m,v.coeffs);
}

/*************************************************
* Name:        indcpa_dec
*
* Description: Decryption function of the CPA-secure
*              public-key encryption scheme underlying OKAI.
*
* Arguments:   - unsigned char *m:        pointer to output decrypted message
*              - const unsigned char *c:  pointer to input ciphertext
*              - const unsigned char *sk: pointer to input secret key
**************************************************/
void indcpa_dec(unsigned char *m,
               const unsigned char *c,
               const unsigned char *sk)
{
  polyvec bp, skpv;
  poly v, mp;
  // int tmp;

  // unpack_ciphertext(&bp, &v, c);
  polyvec_decompress(&bp, c);
  unpack_sk(&skpv, sk);

  polyvec_ntt(&bp);
  polyvec_pointwise_acc(&mp,&skpv,&bp);
  poly_invntt(&mp);

  // poly_sub(&mp, &v, &mp);
  poly_reduce(&mp);

  poly_rec_4(m,c+OKAI_POLYVECCOMPRESSEDBYTES,mp.coeffs);
}
