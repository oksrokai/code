#include "api.h"
#include "randombytes.h"
#include "hal.h"
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <poly.h>
#include <polyvec.h>
#include <ntt.h>
//#include <matacc.h>
#include <indcpa.h>
#define NTESTS 10
// use different names so we can have empty namespaces
#define MUPQ_CRYPTO_BYTES           CRYPTO_BYTES
#define MUPQ_CRYPTO_PUBLICKEYBYTES  CRYPTO_PUBLICKEYBYTES
#define MUPQ_CRYPTO_SECRETKEYBYTES  CRYPTO_SECRETKEYBYTES
#define MUPQ_CRYPTO_CIPHERTEXTBYTES CRYPTO_CIPHERTEXTBYTES
#define MUPQ_CRYPTO_ALGNAME CRYPTO_ALGNAME

#define MUPQ_crypto_kem_keypair crypto_kem_keypair
#define MUPQ_crypto_kem_enc crypto_kem_enc
#define MUPQ_crypto_kem_dec crypto_kem_dec

const uint8_t canary[8] = {
  0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02
};

void print8(uint8_t *x, unsigned long long xlen)
{
  char outs[5*xlen+1];
  unsigned long long i;
  for(i=0;i<xlen;i++)
    sprintf(outs+5*i, "%04d,", x[i]);
  outs[5*xlen] = 0;
  hal_send_str(outs);
}

void print16(int16_t *x, unsigned long long xlen)
{
  char outs[5*xlen+1];
  unsigned long long i;
  for(i=0;i<xlen;i++)
    sprintf(outs+5*i, "%04d,", x[i]);
  outs[5*xlen] = 0;
  hal_send_str(outs);
}

void print32(int32_t *x, unsigned long long xlen)
{
  char outs[9*xlen+1];
  unsigned long long i;
  for(i=0;i<xlen;i++)
    sprintf(outs+9*i, "%08lX,", x[i]);
  outs[9*xlen] = 0;
  hal_send_str(outs);
}
/* allocate a bit more for all keys and messages and
 * make sure it is not touched by the implementations.
 */
static void write_canary(uint8_t *d) {
  for (size_t i = 0; i < 8; i++) {
    d[i] = canary[i];
  }
}

static int check_canary(const uint8_t *d) {
  for (size_t i = 0; i < 8; i++) {
    if (d[i] != canary[i]) {
      return -1;
    }
  }
  return 0;
}

static int test_keys(void)
{
  unsigned char key_a[MUPQ_CRYPTO_BYTES+16], key_b[MUPQ_CRYPTO_BYTES+16];
  unsigned char pk[MUPQ_CRYPTO_PUBLICKEYBYTES+16];
  unsigned char sendb[MUPQ_CRYPTO_CIPHERTEXTBYTES+16];
  unsigned char sk_a[MUPQ_CRYPTO_SECRETKEYBYTES+16];

  write_canary(key_a); write_canary(key_a+sizeof(key_a)-8);
  write_canary(key_b); write_canary(key_b+sizeof(key_b)-8);
  write_canary(pk); write_canary(pk+sizeof(pk)-8);
  write_canary(sendb); write_canary(sendb+sizeof(sendb)-8);
  write_canary(sk_a); write_canary(sk_a+sizeof(sk_a)-8);


  int i;

  for(i=0; i<NTESTS; i++)
  {
    //Alice generates a public key
    MUPQ_crypto_kem_keypair(pk+8, sk_a+8);
    // hal_send_str("DONE key pair generation!");
    //printf("00000000000000000000");

    //Bob derives a secret key and creates a response
    MUPQ_crypto_kem_enc(sendb+8, key_b+8, pk+8);
    // hal_send_str("DONE encapsulation!");

    //Alice uses Bobs response to get her secret key
    MUPQ_crypto_kem_dec(key_a+8, sendb+8, sk_a+8);
    // hal_send_str("DONE decapsulation!");

    if(memcmp(key_a+8, key_b+8, MUPQ_CRYPTO_BYTES))
    {
      hal_send_str("ERROR KEYS\n");
    }
    else if(check_canary(key_a) || check_canary(key_a+sizeof(key_a)-8) ||
            check_canary(key_b) || check_canary(key_b+sizeof(key_b)-8) ||
            check_canary(pk) || check_canary(pk+sizeof(pk)-8) ||
            check_canary(sendb) || check_canary(sendb+sizeof(sendb)-8) ||
            check_canary(sk_a) || check_canary(sk_a+sizeof(sk_a)-8))
    {
      hal_send_str("ERROR canary overwritten\n");
    }
    else
    {
      hal_send_str("OK KEYS\n");
    }
  }

  return 0;
}


static int test_invalid_sk_a(void)
{
  unsigned char sk_a[MUPQ_CRYPTO_SECRETKEYBYTES];
  unsigned char key_a[MUPQ_CRYPTO_BYTES], key_b[MUPQ_CRYPTO_BYTES];
  unsigned char pk[MUPQ_CRYPTO_PUBLICKEYBYTES];
  unsigned char sendb[MUPQ_CRYPTO_CIPHERTEXTBYTES];
  int i;

  for(i=0; i<NTESTS; i++)
  {
    //Alice generates a public key
    MUPQ_crypto_kem_keypair(pk, sk_a);

    //Bob derives a secret key and creates a response
    MUPQ_crypto_kem_enc(sendb, key_b, pk);

    //Replace secret key with random values
    randombytes(sk_a, MUPQ_CRYPTO_SECRETKEYBYTES);

    //Alice uses Bobs response to get her secre key
    MUPQ_crypto_kem_dec(key_a, sendb, sk_a);

    if(!memcmp(key_a, key_b, MUPQ_CRYPTO_BYTES))
    {
      hal_send_str("ERROR invalid sk_a\n");
    }
    else
    {
      hal_send_str("OK invalid sk_a\n");
    }
  }

  return 0;
}


static int test_invalid_ciphertext(void)
{
  unsigned char sk_a[MUPQ_CRYPTO_SECRETKEYBYTES];
  unsigned char key_a[MUPQ_CRYPTO_BYTES], key_b[MUPQ_CRYPTO_BYTES];
  unsigned char pk[MUPQ_CRYPTO_PUBLICKEYBYTES];
  unsigned char sendb[MUPQ_CRYPTO_CIPHERTEXTBYTES];
  int i;
  size_t pos;

  for(i=0; i<NTESTS; i++)
  {
    randombytes((unsigned char *)&pos, sizeof(size_t));

    //Alice generates a public key
    MUPQ_crypto_kem_keypair(pk, sk_a);

    //Bob derives a secret key and creates a response
    MUPQ_crypto_kem_enc(sendb, key_b, pk);

    // Change ciphertext to random value
    randombytes(sendb, sizeof(sendb));

    //Alice uses Bobs response to get her secret key
    MUPQ_crypto_kem_dec(key_a, sendb, sk_a);

    if(!memcmp(key_a, key_b, MUPQ_CRYPTO_BYTES))
    {
      hal_send_str("ERROR invalid ciphertext\n");
    }
    else
    {
      hal_send_str("OK invalid ciphertext\n");
    }
  }

  return 0;
}
//extern void matacc(poly* r, polyvec *b, unsigned char i, const unsigned char *seed, int transposed);

int main(void)
{
  hal_setup(CLOCK_FAST);
  // marker for automated testing
  hal_send_str("==========================");
  test_keys();
  test_invalid_sk_a();
  test_invalid_ciphertext();
  hal_send_str("#");
  while(1);
  return 0;
}
