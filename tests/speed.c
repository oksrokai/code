#include "api.h"
#include "hal.h"
#include "sendfn.h"

#include <stdint.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <poly.h>
#include <polyvec.h>
#include <ntt.h>
#include <indcpa.h>

// use different names so we can have empty namespaces
#define MUPQ_CRYPTO_BYTES           CRYPTO_BYTES
#define MUPQ_CRYPTO_PUBLICKEYBYTES  CRYPTO_PUBLICKEYBYTES
#define MUPQ_CRYPTO_SECRETKEYBYTES  CRYPTO_SECRETKEYBYTES
#define MUPQ_CRYPTO_CIPHERTEXTBYTES CRYPTO_CIPHERTEXTBYTES
#define MUPQ_CRYPTO_ALGNAME CRYPTO_ALGNAME

#define MUPQ_crypto_kem_keypair crypto_kem_keypair
#define MUPQ_crypto_kem_enc crypto_kem_enc
#define MUPQ_crypto_kem_dec crypto_kem_dec

#define printcycles(S, U) send_unsignedll((S), (U))

int main(void)
{
  unsigned char key_a[MUPQ_CRYPTO_BYTES], key_b[MUPQ_CRYPTO_BYTES];
  unsigned char sk[MUPQ_CRYPTO_SECRETKEYBYTES];
  unsigned char pk[MUPQ_CRYPTO_PUBLICKEYBYTES];
  unsigned char ct[MUPQ_CRYPTO_CIPHERTEXTBYTES];
  unsigned long long t0, t1;

  hal_setup(CLOCK_BENCHMARK);

  hal_send_str("==========================");

  // Key-pair generation
  t0 = hal_get_time();
  MUPQ_crypto_kem_keypair(pk, sk);
  t1 = hal_get_time();
  printcycles("keypair cycles: ", t1-t0);

  // Encapsulation
  t0 = hal_get_time();
  MUPQ_crypto_kem_enc(ct, key_a, pk);
  t1 = hal_get_time();
  printcycles("encaps cycles: ", t1-t0);

  // Decapsulation
  t0 = hal_get_time();
  MUPQ_crypto_kem_dec(key_b, ct, sk);
  t1 = hal_get_time();
  printcycles("decaps cycles: ", t1-t0);

  if (memcmp(key_a, key_b, MUPQ_CRYPTO_BYTES)) {
    hal_send_str("ERROR KEYS");
  }
  else {
    hal_send_str("OK KEYS");
  }

  hal_send_str("#");
  while(1);
  return 0;
}