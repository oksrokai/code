#include "api.h"
#include "randombytes.h"
#include "hal.h"
#include "sendfn.h"

#include <string.h>

#define MAX_SIZE 0x1B000

// https://stackoverflow.com/a/1489985/1711232
#define PASTER(x, y) x####y
#define EVALUATOR(x, y) PASTER(x, y)
#define NAMESPACE(fun) EVALUATOR(MUPQ_NAMESPACE, fun)

// use different names so we can have empty namespaces
#define MUPQ_CRYPTO_BYTES           (CRYPTO_BYTES)
#define MUPQ_CRYPTO_PUBLICKEYBYTES  (CRYPTO_PUBLICKEYBYTES)
#define MUPQ_CRYPTO_SECRETKEYBYTES  (CRYPTO_SECRETKEYBYTES)
#define MUPQ_CRYPTO_CIPHERTEXTBYTES (CRYPTO_CIPHERTEXTBYTES)
#define MUPQ_CRYPTO_ALGNAME (CRYPTO_ALGNAME)

#define MUPQ_crypto_kem_keypair (crypto_kem_keypair)
#define MUPQ_crypto_kem_enc (crypto_kem_enc)
#define MUPQ_crypto_kem_dec (crypto_kem_dec)

#define send_stack_usage(S, U) send_unsigned((S), (U))

unsigned int canary_size = MAX_SIZE;
volatile unsigned char *p;
unsigned int c;
uint8_t canary = 0x42;

unsigned char key_a[MUPQ_CRYPTO_BYTES], key_b[MUPQ_CRYPTO_BYTES];
unsigned char pk[MUPQ_CRYPTO_PUBLICKEYBYTES];
unsigned char sendb[MUPQ_CRYPTO_CIPHERTEXTBYTES];
unsigned char sk_a[MUPQ_CRYPTO_SECRETKEYBYTES];
unsigned int stack_key_gen, stack_encaps, stack_decaps;

#define FILL_STACK()                                                           \
  p = &a;                                                                      \
  while (p > &a - canary_size)                                                 \
    *(p--) = canary;
#define CHECK_STACK()                                                          \
  c = canary_size;                                                             \
  p = &a - canary_size + 1;                                                    \
  while (*p == canary && p < &a) {                                             \
    p++;                                                                       \
    c--;                                                                       \
  }

static int test_keys(void) {
  volatile unsigned char a;
  // Alice generates a public key
  FILL_STACK()
  MUPQ_crypto_kem_keypair(pk, sk_a);
  CHECK_STACK()
  if(c >= canary_size) return -1;
  stack_key_gen = c;

  // Bob derives a secret key and creates a response
  FILL_STACK()
  MUPQ_crypto_kem_enc(sendb, key_b, pk);
  CHECK_STACK()
  if(c >= canary_size) return -1;
  stack_encaps = c;

  // Alice uses Bobs response to get her secret key
  FILL_STACK()
  MUPQ_crypto_kem_dec(key_a, sendb, sk_a);
  CHECK_STACK()
  if(c >= canary_size) return -1;
  stack_decaps = c;

  if (memcmp(key_a, key_b, MUPQ_CRYPTO_BYTES)){
    return -1;
  } else {
    send_stack_usage("keypair stack usage:", stack_key_gen);
    send_stack_usage("encaps stack usage:", stack_encaps);
    send_stack_usage("decaps stack usage:", stack_decaps);
    hal_send_str("OK KEYS\n");
    return 0;
  }
}

int main(void) {
  hal_setup(CLOCK_FAST);

  // marker for automated benchmarks
  hal_send_str("==========================");
  canary_size = 0x1000;
  while(test_keys()){
    canary_size += 0x1000;
    if(canary_size >= MAX_SIZE) {
      hal_send_str("failed to measure stack usage.\n");
      break;
    }
  }
  // marker for automated benchmarks
  hal_send_str("#");

  while (1);

  return 0;
}
