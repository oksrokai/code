/********************************************************************************************
* Abstract: benchmarking/testing PKE
*********************************************************************************************/

#include <stdio.h>
#include <string.h>
#include "ds_benchmark.h"
#include "params.h"
#include "indcpa.h"
#include "randombytes.h"

#define CRYPTO_ALGNAME "OKAI-MLWE PKE"
#define PKE_TEST_ITERATIONS 1000
#define PKE_BENCH_SECONDS 1
#define FALSE 0
#define TRUE 1

static int pke_test(const char *named_parameters, int iterations)
{
    uint8_t pk[OKAI_INDCPA_PUBLICKEYBYTES];
    uint8_t sk[OKAI_INDCPA_SECRETKEYBYTES];
    uint8_t c[OKAI_INDCPA_BYTES];
    uint8_t m1[OKAI_INDCPA_MSGBYTES], m2[OKAI_INDCPA_MSGBYTES];

    printf("\n");
    printf("Testing correctness of %s, tests for %d iterations\n", named_parameters, iterations);

    for (int i = 0; i < iterations; i++) {
        indcpa_keypair(pk, sk);
        randombytes(m1, OKAI_INDCPA_MSGBYTES);
        pke_enc(c, m1, pk);
        indcpa_dec(m2, c, sk);
        if (memcmp(m1, m2, OKAI_INDCPA_MSGBYTES) != 0) {
            printf("\n ERROR!\n");
	        return FALSE;
        }
    }
    printf("Tests PASSED. All session keys matched.\n");
    printf("\n");

    return TRUE;
}


static void pke_bench(const int seconds)
{
    uint8_t pk[OKAI_INDCPA_PUBLICKEYBYTES];
    uint8_t sk[OKAI_INDCPA_SECRETKEYBYTES];
    uint8_t c[OKAI_INDCPA_BYTES];
    uint8_t m1[OKAI_INDCPA_MSGBYTES], m2[OKAI_INDCPA_MSGBYTES];

    TIME_OPERATION_SECONDS({ indcpa_keypair(pk, sk); }, "Key generation", seconds);
    indcpa_keypair(pk, sk);

    TIME_OPERATION_SECONDS({ pke_enc(c, m1, pk); }, "PKE encrypt", seconds);
    pke_enc(c, m1, pk);

    TIME_OPERATION_SECONDS({ indcpa_dec(m2, c, sk); }, "PKE decrypt", seconds);
}


int main()
{
    int OK = TRUE;

    OK = pke_test(CRYPTO_ALGNAME, PKE_TEST_ITERATIONS);
    if (OK != TRUE) {
        goto exit;
    }

    PRINT_TIMER_HEADER
    pke_bench(PKE_BENCH_SECONDS);

exit:
    return (OK == TRUE) ? EXIT_SUCCESS : EXIT_FAILURE;
}
