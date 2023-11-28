#ifndef PARAMS_H
#define PARAMS_H

#define OKAI_K 3

#define OKAI_N 256
#define OKAI_Q 7681

#define OKAI_NETA 4
#define OKAI_SETA 1

#define OKAI_SYMBYTES 32   /* size in bytes of shared key, hashes, and seeds */

#define OKAI_POLYBYTES              416
#define OKAI_POLYCOMPRESSEDBYTES    128
#define OKAI_POLYVECBYTES           (OKAI_K * OKAI_POLYBYTES)
#define OKAI_POLYVECCOMPRESSEDBYTES (OKAI_K * 288)

#define OKAI_INDCPA_MSGBYTES       OKAI_SYMBYTES
#define OKAI_INDCPA_PUBLICKEYBYTES (OKAI_POLYVECCOMPRESSEDBYTES + OKAI_SYMBYTES)
#define OKAI_INDCPA_SECRETKEYBYTES (OKAI_POLYVECBYTES)
#define OKAI_INDCPA_BYTES          (OKAI_POLYVECCOMPRESSEDBYTES + OKAI_POLYCOMPRESSEDBYTES)

#define OKAI_PUBLICKEYBYTES  (OKAI_INDCPA_PUBLICKEYBYTES)
#define OKAI_SECRETKEYBYTES  (OKAI_INDCPA_SECRETKEYBYTES +  OKAI_INDCPA_PUBLICKEYBYTES + 2*OKAI_SYMBYTES) /* 32 bytes of additional space to save H(pk) */
#define OKAI_CIPHERTEXTBYTES  OKAI_INDCPA_BYTES

#endif
