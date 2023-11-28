#ifndef PARAMS_H
#define PARAMS_H

#define OSKR_K 2 /* Change this for different security strengths */

/* Don't change parameters below this line */

#define OSKR_N 256
#define OSKR_Q 3329

#define OSKR_ETA1 3
#define OSKR_ETA2 2

#define OSKR_SYMBYTES 32   /* size in bytes of hashes, and seeds */
#define OSKR_SSBYTES  32   /* size in bytes of shared key */

#define OSKR_POLYBYTES              384
#define OSKR_POLYVECBYTES           (OSKR_K * OSKR_POLYBYTES)

#define OSKR_POLYCOMPRESSEDBYTES    128
#define OSKR_POLYVECCOMPRESSEDBYTES (OSKR_K * 320)

#define OSKR_INDCPA_MSGBYTES       OSKR_SYMBYTES
#define OSKR_INDCPA_PUBLICKEYBYTES (OSKR_POLYVECBYTES + OSKR_SYMBYTES)
#define OSKR_INDCPA_SECRETKEYBYTES (OSKR_POLYVECBYTES)
#define OSKR_INDCPA_BYTES          (OSKR_POLYVECCOMPRESSEDBYTES + OSKR_POLYCOMPRESSEDBYTES)

#define OSKR_PUBLICKEYBYTES  (OSKR_INDCPA_PUBLICKEYBYTES)
#define OSKR_SECRETKEYBYTES  (OSKR_INDCPA_SECRETKEYBYTES +  OSKR_INDCPA_PUBLICKEYBYTES + 2*OSKR_SYMBYTES) /* 32 bytes of additional space to save H(pk) */
#define OSKR_CIPHERTEXTBYTES  OSKR_INDCPA_BYTES

#endif
