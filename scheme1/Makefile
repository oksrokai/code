CC ?= /usr/bin/cc
CFLAGS += -Wall -Wextra -Wpedantic -Wmissing-prototypes -Wredundant-decls \
  -Wshadow -Wpointer-arith -mavx2 -mbmi2 -mpopcnt -maes \
  -march=native -mtune=native -O3 -fomit-frame-pointer
NISTFLAGS += -Wno-unused-result -mavx2 -mbmi2 -mpopcnt -maes \
  -march=native -mtune=native -O3 -fomit-frame-pointer
RM = /bin/rm

SOURCES = kem.c indcpa.c polyvec.c poly.c fq.S shuffle.S ntt.S invntt.S \
  basemul.S basemul_multi.S consts.c rejsample.c cbd.c verify.c con45.S rec45.S poly_pack.S poly_unpack.S
SOURCESKECCAK   = $(SOURCES) fips202.c fips202x4.c symmetric-shake.c \
  keccak4x/KeccakP-1600-times4-SIMD256.o
SOURCESNINETIES = $(SOURCES) aes256ctr.c
HEADERS = params.h align.h kem.h indcpa.h polyvec.h poly.h reduce.h fq.inc shuffle.inc \
  ntt.h consts.h rejsample.h cbd.h verify.h symmetric.h randombytes.h
HEADERSKECCAK   = $(HEADERS) fips202.h fips202x4.h
HEADERSNINETIES = $(HEADERS) aes256ctr.h

.PHONY: all clean

all:test_OSKR512_v1 test_OSKR768_v1 test_OSKR1024_v1 speed512_v1 speed768_v1 speed1024_v1 test_OSKR512_v2 test_OSKR768_v2 test_OSKR1024_v2 speed512_v2 speed768_v2 speed1024_v2

test_OSKR512_v1: $(SOURCESKECCAK) $(HEADERSKECCAK) test_OSKR.c randombytes.c
	$(CC) $(CFLAGS) -DOSKR_K=2 -DOSKR_N=256 $(SOURCESKECCAK) randombytes.c test_OSKR.c -o test_OSKR512_v1

test_OSKR768_v1: $(SOURCESKECCAK) $(HEADERSKECCAK) test_OSKR.c randombytes.c
	$(CC) $(CFLAGS) -DOSKR_K=3 -DOSKR_N=256 $(SOURCESKECCAK) randombytes.c test_OSKR.c -o test_OSKR768_v1

test_OSKR1024_v1: $(SOURCESKECCAK) $(HEADERSKECCAK) test_OSKR.c randombytes.c
	$(CC) $(CFLAGS) -DOSKR_K=2 -DOSKR_N=512 $(SOURCESKECCAK) randombytes.c test_OSKR.c -o test_OSKR1024_v1

speed512_v1: $(SOURCESKECCAK) $(HEADERSKECCAK) cpucycles.h cpucycles.c speed_print.h speed_print.c test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DOSKR_K=2 -DOSKR_N=256 $(SOURCESKECCAK) randombytes.c cpucycles.c speed_print.c test_speed.c -o speed512_v1

speed768_v1: $(SOURCESKECCAK) $(HEADERSKECCAK) cpucycles.h cpucycles.c speed_print.h speed_print.c test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DOSKR_K=3 -DOSKR_N=256 $(SOURCESKECCAK) randombytes.c cpucycles.c speed_print.c test_speed.c -o speed768_v1

speed1024_v1: $(SOURCESKECCAK) $(HEADERSKECCAK) cpucycles.h cpucycles.c speed_print.h speed_print.c test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DOSKR_K=2 -DOSKR_N=512 $(SOURCESKECCAK) randombytes.c cpucycles.c speed_print.c test_speed.c -o speed1024_v1

test_OSKR512_v2: $(SOURCESKECCAK) $(HEADERSKECCAK) test_OSKR.c randombytes.c
	$(CC) $(CFLAGS) -DOSKR_K=2 -DOSKR_N=256 -DFO=1 $(SOURCESKECCAK) randombytes.c test_OSKR.c -o test_OSKR512_v2

test_OSKR768_v2: $(SOURCESKECCAK) $(HEADERSKECCAK) test_OSKR.c randombytes.c
	$(CC) $(CFLAGS) -DOSKR_K=3 -DOSKR_N=256 -DFO=1 $(SOURCESKECCAK) randombytes.c test_OSKR.c -o test_OSKR768_v2

test_OSKR1024_v2: $(SOURCESKECCAK) $(HEADERSKECCAK) test_OSKR.c randombytes.c
	$(CC) $(CFLAGS) -DOSKR_K=2 -DOSKR_N=512 -DFO=1 $(SOURCESKECCAK) randombytes.c test_OSKR.c -o test_OSKR1024_v2

speed512_v2: $(SOURCESKECCAK) $(HEADERSKECCAK) cpucycles.h cpucycles.c speed_print.h speed_print.c test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DOSKR_K=2 -DOSKR_N=256 -DFO=1 $(SOURCESKECCAK) randombytes.c cpucycles.c speed_print.c test_speed.c -o speed512_v2

speed768_v2: $(SOURCESKECCAK) $(HEADERSKECCAK) cpucycles.h cpucycles.c speed_print.h speed_print.c test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DOSKR_K=3 -DOSKR_N=256 -DFO=1 $(SOURCESKECCAK) randombytes.c cpucycles.c speed_print.c test_speed.c -o speed768_v2

speed1024_v2: $(SOURCESKECCAK) $(HEADERSKECCAK) cpucycles.h cpucycles.c speed_print.h speed_print.c test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DOSKR_K=2 -DOSKR_N=512 -DFO=1 $(SOURCESKECCAK) randombytes.c cpucycles.c speed_print.c test_speed.c -o speed1024_v2

clean:
	-$(RM) -rf *.o 
	-$(RM) -rf test_OSKR512_v1
	-$(RM) -rf test_OSKR768_v1
	-$(RM) -rf test_OSKR1024_v1
	-$(RM) -rf speed512_v1
	-$(RM) -rf speed768_v1
	-$(RM) -rf speed1024_v1
	-$(RM) -rf test_OSKR512_v2
	-$(RM) -rf test_OSKR768_v2
	-$(RM) -rf test_OSKR1024_v2
	-$(RM) -rf speed512_v2
	-$(RM) -rf speed768_v2
	-$(RM) -rf speed1024_v2