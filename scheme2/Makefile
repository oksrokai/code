CC ?= /usr/bin/cc
CFLAGS += -Wall -Wextra -Wpedantic -Wmissing-prototypes -Wredundant-decls \
  -Wshadow -Wpointer-arith -mavx2 -mbmi2 -mpopcnt -maes \
  -march=native -mtune=native -O3 -fomit-frame-pointer -fno-stack-check
NISTFLAGS += -Wno-unused-result -mavx2 -mbmi -mpopcnt -maes \
  -march=native -mtune=native -O3
RM = /bin/rm

SOURCES = kem.c indcpa.c polyvec.c poly.c fq.S shuffle.S ntt.S invntt.S \
  basemul.S consts.c rejsample.c cbd.c verify.c con.S rec.S poly_pack.S poly_unpack.S basemul_multi.S 
SOURCESKECCAK   = $(SOURCES) fips202.c fips202x4.c symmetric-shake.c \
  keccak4x/KeccakP-1600-times4-SIMD256.o
SOURCESNINETIES = $(SOURCES) aes256ctr.c
HEADERS = params.h api.h indcpa.h polyvec.h poly.h reduce.h fq.inc shuffle.inc \
  ntt.h rejsample.h cbd.h verify.h symmetric.h randombytes.h consts.h align.h
HEADERSKECCAK   = $(HEADERS) fips202.h fips202x4.h

all: test_OKAI512_v1 test_OKAI768_v1 test_OKAI1024_v1 speed512_v1 speed768_v1 speed1024_v1 test_OKAI512_v2 test_OKAI768_v2 test_OKAI1024_v2 speed512_v2 speed768_v2 speed1024_v2

test_OKAI512_v1: $(SOURCESKECCAK) $(HEADERSKECCAK) test_OKAI.c randombytes.c
	$(CC) $(CFLAGS) -DOKAI_K=2 -DOKAI_N=256 $(SOURCESKECCAK) randombytes.c test_OKAI.c -o test_OKAI512_v1

test_OKAI768_v1: $(SOURCESKECCAK) $(HEADERSKECCAK) test_OKAI.c randombytes.c
	$(CC) $(CFLAGS) -DOKAI_K=3 -DOKAI_N=256 $(SOURCESKECCAK) randombytes.c test_OKAI.c -o test_OKAI768_v1

test_OKAI1024_v1: $(SOURCESKECCAK) $(HEADERSKECCAK) test_OKAI.c randombytes.c
	$(CC) $(CFLAGS) -DOKAI_K=2 -DOKAI_N=512 $(SOURCESKECCAK) randombytes.c test_OKAI.c -o test_OKAI1024_v1

speed512_v1: $(SOURCESKECCAK) $(HEADERSKECCAK) cpucycles.h cpucycles.c speed_print.h speed_print.c test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DOKAI_K=2 -DOKAI_N=256 $(SOURCESKECCAK) randombytes.c cpucycles.c speed_print.c test_speed.c -o speed512_v1

speed768_v1: $(SOURCESKECCAK) $(HEADERSKECCAK) cpucycles.h cpucycles.c speed_print.h speed_print.c test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DOKAI_K=3 -DOKAI_N=256 $(SOURCESKECCAK) randombytes.c cpucycles.c speed_print.c test_speed.c -o speed768_v1

speed1024_v1: $(SOURCESKECCAK) $(HEADERSKECCAK) cpucycles.h cpucycles.c speed_print.h speed_print.c test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DOKAI_K=2 -DOKAI_N=512 $(SOURCESKECCAK) randombytes.c cpucycles.c speed_print.c test_speed.c -o speed1024_v1

test_OKAI512_v2: $(SOURCESKECCAK) $(HEADERSKECCAK) test_OKAI.c randombytes.c
	$(CC) $(CFLAGS) -DOKAI_K=2 -DOKAI_N=256 -DFO=1 $(SOURCESKECCAK) randombytes.c test_OKAI.c -o test_OKAI512_v2

test_OKAI768_v2: $(SOURCESKECCAK) $(HEADERSKECCAK) test_OKAI.c randombytes.c
	$(CC) $(CFLAGS) -DOKAI_K=3 -DOKAI_N=256 -DFO=1 $(SOURCESKECCAK) randombytes.c test_OKAI.c -o test_OKAI768_v2

test_OKAI1024_v2: $(SOURCESKECCAK) $(HEADERSKECCAK) test_OKAI.c randombytes.c
	$(CC) $(CFLAGS) -DOKAI_K=2 -DOKAI_N=512 -DFO=1 $(SOURCESKECCAK) randombytes.c test_OKAI.c -o test_OKAI1024_v2

speed512_v2: $(SOURCESKECCAK) $(HEADERSKECCAK) cpucycles.h cpucycles.c speed_print.h speed_print.c test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DOKAI_K=2 -DOKAI_N=256 -DFO=1 $(SOURCESKECCAK) randombytes.c cpucycles.c speed_print.c test_speed.c -o speed512_v2

speed768_v2: $(SOURCESKECCAK) $(HEADERSKECCAK) cpucycles.h cpucycles.c speed_print.h speed_print.c test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DOKAI_K=3 -DOKAI_N=256 -DFO=1 $(SOURCESKECCAK) randombytes.c cpucycles.c speed_print.c test_speed.c -o speed768_v2

speed1024_v2: $(SOURCESKECCAK) $(HEADERSKECCAK) cpucycles.h cpucycles.c speed_print.h speed_print.c test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DOKAI_K=2 -DOKAI_N=512 -DFO=1 $(SOURCESKECCAK) randombytes.c cpucycles.c speed_print.c test_speed.c -o speed1024_v2

.PHONY: clean

clean:
	$(RM) -rf test_OKAI512_v1 test_OKAI768_v1 test_OKAI1024_v1 speed512_v1 speed768_v1 speed1024_v1
	$(RM) -rf test_OKAI512_v2 test_OKAI768_v2 test_OKAI1024_v2 speed512_v2 speed768_v2 speed1024_v2
	$(RM) -rf keccak4x/KeccakP-1600-times4-SIMD256.o