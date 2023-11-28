CC=/usr/bin/gcc
LDFLAGS=-L/usr/local/opt/openssl@1.1/lib -I/usr/local/opt/openssl@1.1/include -lcrypto -lm

SOURCES= cbd.c fips202.c indcpa.c kem.c ntt.c poly.c polyvec.c reduce.c rng.c symmetric-fips202.c verify.c
HEADERS= api.h cbd.h fips202.h indcpa.h ntt.h params.h poly.h polyvec.h reduce.h rng.h symmetric.h verify.h

all: testbench_kem maintest #testbench_pke PQCgenKAT_kem

# PQCgenKAT_kem: $(HEADERS) $(SOURCES)
# 	$(CC) -o $@ $(SOURCES) PQCgenKAT_kem.c $(LDFLAGS)

testbench_kem: $(HEADERS) $(SOURCES)
	$(CC) -o $@ $(SOURCES) testbench_kem.c $(LDFLAGS)

maintest: $(HEADERS) $(SOURCES)
	$(CC) -o $@ $(SOURCES) maintest.c $(LDFLAGS)

# testbench_pke: $(HEADERS) $(SOURCES)
# 	$(CC) -o $@ $(SOURCES) testbench_pke.c $(LDFLAGS)

.PHONY: clean

clean:
	rm -rf testbench_kem maintest
