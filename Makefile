OPENCM3DIR  = ./libopencm3
OPENCM3NAME = opencm3_stm32f4
OPENCM3FILE = $(OPENCM3DIR)/lib/lib$(OPENCM3NAME).a
LDSCRIPT    = common/stm32f407.ld

PREFIX     ?= arm-none-eabi
CC          = $(PREFIX)-gcc
LD          = $(PREFIX)-gcc
OBJCOPY     = $(PREFIX)-objcopy

ARCH_FLAGS  = -mthumb -mcpu=cortex-m4 -mfloat-abi=hard -mfpu=fpv4-sp-d16
DEFINES     = -DSTM32F4

CFLAGS     += -O3 -std=gnu99 \
              -Wall -Wextra -Wimplicit-function-declaration \
              -Wredundant-decls -Wmissing-prototypes -Wstrict-prototypes \
              -Wundef -Wshadow \
              -I$(OPENCM3DIR)/include \
              -fno-common $(ARCH_FLAGS) -MD $(DEFINES)

CC_HOST    = gcc
LD_HOST    = gcc

CFLAGS_HOST = -O3 -Wall -Wextra -Wpedantic
LDFLAGS_HOST = -lm

COMMONSOURCES=common/fips202.c common/sp800-185.c common/nistseedexpander.c
COMMONSOURCES_HOST=$(COMMONSOURCES) common/keccakf1600.c common/host/aes.c common/host/sha2.c
COMMONSOURCES_M4=$(COMMONSOURCES) common/keccakf1600.S common/aes.c common/aes.S common/sha2.c common/crypto_hashblocks_sha512.c common/crypto_hashblocks_sha512_inner32.s

COMMONINCLUDES=-I"common"
COMMONINCLUDES_M4=$(COMMONINCLUDES)

RANDOMBYTES_M4=common/randombytes.c

DEST_HOST=bin/host
DEST=bin

TARGET_NAME = $(shell echo $(IMPLEMENTATION_PATH) | sed 's@/@_@g')
IMPLEMENTATION_SOURCES = $(wildcard $(IMPLEMENTATION_PATH)/*.c) $(wildcard $(IMPLEMENTATION_PATH)/*.s) $(wildcard $(IMPLEMENTATION_PATH)/*.S)
IMPLEMENTATION_HEADERS = $(IMPLEMENTATION_PATH)/*.h

LDFLAGS    += --static -Wl,--start-group -lc -lgcc -lnosys -Wl,--end-group \
              -T$(LDSCRIPT) -nostartfiles -Wl,--gc-sections \
               $(ARCH_FLAGS) -L$(OPENCM3DIR)/lib -lm -l$(OPENCM3NAME)

.PHONY: all
.SECONDARY:
all:

$(DEST_HOST)/%_testvectors: $(COMMONSOURCES_HOST) $(IMPLEMENTATION_SOURCES) $(IMPLEMENTATION_HEADERS)
	mkdir -p $(DEST_HOST)
	$(CC_HOST) -o $@ \
		$(CFLAGS_HOST) \
		tests/testvectors-host.c \
		$(COMMONSOURCES_HOST) \
		$(IMPLEMENTATION_SOURCES) \
		-I$(IMPLEMENTATION_PATH) \
		$(COMMONINCLUDES) \
		$(LDFLAGS_HOST)

$(DEST)/%.bin: elf/%.elf
	mkdir -p $(DEST)
	$(OBJCOPY) -Obinary $^ $@


# pattern rules, intended to match % to the type of test (i.e. test, speed, stack)
# note that this excludes testvectors, as that is a special case that provides its own randombytes
# TODO use notrandombytes more generically rather than included in testvectors.c
elf/$(TARGET_NAME)_%.elf: tests/%.c $(COMMONSOURCES_M4) $(RANDOMBYTES_M4) $(IMPLEMENTATION_SOURCES) $(IMPLEMENTATION_HEADERS) $(OPENCM3FILE) common/hal-stm32f4.c
	mkdir -p elf
	$(CC) -o $@ $(CFLAGS) \
		$< $(COMMONSOURCES_M4) $(RANDOMBYTES_M4) $(IMPLEMENTATION_SOURCES) common/hal-stm32f4.c \
		-I$(IMPLEMENTATION_PATH) $(COMMONINCLUDES_M4) $(LDFLAGS)

elf/$(TARGET_NAME)_testvectors.elf: tests/testvectors.c $(COMMONSOURCES_M4) $(IMPLEMENTATION_SOURCES) $(IMPLEMENTATION_HEADERS) $(OPENCM3FILE) common/hal-stm32f4.c
	mkdir -p elf
	$(CC) -o $@ $(CFLAGS) \
		$< $(COMMONSOURCES_M4) $(IMPLEMENTATION_SOURCES) common/hal-stm32f4.c \
		-I$(IMPLEMENTATION_PATH) $(COMMONINCLUDES_M4) $(LDFLAGS)

elf/$(TARGET_NAME)_hashing.elf: tests/hashing.c $(COMMONSOURCES_M4) $(IMPLEMENTATION_SOURCES) $(IMPLEMENTATION_HEADERS) $(OPENCM3FILE) common/hal-stm32f4.c
	mkdir -p elf
	$(CC) -o $@ $(CFLAGS) -DPROFILE_HASHING \
		$< $(COMMONSOURCES_M4) $(RANDOMBYTES_M4) $(IMPLEMENTATION_SOURCES) common/hal-stm32f4.c \
		-I$(IMPLEMENTATION_PATH) $(COMMONINCLUDES_M4) $(LDFLAGS)

obj/$(TARGET_NAME)_%.o: $(IMPLEMENTATION_PATH)/%.c $(IMPLEMENTATION_HEADERS)
	mkdir -p obj
	$(CC) -o $@ -c $(CFLAGS) \
		-I$(IMPLEMENTATION_PATH) $(COMMONINCLUDES_M4) $<

obj/$(TARGET_NAME)_%.o: $(IMPLEMENTATION_PATH)/%.s $(IMPLEMENTATION_HEADERS)
	mkdir -p obj
	$(CC) -o $@ -c $(CFLAGS) \
		-I$(IMPLEMENTATION_PATH) $(COMMONINCLUDES_M4) $<

obj/$(TARGET_NAME)_%.o: $(IMPLEMENTATION_PATH)/%.S $(IMPLEMENTATION_HEADERS)
	mkdir -p obj
	$(CC) -o $@ -c $(CFLAGS) \
		-I$(IMPLEMENTATION_PATH) $(COMMONINCLUDES_M4) $<

$(OPENCM3FILE):
	@if [ ! "`ls -A $(OPENCM3_DIR)`" ] ; then \
		printf "######## ERROR ########\n"; \
		printf "\tlibopencm3 is not initialized.\n"; \
		printf "\tPlease run (in the root directory):\n"; \
		printf "\t$$ git submodule init\n"; \
		printf "\t$$ git submodule update\n"; \
		printf "\tbefore running make.\n"; \
		printf "######## ERROR ########\n"; \
		exit 1; \
		fi
	make -C $(OPENCM3DIR)

.PHONY: clean libclean

clean:
	rm -rf elf/
	rm -rf bin/
	rm -rf obj/

libclean:
	make -C $(OPENCM3DIR) clean
