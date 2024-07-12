# SPDX-License-Identifier: MIT

CC = aarch64-linux-gnu-gcc
CFLAGS += -I. -O3 -Wall -Wextra -march=armv8-a+crypto -fomit-frame-pointer
NISTFLAGS = -Wno-sign-compare -Wno-unused-but-set-variable -Wno-unused-parameter -Wno-unused-result
SHAKE_PATH = shake
SHAKE_LIB = libshake.a
LDFLAGS = $(SHAKE_PATH)/$(SHAKE_LIB)

EXECUTABLE_TESTAIM2 = tests/test_aim2
EXECUTABLE_TESTSIGN = tests/test_sign
EXECUTABLE_KAT      = PQCgenKAT_sign

.PHONY: all

all: $(SHAKE_LIB) $(EXECUTABLE_TESTAIM2) $(EXECUTABLE_TESTSIGN) $(EXECUTABLE_KAT)

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

$(SHAKE_LIB):
	$(MAKE) -C $(SHAKE_PATH)

$(EXECUTABLE_TESTAIM2): $(EXECUTABLE_TESTAIM2).c field.c field.s aim2.c hash.c
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

$(EXECUTABLE_TESTSIGN): $(EXECUTABLE_TESTSIGN).c field.c field.s aim2.c rng.c aes.c hash.c tree.c aimer.c api.c
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

$(EXECUTABLE_KAT): $(EXECUTABLE_KAT).c field.c field.s aim2.c rng.c aes.c hash.c tree.c aimer.c api.c
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

clean:
	rm -f $(wildcard *.o) $(EXECUTABLE_TESTAIM2) $(EXECUTABLE_TESTSIGN) $(EXECUTABLE_KAT) PQCsignKAT_*.rsp PQCsignKAT_*.req
	$(MAKE) -C $(SHAKE_PATH) clean
