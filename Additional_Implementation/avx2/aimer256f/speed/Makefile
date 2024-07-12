# SPDX-License-Identifier: MIT

CC = gcc
CFLAGS += -I. -O3 -Wall -Wextra -march=native -fomit-frame-pointer
NISTFLAGS = -Wno-sign-compare -Wno-unused-but-set-variable -Wno-unused-parameter -Wno-unused-result
AVX2FLAGS = -mavx2 -mpclmul
SHAKE_PATH = ../shake
SHAKE_LIB = libshake.a
LDFLAGS = $(SHAKE_PATH)/$(SHAKE_LIB)

EXECUTABLE_TESTSPEED = test_speed
EXECUTABLE_TESTSIGN = test_sign

.PHONY: all

all: $(SHAKE_LIB) $(EXECUTABLE_TESTSPEED) $(EXECUTABLE_TESTSIGN)

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

$(SHAKE_LIB):
	$(MAKE) -C $(SHAKE_PATH)

$(EXECUTABLE_TESTSPEED): $(EXECUTABLE_TESTSPEED).c speed_print.c cpucycles.c ../field.c ../aim2.c ../rng.c ../aes.c ../hash.c ../tree.c ../aimer.c ../api.c
	$(CC) $(CFLAGS) $(AVX2FLAGS) $^ $(LDFLAGS) -o $@

$(EXECUTABLE_TESTSIGN): $(EXECUTABLE_TESTSIGN).c ../field.c ../aim2.c ../rng.c ../aes.c ../hash.c ../tree.c ../aimer.c ../api.c
	$(CC) $(CFLAGS) $(AVX2FLAGS) $^ $(LDFLAGS) -o $@

clean:
	rm -f $(wildcard *.o) $(EXECUTABLE_TESTSPEED) $(EXECUTABLE_TESTSIGN)
	$(MAKE) -C $(SHAKE_PATH) clean
