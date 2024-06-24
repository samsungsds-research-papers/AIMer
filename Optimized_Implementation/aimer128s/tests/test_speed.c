// SPDX-License-Identifier: MIT

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../api.h"
#include "../common/cpucycles.h"
#include "../common/speed_print.h"
#include "../common/rng.h"

#define NTESTS 10000

uint64_t t[NTESTS];

int main()
{
  int i;
  int ret = 0;

  size_t mlen = 59;
  size_t smlen = 0;
  uint8_t message[mlen];
  uint8_t m[mlen];
  uint8_t pk[CRYPTO_PUBLICKEYBYTES];
  uint8_t sk[CRYPTO_SECRETKEYBYTES];
  uint8_t sm[mlen + CRYPTO_BYTES];

  uint64_t time_begin, time_end;

  randombytes((unsigned char*)message, 59);

  for(i = 0; i < NTESTS; i++)
  {
    time_begin = cpucycles();
    ret = crypto_sign_keypair(pk, sk);
    time_end   = cpucycles();
    t[i] = time_end - time_begin;
  }
  print_results("aimer_keypair: ", t, NTESTS);

  for(i = 0; i < NTESTS; i++)
  {
    time_begin = cpucycles();
    ret = crypto_sign(sm, &smlen, message, mlen, sk);
    time_end   = cpucycles();
    t[i] = time_end - time_begin;
  }
  print_results("aimer_sign   : ", t, NTESTS);

  for(i = 0; i < NTESTS; i++)
  {
    time_begin = cpucycles();
    ret = crypto_sign_open(m, &mlen, sm, smlen, pk);
    time_end   = cpucycles();
    t[i] = time_end - time_begin;
  }
  print_results("aimer_verify : ", t, NTESTS);

  return ret;
}
