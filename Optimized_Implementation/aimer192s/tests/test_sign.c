// SPDX-License-Identifier: MIT

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include "../api.h"
#include "../common/rng.h"

#define MLEN 59
#define NTESTS 5

int main(void)
{
  size_t i, j;
  size_t mlen, smlen;
  uint8_t b;
  uint8_t m[MLEN + CRYPTO_BYTES];
  uint8_t m2[MLEN + CRYPTO_BYTES];
  uint8_t sm[MLEN + CRYPTO_BYTES];
  uint8_t pk[CRYPTO_PUBLICKEYBYTES];
  uint8_t sk[CRYPTO_SECRETKEYBYTES];

  for (i = 0; i < NTESTS; ++i)
  {
    int ret = 0;

    randombytes(m, MLEN);

    crypto_sign_keypair(pk, sk);
    crypto_sign(sm, &smlen, m, MLEN, sk);

    ret = crypto_sign_open(m2, &mlen, sm, smlen, pk);

    if (ret)
    {
      printf("Verification failed\n");
      return -1;
    }
    if (smlen != MLEN + CRYPTO_BYTES)
    {
      printf("Signed message lengths wrong\n");
      return -1;
    }
    if (mlen != MLEN)
    {
      printf("Message lengths wrong\n");
      return -1;
    }
    for (j = 0; j < MLEN; ++j)
    {
      if (m2[j] != m[j])
      {
        printf("Messages don't match\n");
        return -1;
      }
    }

    randombytes((uint8_t *)&j, sizeof(j));
    do
    {
      randombytes(&b, 1);
    }
    while (!b);
    sm[j % (MLEN + CRYPTO_BYTES)] += b;
    ret = crypto_sign_open(m2, &mlen, sm, smlen, pk);
    if (!ret)
    {
      printf("Trivial forgeries possible\n");
      return -1;
    }
  }

  printf("CRYPTO_PUBLICKEYBYTES = %d\n", CRYPTO_PUBLICKEYBYTES);
  printf("CRYPTO_SECRETKEYBYTES = %d\n", CRYPTO_SECRETKEYBYTES);
  printf("CRYPTO_BYTES = %d\n", CRYPTO_BYTES);

  return 0;
}
