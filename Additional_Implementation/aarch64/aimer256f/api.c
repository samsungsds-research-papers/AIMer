// SPDX-License-Identifier: MIT

#include <string.h>
#include "api.h"
#include "aim2.h"
#include "aimer.h"
#include "params.h"
#include "rng.h"

int
crypto_sign_keypair(unsigned char *pk, unsigned char *sk)
{
  int ret = 0;

  if (!pk || !sk)
  {
    return -1;
  }

  randombytes(sk, AIM2_NUM_BYTES_FIELD);
  randombytes(pk, AIM2_IV_SIZE);

  aim2(sk, pk, pk + AIM2_IV_SIZE);
  memcpy(sk + AIM2_NUM_BYTES_FIELD, pk, AIM2_IV_SIZE + AIM2_NUM_BYTES_FIELD);

  return ret;
}

int
crypto_sign(unsigned char *sm, unsigned long long *smlen,
            const unsigned char *m, unsigned long long mlen,
            const unsigned char *sk)
{
  int ret = 0;
  uint8_t sig[CRYPTO_BYTES];

  ret = aimer_sign(sk + AIM2_NUM_BYTES_FIELD, sk, m, (size_t)mlen, sig);
  if (ret != 0)
  {
    return ret;
  }

  *smlen = mlen + CRYPTO_BYTES;
  memcpy(sm, m, mlen);
  memcpy(sm + mlen, sig, CRYPTO_BYTES);

  return ret;
}

int
crypto_sign_open(unsigned char *m, unsigned long long *mlen,
                 const unsigned char *sm, unsigned long long smlen,
                 const unsigned char *pk)
{
  int ret = 0;

  const size_t message_len = smlen - CRYPTO_BYTES;
  const uint8_t* message = sm;
  const uint8_t* signature = sm + message_len;

  ret = aimer_verify(pk, signature, message, message_len);
  if (ret != 0)
  {
    return ret;
  }

  memcpy(m, message, message_len);
  *mlen = message_len;

  return ret;
}
