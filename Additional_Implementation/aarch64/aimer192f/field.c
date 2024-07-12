// SPDX-License-Identifier: MIT

#include <stddef.h>
#include <string.h>
#include "field.h"
#include "portable_endian.h"

void GF_to_bytes(const GF in, uint8_t* out)
{
  uint64_t temp = htole64(in[0]);
  size_t num_bytes = sizeof(uint64_t);

  memcpy(out, (uint8_t*)(&temp), num_bytes);
  temp = htole64(in[1]);
  memcpy(out + num_bytes, (uint8_t*)(&temp), num_bytes);
  temp = htole64(in[2]);
  memcpy(out + 2 * num_bytes, (uint8_t*)(&temp), num_bytes);
}

void GF_from_bytes(const uint8_t* in, GF out)
{
  uint64_t temp;
  size_t num_bytes = sizeof(uint64_t);

  memcpy((uint8_t*)(&temp), in, num_bytes);
  out[0] = le64toh(temp);
  memcpy((uint8_t*)(&temp), in + num_bytes, num_bytes);
  out[1] = le64toh(temp);
  memcpy((uint8_t*)(&temp), in + 2 * num_bytes, num_bytes);
  out[2] = le64toh(temp);
}

void GF_transposed_matmul(const GF a, const GF b[AIM2_NUM_BITS_FIELD], GF c)
{
  unsigned int i, j;
  const uint64_t* a_ptr = a;
  const GF* b_ptr = b;

  uint64_t temp_c0 = 0;
  uint64_t temp_c1 = 0;
  uint64_t temp_c2 = 0;
  uint64_t mask;
  for (i = AIM2_NUM_WORDS_FIELD; i; --i, ++a_ptr)
  {
    uint64_t index = *a_ptr;
    for (j = AIM2_NUM_BITS_WORD; j; j -= 4, index >>= 4, b_ptr += 4)
    {
      mask = -(index & 1);
      temp_c0 ^= (b_ptr[0][0] & mask);
      temp_c1 ^= (b_ptr[0][1] & mask);
      temp_c2 ^= (b_ptr[0][2] & mask);

      mask = -((index >> 1) & 1);
      temp_c0 ^= (b_ptr[1][0] & mask);
      temp_c1 ^= (b_ptr[1][1] & mask);
      temp_c2 ^= (b_ptr[1][2] & mask);

      mask = -((index >> 2) & 1);
      temp_c0 ^= (b_ptr[2][0] & mask);
      temp_c1 ^= (b_ptr[2][1] & mask);
      temp_c2 ^= (b_ptr[2][2] & mask);

      mask = -((index >> 3) & 1);
      temp_c0 ^= (b_ptr[3][0] & mask);
      temp_c1 ^= (b_ptr[3][1] & mask);
      temp_c2 ^= (b_ptr[3][2] & mask);
    }
  }
  c[0] = temp_c0;
  c[1] = temp_c1;
  c[2] = temp_c2;
}

void GF_transposed_matmul_add_N(const GF a[AIMER_NUM_MPC_PARTIES],
                                const GF b[AIM2_NUM_BITS_FIELD],
                                GF c[AIMER_NUM_MPC_PARTIES])
{
  unsigned int i, j;

  for (size_t party = 0; party < AIMER_NUM_MPC_PARTIES; party ++)
  {
    const uint64_t* a_ptr = a[party];
    const GF* b_ptr = b;

    uint64_t temp_c0 = 0;
    uint64_t temp_c1 = 0;
    uint64_t temp_c2 = 0;
    uint64_t mask;
    for (i = AIM2_NUM_WORDS_FIELD; i; --i, ++a_ptr)
    {
      uint64_t index = *a_ptr;
      for (j = AIM2_NUM_BITS_WORD; j; j -= 4, index >>= 4, b_ptr += 4)
      {
        mask = -(index & 1);
        temp_c0 ^= (b_ptr[0][0] & mask);
        temp_c1 ^= (b_ptr[0][1] & mask);
        temp_c2 ^= (b_ptr[0][2] & mask);

        mask = -((index >> 1) & 1);
        temp_c0 ^= (b_ptr[1][0] & mask);
        temp_c1 ^= (b_ptr[1][1] & mask);
        temp_c2 ^= (b_ptr[1][2] & mask);

        mask = -((index >> 2) & 1);
        temp_c0 ^= (b_ptr[2][0] & mask);
        temp_c1 ^= (b_ptr[2][1] & mask);
        temp_c2 ^= (b_ptr[2][2] & mask);

        mask = -((index >> 3) & 1);
        temp_c0 ^= (b_ptr[3][0] & mask);
        temp_c1 ^= (b_ptr[3][1] & mask);
        temp_c2 ^= (b_ptr[3][2] & mask);
      }
    }
    c[party][0] ^= temp_c0;
    c[party][1] ^= temp_c1;
    c[party][2] ^= temp_c2;
  }
}
