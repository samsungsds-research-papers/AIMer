// SPDX-License-Identifier: MIT

#include "field.h"
#include "portable_endian.h"

// square the lower 32-bit of the input
#define SQR_LOW(x) \
  sqr_table[((x) >> 28) & 0xf] << 56 | sqr_table[((x) >> 24) & 0xf] << 48 | \
  sqr_table[((x) >> 20) & 0xf] << 40 | sqr_table[((x) >> 16) & 0xf] << 32 | \
  sqr_table[((x) >> 12) & 0xf] << 24 | sqr_table[((x) >>  8) & 0xf] << 16 | \
  sqr_table[((x) >>  4) & 0xf] <<  8 | sqr_table[((x)      ) & 0xf]

// square the upper 32-bit of the input
#define SQR_HIGH(x) \
  sqr_table[((x) >> 60)      ] << 56 | sqr_table[((x) >> 56) & 0xf] << 48 | \
  sqr_table[((x) >> 52) & 0xf] << 40 | sqr_table[((x) >> 48) & 0xf] << 32 | \
  sqr_table[((x) >> 44) & 0xf] << 24 | sqr_table[((x) >> 40) & 0xf] << 16 | \
  sqr_table[((x) >> 36) & 0xf] <<  8 | sqr_table[((x) >> 32) & 0xf]

const uint64_t sqr_table[16] = {0x00, 0x01, 0x04, 0x05, 0x10, 0x11, 0x14, 0x15,
                                0x40, 0x41, 0x44, 0x45, 0x50, 0x51, 0x54, 0x55};

void poly64_mul(const uint64_t a, const uint64_t b, uint64_t *c1, uint64_t *c0);

unsigned GF_getbit(const GF a, unsigned i)
{
  return (a[i / AIM2_NUM_BITS_WORD] >> (i % AIM2_NUM_BITS_WORD)) & 1;
}

void GF_set0(GF a)
{
  for (unsigned i = 0; i < AIM2_NUM_WORDS_FIELD; i++)
  {
    a[i] = 0;
  }
}

void GF_to_bytes(const GF in, uint8_t* out)
{
  uint64_t* temp = (uint64_t*)out;
  for (unsigned i = 0; i < AIM2_NUM_WORDS_FIELD; i++)
  {
    temp[i] = htole64(in[i]);
  }
}

void GF_from_bytes(const uint8_t* in, GF out)
{
  uint64_t* temp = (uint64_t*)in;
  for (unsigned i = 0; i < AIM2_NUM_WORDS_FIELD; i++)
  {
    out[i] = le64toh(temp[i]);
  }
}

void GF_copy(const GF in, GF out)
{
  if (in == out)
    return;

  for (unsigned i = 0; i < AIM2_NUM_WORDS_FIELD; i++)
  {
    out[i] = in[i];
  }
}

void GF_add(const GF a, const GF b, GF c)
{
  for (unsigned i = 0; i < AIM2_NUM_WORDS_FIELD; i++)
  {
    c[i] = a[i] ^ b[i];
  }
}

void poly64_mul(const uint64_t a, const uint64_t b, uint64_t *c1, uint64_t *c0)
{
  uint64_t table[16];
  uint64_t temp, mask, high, low;
  uint64_t top3 = a >> 61;

  table[0] = 0;
  table[1] = a & 0x1fffffffffffffffULL;
  table[2] = table[1] << 1;
  table[4] = table[2] << 1;
  table[8] = table[4] << 1;

  table[3] = table[1] ^ table[2];

  table[5] = table[1] ^ table[4];
  table[6] = table[2] ^ table[4];
  table[7] = table[1] ^ table[6];

  table[9] = table[1] ^ table[8];
  table[10] = table[2] ^ table[8];
  table[11] = table[3] ^ table[8];
  table[12] = table[4] ^ table[8];
  table[13] = table[5] ^ table[8];
  table[14] = table[6] ^ table[8];
  table[15] = table[7] ^ table[8];

  low = table[b & 0xf];
  temp = table[(b >> 4) & 0xf];
  low ^= temp << 4;
  high = temp >> 60;
  temp = table[(b >> 8) & 0xf];
  low ^= temp << 8;
  high ^= temp >> 56;
  temp = table[(b >> 12) & 0xf];
  low ^= temp << 12;
  high ^= temp >> 52;
  temp = table[(b >> 16) & 0xf];
  low ^= temp << 16;
  high ^= temp >> 48;
  temp = table[(b >> 20) & 0xf];
  low ^= temp << 20;
  high ^= temp >> 44;
  temp = table[(b >> 24) & 0xf];
  low ^= temp << 24;
  high ^= temp >> 40;
  temp = table[(b >> 28) & 0xf];
  low ^= temp << 28;
  high ^= temp >> 36;
  temp = table[(b >> 32) & 0xf];
  low ^= temp << 32;
  high ^= temp >> 32;
  temp = table[(b >> 36) & 0xf];
  low ^= temp << 36;
  high ^= temp >> 28;
  temp = table[(b >> 40) & 0xf];
  low ^= temp << 40;
  high ^= temp >> 24;
  temp = table[(b >> 44) & 0xf];
  low ^= temp << 44;
  high ^= temp >> 20;
  temp = table[(b >> 48) & 0xf];
  low ^= temp << 48;
  high ^= temp >> 16;
  temp = table[(b >> 52) & 0xf];
  low ^= temp << 52;
  high ^= temp >> 12;
  temp = table[(b >> 56) & 0xf];
  low ^= temp << 56;
  high ^= temp >> 8;
  temp = table[b >> 60];
  low ^= temp << 60;
  high ^= temp >> 4;

  mask = -(int64_t)(top3 & 0x1);
  low ^= mask & (b << 61);
  high ^= mask & (b >> 3);
  mask = -(int64_t)((top3 >> 1) & 0x1);
  low ^= mask & (b << 62);
  high ^= mask & (b >> 2);
  mask = -(int64_t)((top3 >> 2) & 0x1);
  low ^= mask & (b << 63);
  high ^= mask & (b >> 1);

  *c0 = low;
  *c1 = high;
}

void GF_mul(const GF a, const GF b, GF c)
{
  uint64_t t[2] = {0,};
  uint64_t temp[4] = {0,};

  poly64_mul(a[1], b[1], &temp[3], &temp[2]);
  poly64_mul(a[0], b[0], &temp[1], &temp[0]);

  poly64_mul((a[0] ^ a[1]), (b[0] ^ b[1]), &t[1], &t[0]);
  temp[1] ^= t[0] ^ temp[0] ^ temp[2];
  temp[2] = t[0] ^ t[1] ^ temp[0] ^ temp[1] ^ temp[3];

  t[0] = temp[2] ^ ((temp[3] >> 57) ^ (temp[3] >> 62) ^ (temp[3] >> 63));

  c[1] = temp[1] ^ temp[3];
  c[1] ^= (temp[3] << 7) | (t[0] >> 57);
  c[1] ^= (temp[3] << 2) | (t[0] >> 62);
  c[1] ^= (temp[3] << 1) | (t[0] >> 63);

  c[0] = temp[0] ^ t[0];
  c[0] ^= (t[0] << 7);
  c[0] ^= (t[0] << 2);
  c[0] ^= (t[0] << 1);
}

void GF_mul_add(const GF a, const GF b, GF c)
{
  uint64_t t[2] = {0,};
  uint64_t temp[4] = {0,};

  poly64_mul(a[1], b[1], &temp[3], &temp[2]);
  poly64_mul(a[0], b[0], &temp[1], &temp[0]);

  poly64_mul((a[0] ^ a[1]), (b[0] ^ b[1]), &t[1], &t[0]);
  temp[1] ^= t[0] ^ temp[0] ^ temp[2];
  temp[2] = t[0] ^ t[1] ^ temp[0] ^ temp[1] ^ temp[3];

  t[0] = temp[2] ^ ((temp[3] >> 57) ^ (temp[3] >> 62) ^ (temp[3] >> 63));

  c[1] ^= temp[1] ^ temp[3];
  c[1] ^= (temp[3] << 7) | (t[0] >> 57);
  c[1] ^= (temp[3] << 2) | (t[0] >> 62);
  c[1] ^= (temp[3] << 1) | (t[0] >> 63);

  c[0] ^= temp[0] ^ t[0];
  c[0] ^= (t[0] << 7);
  c[0] ^= (t[0] << 2);
  c[0] ^= (t[0] << 1);
}

void GF_sqr(const GF a, GF c)
{
  uint64_t t = 0;
  uint64_t temp[4] = {0,};

  temp[0] = SQR_LOW(a[0]);
  temp[1] = SQR_HIGH(a[0]);
  temp[2] = SQR_LOW(a[1]);
  temp[3] = SQR_HIGH(a[1]);

  t = temp[2] ^ ((temp[3] >> 57) ^ (temp[3] >> 62) ^ (temp[3] >> 63));

  c[1] = temp[1] ^ temp[3];
  c[1] ^= (temp[3] << 7) | (t >> 57);
  c[1] ^= (temp[3] << 2) | (t >> 62);
  c[1] ^= (temp[3] << 1) | (t >> 63);

  c[0] = temp[0] ^ t;
  c[0] ^= (t << 7);
  c[0] ^= (t << 2);
  c[0] ^= (t << 1);
}

void GF_exp(const GF in, GF out, const uint64_t* exp)
{
  GF in_;
  GF_copy(in, in_);
  GF_set0(out);
  out[0] = 1;
  for (unsigned i = 0; i < AIM2_NUM_WORDS_FIELD; i++)
  {
    uint64_t e = exp[i];
    for (unsigned j = 0; j < AIM2_NUM_BITS_WORD; j++, e >>= 1)
    {
      if (e & 1)
      {
        GF_mul(out, in_, out);
      }
      GF_sqr(in_, in_);
    }
  }
}

void GF_exp_power_of_2(const GF in, GF out, const unsigned exp_power_of_2)
{
  GF_copy(in, out);
  for (unsigned i = 0; i < exp_power_of_2; i++)
  {
    GF_sqr(out, out);
  }
}

// c = sum_i a[i] * b[i]
void GF_transposed_matmul(const GF a, const GF b[AIM2_NUM_BITS_FIELD], GF c)
{
  GF tmp;
  GF_set0(tmp);
  GF_transposed_matmul_add(a, b, tmp);
  GF_copy(tmp, c);
}

// c += sum_i a[i] * b[i]
void GF_transposed_matmul_add(const GF a, const GF b[AIM2_NUM_BITS_FIELD], GF c)
{
  GF a_;
  GF_copy(a, a_);
  for (unsigned i = 0; i < AIM2_NUM_BITS_FIELD; i++)
  {
    if (GF_getbit(a_, i))
    {
      GF_add(c, b[i], c);
    }
  }
}
