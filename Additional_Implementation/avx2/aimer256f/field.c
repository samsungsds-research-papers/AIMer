// SPDX-License-Identifier: MIT

#include <immintrin.h>
#include <emmintrin.h>
#include "field.h"
#include "portable_endian.h"
#include "params.h"

void GF_to_bytes(const GF in, uint8_t* out)
{
  ((uint64_t*)out)[0] = htole64(in[0]);
  ((uint64_t*)out)[1] = htole64(in[1]);
  ((uint64_t*)out)[2] = htole64(in[2]);
  ((uint64_t*)out)[3] = htole64(in[3]);
}

void GF_from_bytes(const uint8_t* in, GF out)
{
  out[0] = le64toh(((uint64_t*)in)[0]);
  out[1] = le64toh(((uint64_t*)in)[1]);
  out[2] = le64toh(((uint64_t*)in)[2]);
  out[3] = le64toh(((uint64_t*)in)[3]);
}

void GF_set0(GF a)
{
  a[0] = 0;
  a[1] = 0;
  a[2] = 0;
  a[3] = 0;
}

void GF_copy(const GF in, GF out)
{
  out[0] = in[0];
  out[1] = in[1];
  out[2] = in[2];
  out[3] = in[3];
}

void GF_add(const GF a, const GF b, GF c)
{
  c[0] = a[0] ^ b[0];
  c[1] = a[1] ^ b[1];
  c[2] = a[2] ^ b[2];
  c[3] = a[3] ^ b[3];
}

void GF_mul(const GF a, const GF b, GF c)
{
  __m128i x[2], y[2], z[4], t[4];
  __m128i irr = _mm_set_epi64x(0x0, 0x425);

  // polynomial multiplication
  x[0] = _mm_load_si128((const __m128i*)&a[0]); // a0 a1
  x[1] = _mm_load_si128((const __m128i*)&a[2]); // a2 a3
  y[0] = _mm_load_si128((const __m128i*)&b[0]); // b0 b1
  y[1] = _mm_load_si128((const __m128i*)&b[2]); // b2 b3

  t[0] = _mm_clmulepi64_si128(x[0], y[0], 0x10);
  t[1] = _mm_clmulepi64_si128(x[0], y[0], 0x01);
  z[0] = _mm_clmulepi64_si128(x[0], y[0], 0x00);
  z[1] = _mm_clmulepi64_si128(x[0], y[0], 0x11);

  t[0] = _mm_xor_si128(t[0], t[1]);
  t[1] = _mm_srli_si128(t[0], 8);
  t[0] = _mm_slli_si128(t[0], 8);
  z[0] = _mm_xor_si128(z[0], t[0]);
  z[1] = _mm_xor_si128(z[1], t[1]);

  t[2] = _mm_clmulepi64_si128(x[1], y[1], 0x10);
  t[3] = _mm_clmulepi64_si128(x[1], y[1], 0x01);
  z[2] = _mm_clmulepi64_si128(x[1], y[1], 0x00);
  z[3] = _mm_clmulepi64_si128(x[1], y[1], 0x11);

  t[2] = _mm_xor_si128(t[2], t[3]);
  t[3] = _mm_srli_si128(t[2], 8);
  t[2] = _mm_slli_si128(t[2], 8);
  z[2] = _mm_xor_si128(z[2], t[2]);
  z[3] = _mm_xor_si128(z[3], t[3]);

  // Start Karat
  x[0] = _mm_xor_si128(x[0], x[1]);
  y[0] = _mm_xor_si128(y[0], y[1]);

  t[0] = _mm_clmulepi64_si128(x[0], y[0], 0x00);
  t[1] = _mm_clmulepi64_si128(x[0], y[0], 0x11);
  t[2] = _mm_clmulepi64_si128(x[0], y[0], 0x01);
  t[3] = _mm_clmulepi64_si128(x[0], y[0], 0x10);

  t[2] = _mm_xor_si128(t[2], t[3]);

  t[3] = _mm_srli_si128(t[2], 8);
  t[2] = _mm_slli_si128(t[2], 8);

  t[0] = _mm_xor_si128(t[0], z[0]);
  t[1] = _mm_xor_si128(t[1], z[1]);
  t[2] = _mm_xor_si128(z[2], t[2]);
  t[3] = _mm_xor_si128(z[3], t[3]);

  t[0] = _mm_xor_si128(t[0], t[2]); // t[0] = z[0] + z[2] + t[2]
  t[1] = _mm_xor_si128(t[1], t[3]); // t[1] = z[0] + z[2] + t[3]

  z[1] = _mm_xor_si128(z[1], t[0]);
  z[2] = _mm_xor_si128(z[2], t[1]);

  // modular reduction
  t[0] = _mm_clmulepi64_si128(z[2], irr, 0x01); // 2 ^ 64
  t[1] = _mm_clmulepi64_si128(z[3], irr, 0x00); // 2 ^ 128
  t[2] = _mm_clmulepi64_si128(z[3], irr, 0x01); // 2 ^ 192

  z[0] = _mm_xor_si128(z[0], _mm_slli_si128(t[0], 8));
  z[1] = _mm_xor_si128(z[1], _mm_srli_si128(t[0], 8));
  z[1] = _mm_xor_si128(z[1], t[1]);
  z[1] = _mm_xor_si128(z[1], _mm_slli_si128(t[2], 8));
  z[2] = _mm_xor_si128(z[2], _mm_srli_si128(t[2], 8));

  t[0] = _mm_clmulepi64_si128(z[2], irr, 0x00); // 2 ^ 0
  z[0] = _mm_xor_si128(z[0], t[0]);

  _mm_store_si128((__m128i*)&c[0], z[0]);
  _mm_store_si128((__m128i*)&c[2], z[1]);
}

void GF_mul_N(const GF a[AIMER_NUM_MPC_PARTIES], const GF b,
                    GF c[AIMER_NUM_MPC_PARTIES])
{
  __m128i x[2], y[3], z[4], t[4];
  __m128i irr = _mm_set_epi64x(0x0, 0x425);

  y[0] = _mm_load_si128((const __m128i*)&b[0]); // b0 b1
  y[1] = _mm_load_si128((const __m128i*)&b[2]); // b2 b3
  y[2] = _mm_xor_si128(y[0], y[1]);

  for (size_t party = 0; party < AIMER_NUM_MPC_PARTIES; party++)
  {
    // polynomial multiplication
    x[0] = _mm_load_si128((const __m128i*)&a[party][0]); // a0 a1
    x[1] = _mm_load_si128((const __m128i*)&a[party][2]); // a2 a3

    t[0] = _mm_clmulepi64_si128(x[0], y[0], 0x10);
    t[1] = _mm_clmulepi64_si128(x[0], y[0], 0x01);
    z[0] = _mm_clmulepi64_si128(x[0], y[0], 0x00);
    z[1] = _mm_clmulepi64_si128(x[0], y[0], 0x11);

    t[0] = _mm_xor_si128(t[0], t[1]);
    t[1] = _mm_srli_si128(t[0], 8);
    t[0] = _mm_slli_si128(t[0], 8);
    z[0] = _mm_xor_si128(z[0], t[0]);
    z[1] = _mm_xor_si128(z[1], t[1]);

    t[2] = _mm_clmulepi64_si128(x[1], y[1], 0x10);
    t[3] = _mm_clmulepi64_si128(x[1], y[1], 0x01);
    z[2] = _mm_clmulepi64_si128(x[1], y[1], 0x00);
    z[3] = _mm_clmulepi64_si128(x[1], y[1], 0x11);

    t[2] = _mm_xor_si128(t[2], t[3]);
    t[3] = _mm_srli_si128(t[2], 8);
    t[2] = _mm_slli_si128(t[2], 8);
    z[2] = _mm_xor_si128(z[2], t[2]);
    z[3] = _mm_xor_si128(z[3], t[3]);

    // Start Karat
    x[0] = _mm_xor_si128(x[0], x[1]);

    t[0] = _mm_clmulepi64_si128(x[0], y[2], 0x00);
    t[1] = _mm_clmulepi64_si128(x[0], y[2], 0x11);
    t[2] = _mm_clmulepi64_si128(x[0], y[2], 0x01);
    t[3] = _mm_clmulepi64_si128(x[0], y[2], 0x10);

    t[2] = _mm_xor_si128(t[2], t[3]);

    t[3] = _mm_srli_si128(t[2], 8);
    t[2] = _mm_slli_si128(t[2], 8);

    t[0] = _mm_xor_si128(t[0], z[0]);
    t[1] = _mm_xor_si128(t[1], z[1]);
    t[2] = _mm_xor_si128(z[2], t[2]);
    t[3] = _mm_xor_si128(z[3], t[3]);

    t[0] = _mm_xor_si128(t[0], t[2]); // t[0] = z[0] + z[2] + t[2]
    t[1] = _mm_xor_si128(t[1], t[3]); // t[1] = z[0] + z[2] + t[3]

    z[1] = _mm_xor_si128(z[1], t[0]);
    z[2] = _mm_xor_si128(z[2], t[1]);

    // modular reduction
    t[0] = _mm_clmulepi64_si128(z[2], irr, 0x01); // 2 ^ 64
    t[1] = _mm_clmulepi64_si128(z[3], irr, 0x00); // 2 ^ 128
    t[2] = _mm_clmulepi64_si128(z[3], irr, 0x01); // 2 ^ 192

    z[0] = _mm_xor_si128(z[0], _mm_slli_si128(t[0], 8));
    z[1] = _mm_xor_si128(z[1], _mm_srli_si128(t[0], 8));
    z[1] = _mm_xor_si128(z[1], t[1]);
    z[1] = _mm_xor_si128(z[1], _mm_slli_si128(t[2], 8));
    z[2] = _mm_xor_si128(z[2], _mm_srli_si128(t[2], 8));

    t[0] = _mm_clmulepi64_si128(z[2], irr, 0x00); // 2 ^ 0
    z[0] = _mm_xor_si128(z[0], t[0]);

    _mm_store_si128((__m128i*)&c[party][0], z[0]);
    _mm_store_si128((__m128i*)&c[party][2], z[1]);
  }
}

void GF_mul_add(const GF a, const GF b, GF c)
{
  __m128i x[2], y[2], z[4], t[4];
  __m128i irr = _mm_set_epi64x(0x0, 0x425);

  // polynomial multiplication
  x[0] = _mm_load_si128((const __m128i*)&a[0]); // a0 a1
  x[1] = _mm_load_si128((const __m128i*)&a[2]); // a2 a3
  y[0] = _mm_load_si128((const __m128i*)&b[0]); // b0 b1
  y[1] = _mm_load_si128((const __m128i*)&b[2]); // b2 b3
  z[0] = _mm_load_si128((const __m128i*)&c[0]);
  z[1] = _mm_load_si128((const __m128i*)&c[2]);
  z[2] = _mm_setzero_si128();
  z[3] = _mm_setzero_si128();

  // [t2 t3] = x[0] * y[0]
  t[0] = _mm_clmulepi64_si128(x[0], y[0], 0x10);
  t[1] = _mm_clmulepi64_si128(x[0], y[0], 0x01);
  t[2] = _mm_clmulepi64_si128(x[0], y[0], 0x00);
  t[3] = _mm_clmulepi64_si128(x[0], y[0], 0x11);

  t[0] = _mm_xor_si128(t[0], t[1]);
  t[1] = _mm_srli_si128(t[0], 8);
  t[0] = _mm_slli_si128(t[0], 8);
  t[2] = _mm_xor_si128(t[2], t[0]);
  t[3] = _mm_xor_si128(t[3], t[1]);

  // [z0 z1 z2 z3] += [t2 t3 0 0] + [0 t2 t3 0]
  z[0] = _mm_xor_si128(z[0], t[2]);
  z[1] = _mm_xor_si128(z[1], t[3]);
  z[1] = _mm_xor_si128(z[1], t[2]);
  z[2] = _mm_xor_si128(z[2], t[3]);

  // [t2 t3] = x[1] * y[1]
  t[0] = _mm_clmulepi64_si128(x[1], y[1], 0x10);
  t[1] = _mm_clmulepi64_si128(x[1], y[1], 0x01);
  t[2] = _mm_clmulepi64_si128(x[1], y[1], 0x00);
  t[3] = _mm_clmulepi64_si128(x[1], y[1], 0x11);

  t[0] = _mm_xor_si128(t[0], t[1]);
  t[1] = _mm_srli_si128(t[0], 8);
  t[0] = _mm_slli_si128(t[0], 8);
  t[2] = _mm_xor_si128(t[2], t[0]);
  t[3] = _mm_xor_si128(t[3], t[1]);

  // [z0 z1 z2 z3] += [0 t2 t3 0] + [0 0 t2 t3]
  z[1] = _mm_xor_si128(z[1], t[2]);
  z[2] = _mm_xor_si128(z[2], t[3]);
  z[2] = _mm_xor_si128(z[2], t[2]);
  z[3] = _mm_xor_si128(z[3], t[3]);

  // [t2 t3] = (x[0] + x[1]) * (y[0] + y[1])
  x[0] = _mm_xor_si128(x[0], x[1]);
  y[0] = _mm_xor_si128(y[0], y[1]);

  t[0] = _mm_clmulepi64_si128(x[0], y[0], 0x10);
  t[1] = _mm_clmulepi64_si128(x[0], y[0], 0x01);
  t[2] = _mm_clmulepi64_si128(x[0], y[0], 0x00);
  t[3] = _mm_clmulepi64_si128(x[0], y[0], 0x11);

  t[0] = _mm_xor_si128(t[0], t[1]);
  t[1] = _mm_srli_si128(t[0], 8);
  t[0] = _mm_slli_si128(t[0], 8);
  t[2] = _mm_xor_si128(t[2], t[0]);
  t[3] = _mm_xor_si128(t[3], t[1]);

  // [z0 z1 z2 z3] += [0 t2 t3 0]
  z[1] = _mm_xor_si128(z[1], t[2]);
  z[2] = _mm_xor_si128(z[2], t[3]);

  // modular reduction
  t[0] = _mm_clmulepi64_si128(z[2], irr, 0x01); // 2 ^ 64
  t[1] = _mm_clmulepi64_si128(z[3], irr, 0x00); // 2 ^ 128
  t[2] = _mm_clmulepi64_si128(z[3], irr, 0x01); // 2 ^ 192

  z[0] = _mm_xor_si128(z[0], _mm_slli_si128(t[0], 8));
  z[1] = _mm_xor_si128(z[1], _mm_srli_si128(t[0], 8));
  z[1] = _mm_xor_si128(z[1], t[1]);
  z[1] = _mm_xor_si128(z[1], _mm_slli_si128(t[2], 8));
  z[2] = _mm_xor_si128(z[2], _mm_srli_si128(t[2], 8));

  t[0] = _mm_clmulepi64_si128(z[2], irr, 0x00); // 2 ^ 0
  z[0] = _mm_xor_si128(z[0], t[0]);

  _mm_store_si128((__m128i*)&c[0], z[0]);
  _mm_store_si128((__m128i*)&c[2], z[1]);
}

void GF_mul_add_N(const GF a[AIMER_NUM_MPC_PARTIES], const GF b,
                        GF c[AIMER_NUM_MPC_PARTIES])
{
  __m128i x[2], y[3], z[4], t[4];
  __m128i irr = _mm_set_epi64x(0x0, 0x425);

  y[0] = _mm_load_si128((const __m128i*)&b[0]); // b0 b1
  y[1] = _mm_load_si128((const __m128i*)&b[2]); // b2 b3
  y[2] = _mm_xor_si128(y[0], y[1]);

  for (size_t party = 0; party < AIMER_NUM_MPC_PARTIES; party++)
  {
    // polynomial multiplication
    x[0] = _mm_load_si128((const __m128i*)&a[party][0]); // a0 a1
    x[1] = _mm_load_si128((const __m128i*)&a[party][2]); // a2 a3
    z[0] = _mm_load_si128((const __m128i*)&c[party][0]);
    z[1] = _mm_load_si128((const __m128i*)&c[party][2]);
    z[2] = _mm_setzero_si128();
    z[3] = _mm_setzero_si128();

    // [t2 t3] = x[0] * y[0]
    t[0] = _mm_clmulepi64_si128(x[0], y[0], 0x10);
    t[1] = _mm_clmulepi64_si128(x[0], y[0], 0x01);
    t[2] = _mm_clmulepi64_si128(x[0], y[0], 0x00);
    t[3] = _mm_clmulepi64_si128(x[0], y[0], 0x11);

    t[0] = _mm_xor_si128(t[0], t[1]);
    t[1] = _mm_srli_si128(t[0], 8);
    t[0] = _mm_slli_si128(t[0], 8);
    t[2] = _mm_xor_si128(t[2], t[0]);
    t[3] = _mm_xor_si128(t[3], t[1]);

    // [z0 z1 z2 z3] += [t2 t3 0 0] + [0 t2 t3 0]
    z[0] = _mm_xor_si128(z[0], t[2]);
    z[1] = _mm_xor_si128(z[1], t[3]);
    z[1] = _mm_xor_si128(z[1], t[2]);
    z[2] = _mm_xor_si128(z[2], t[3]);

    // [t2 t3] = x[1] * y[1]
    t[0] = _mm_clmulepi64_si128(x[1], y[1], 0x10);
    t[1] = _mm_clmulepi64_si128(x[1], y[1], 0x01);
    t[2] = _mm_clmulepi64_si128(x[1], y[1], 0x00);
    t[3] = _mm_clmulepi64_si128(x[1], y[1], 0x11);

    t[0] = _mm_xor_si128(t[0], t[1]);
    t[1] = _mm_srli_si128(t[0], 8);
    t[0] = _mm_slli_si128(t[0], 8);
    t[2] = _mm_xor_si128(t[2], t[0]);
    t[3] = _mm_xor_si128(t[3], t[1]);

    // [z0 z1 z2 z3] += [0 t2 t3 0] + [0 0 t2 t3]
    z[1] = _mm_xor_si128(z[1], t[2]);
    z[2] = _mm_xor_si128(z[2], t[3]);
    z[2] = _mm_xor_si128(z[2], t[2]);
    z[3] = _mm_xor_si128(z[3], t[3]);

    // [t2 t3] = (x[0] + x[1]) * (y[0] + y[1])
    x[0] = _mm_xor_si128(x[0], x[1]);

    t[0] = _mm_clmulepi64_si128(x[0], y[2], 0x10);
    t[1] = _mm_clmulepi64_si128(x[0], y[2], 0x01);
    t[2] = _mm_clmulepi64_si128(x[0], y[2], 0x00);
    t[3] = _mm_clmulepi64_si128(x[0], y[2], 0x11);

    t[0] = _mm_xor_si128(t[0], t[1]);
    t[1] = _mm_srli_si128(t[0], 8);
    t[0] = _mm_slli_si128(t[0], 8);
    t[2] = _mm_xor_si128(t[2], t[0]);
    t[3] = _mm_xor_si128(t[3], t[1]);

    // [z0 z1 z2 z3] += [0 t2 t3 0]
    z[1] = _mm_xor_si128(z[1], t[2]);
    z[2] = _mm_xor_si128(z[2], t[3]);

    // modular reduction
    t[0] = _mm_clmulepi64_si128(z[2], irr, 0x01); // 2 ^ 64
    t[1] = _mm_clmulepi64_si128(z[3], irr, 0x00); // 2 ^ 128
    t[2] = _mm_clmulepi64_si128(z[3], irr, 0x01); // 2 ^ 192

    z[0] = _mm_xor_si128(z[0], _mm_slli_si128(t[0], 8));
    z[1] = _mm_xor_si128(z[1], _mm_srli_si128(t[0], 8));
    z[1] = _mm_xor_si128(z[1], t[1]);
    z[1] = _mm_xor_si128(z[1], _mm_slli_si128(t[2], 8));
    z[2] = _mm_xor_si128(z[2], _mm_srli_si128(t[2], 8));

    t[0] = _mm_clmulepi64_si128(z[2], irr, 0x00); // 2 ^ 0
    z[0] = _mm_xor_si128(z[0], t[0]);

    _mm_store_si128((__m128i*)&c[party][0], z[0]);
    _mm_store_si128((__m128i*)&c[party][2], z[1]);
  }
}

void GF_sqr(const GF a, GF c)
{
  __m128i x[3], z[4];
  __m128i irr = _mm_set_epi64x(0x0, 0x425);

  // polynomial multiplication
  x[0] = _mm_load_si128((const __m128i*)&a[0]); // a0 a1
  x[1] = _mm_load_si128((const __m128i*)&a[2]); // a2 a3

  z[0] = _mm_clmulepi64_si128(x[0], x[0], 0x00);
  z[1] = _mm_clmulepi64_si128(x[0], x[0], 0x11);
  z[2] = _mm_clmulepi64_si128(x[1], x[1], 0x00);
  z[3] = _mm_clmulepi64_si128(x[1], x[1], 0x11);

  // modular reduction
  x[0] = _mm_clmulepi64_si128(z[2], irr, 0x01); // 2 ^ 64
  x[1] = _mm_clmulepi64_si128(z[3], irr, 0x00); // 2 ^ 128
  x[2] = _mm_clmulepi64_si128(z[3], irr, 0x01); // 2 ^ 192

  z[0] = _mm_xor_si128(z[0], _mm_slli_si128(x[0], 8));
  z[1] = _mm_xor_si128(z[1], _mm_srli_si128(x[0], 8));
  z[1] = _mm_xor_si128(z[1], x[1]);
  z[1] = _mm_xor_si128(z[1], _mm_slli_si128(x[2], 8));
  z[2] = _mm_xor_si128(z[2], _mm_srli_si128(x[2], 8));

  x[0] = _mm_clmulepi64_si128(z[2], irr, 0x00); // 2 ^ 0
  z[0] = _mm_xor_si128(z[0], x[0]);

  _mm_store_si128((__m128i*)&c[0], z[0]);
  _mm_store_si128((__m128i*)&c[2], z[1]);
}

void GF_sqr_N(const GF a[AIMER_NUM_MPC_PARTIES],
                    GF c[AIMER_NUM_MPC_PARTIES])
{
  __m128i x[6], z[8];
  __m128i irr = _mm_set_epi64x(0x0, 0x425);
  for (size_t party = 0; party < AIMER_NUM_MPC_PARTIES; party += 2)
  {
    // polynomial multiplication
    x[0] = _mm_load_si128((const __m128i*)&a[party][0]);     // a0 a1
    x[1] = _mm_load_si128((const __m128i*)&a[party][2]);     // a2 a3
    x[2] = _mm_load_si128((const __m128i*)&a[party + 1][0]); // a0 a1
    x[3] = _mm_load_si128((const __m128i*)&a[party + 1][2]); // a2 a3

    z[0] = _mm_clmulepi64_si128(x[0], x[0], 0x00);
    z[1] = _mm_clmulepi64_si128(x[0], x[0], 0x11);
    z[2] = _mm_clmulepi64_si128(x[1], x[1], 0x00);
    z[3] = _mm_clmulepi64_si128(x[1], x[1], 0x11);

    z[4] = _mm_clmulepi64_si128(x[2], x[2], 0x00);
    z[5] = _mm_clmulepi64_si128(x[2], x[2], 0x11);
    z[6] = _mm_clmulepi64_si128(x[3], x[3], 0x00);
    z[7] = _mm_clmulepi64_si128(x[3], x[3], 0x11);

    // modular reduction
    x[0] = _mm_clmulepi64_si128(z[2], irr, 0x01); // 2 ^ 64
    x[1] = _mm_clmulepi64_si128(z[3], irr, 0x00); // 2 ^ 128
    x[2] = _mm_clmulepi64_si128(z[3], irr, 0x01); // 2 ^ 192

    x[3] = _mm_clmulepi64_si128(z[6], irr, 0x01); // 2 ^ 64
    x[4] = _mm_clmulepi64_si128(z[7], irr, 0x00); // 2 ^ 128
    x[5] = _mm_clmulepi64_si128(z[7], irr, 0x01); // 2 ^ 192

    z[0] = _mm_xor_si128(z[0], _mm_slli_si128(x[0], 8));
    z[1] = _mm_xor_si128(z[1], _mm_srli_si128(x[0], 8));
    z[1] = _mm_xor_si128(z[1], x[1]);
    z[1] = _mm_xor_si128(z[1], _mm_slli_si128(x[2], 8));
    z[2] = _mm_xor_si128(z[2], _mm_srli_si128(x[2], 8));

    z[4] = _mm_xor_si128(z[4], _mm_slli_si128(x[3], 8));
    z[5] = _mm_xor_si128(z[5], _mm_srli_si128(x[3], 8));
    z[5] = _mm_xor_si128(z[5], x[4]);
    z[5] = _mm_xor_si128(z[5], _mm_slli_si128(x[5], 8));
    z[6] = _mm_xor_si128(z[6], _mm_srli_si128(x[5], 8));

    x[0] = _mm_clmulepi64_si128(z[2], irr, 0x00); // 2 ^ 0
    x[1] = _mm_clmulepi64_si128(z[6], irr, 0x00); // 2 ^ 0
    z[0] = _mm_xor_si128(z[0], x[0]);
    z[4] = _mm_xor_si128(z[4], x[1]);

    _mm_store_si128((__m128i*)&c[party][0], z[0]);
    _mm_store_si128((__m128i*)&c[party][2], z[1]);
    _mm_store_si128((__m128i*)&c[party + 1][0], z[4]);
    _mm_store_si128((__m128i*)&c[party + 1][2], z[5]);
  }
}

void GF_transposed_matmul(const GF a, const GF b[AIM2_NUM_BITS_FIELD], GF c)
{
  const __m256i shift = _mm256_set_epi64x(0, 1, 2, 3);
  const __m256i zero = _mm256_setzero_si256();

  __m256i c0 = _mm256_setzero_si256();
  __m256i c1 = _mm256_setzero_si256();
  __m256i c2 = _mm256_setzero_si256();
  __m256i c3 = _mm256_setzero_si256();
  __m256i m0, m1, m2, m3, a0, a1, a2, a3;
  __m256i mask;

  for (int i = 3; i >= 0; i--)
  {
    mask = _mm256_set1_epi64x(a[i]);
    mask = _mm256_sllv_epi64(mask, shift);
    for (int row = 64 * (i + 1); row > 64 * i; row -= 4)
    {
      m0 = _mm256_load_si256((const __m256i*)b[row - 4]);
      m1 = _mm256_load_si256((const __m256i*)b[row - 3]);
      m2 = _mm256_load_si256((const __m256i*)b[row - 2]);
      m3 = _mm256_load_si256((const __m256i*)b[row - 1]);

      a0 = _mm256_permute4x64_epi64(mask, 0b00000000);
      a1 = _mm256_permute4x64_epi64(mask, 0b01010101);
      a2 = _mm256_permute4x64_epi64(mask, 0b10101010);
      a3 = _mm256_permute4x64_epi64(mask, 0b11111111);

      a0 = _mm256_cmpgt_epi64(zero, a0);
      a1 = _mm256_cmpgt_epi64(zero, a1);
      a2 = _mm256_cmpgt_epi64(zero, a2);
      a3 = _mm256_cmpgt_epi64(zero, a3);

      c0 = _mm256_xor_si256(c0, _mm256_and_si256(m0, a0));
      c1 = _mm256_xor_si256(c1, _mm256_and_si256(m1, a1));
      c2 = _mm256_xor_si256(c2, _mm256_and_si256(m2, a2));
      c3 = _mm256_xor_si256(c3, _mm256_and_si256(m3, a3));

      mask = _mm256_slli_epi64(mask, 4);
    }
  }
  c0 = _mm256_xor_si256(c0, c1);
  c2 = _mm256_xor_si256(c2, c3);
  c0 = _mm256_xor_si256(c0, c2);
  _mm256_store_si256((__m256i*)c, c0);
}

void GF_transposed_matmul_add_N(const GF a[AIMER_NUM_MPC_PARTIES],
                                const GF b[AIM2_NUM_BITS_FIELD],
                                GF c[AIMER_NUM_MPC_PARTIES])
{
  const __m256i shift = _mm256_set_epi64x(0, 1, 2, 3);
  const __m256i zero = _mm256_setzero_si256();

  __m256i m0, m1, m2, m3, a0, a1, a2, a3, c0, c1, c2, c3;
  __m256i mask1, mask2;

  for (size_t party = 0; party < AIMER_NUM_MPC_PARTIES; party += 2)
  {
    c0 = _mm256_load_si256((const __m256i*)c[party]);
    c1 = _mm256_setzero_si256();
    c2 = _mm256_load_si256((const __m256i*)c[party + 1]);
    c3 = _mm256_setzero_si256();

    for (int i = 3; i >= 0; i--)
    {
      mask1 = _mm256_set1_epi64x(a[party][i]);
      mask2 = _mm256_set1_epi64x(a[party + 1][i]);
      mask1 = _mm256_sllv_epi64(mask1, shift);
      mask2 = _mm256_sllv_epi64(mask2, shift);
      for (int row = 64 * (i + 1); row > 64 * i; row -= 4)
      {
        m0 = _mm256_load_si256((const __m256i*)b[row - 4]);
        m1 = _mm256_load_si256((const __m256i*)b[row - 3]);
        m2 = _mm256_load_si256((const __m256i*)b[row - 2]);
        m3 = _mm256_load_si256((const __m256i*)b[row - 1]);

        a0 = _mm256_permute4x64_epi64(mask1, 0b00000000);
        a1 = _mm256_permute4x64_epi64(mask1, 0b01010101);
        a2 = _mm256_permute4x64_epi64(mask1, 0b10101010);
        a3 = _mm256_permute4x64_epi64(mask1, 0b11111111);

        a0 = _mm256_cmpgt_epi64(zero, a0);
        a1 = _mm256_cmpgt_epi64(zero, a1);
        a2 = _mm256_cmpgt_epi64(zero, a2);
        a3 = _mm256_cmpgt_epi64(zero, a3);

        a0 = _mm256_and_si256(m0, a0);
        a1 = _mm256_and_si256(m1, a1);
        a2 = _mm256_and_si256(m2, a2);
        a3 = _mm256_and_si256(m3, a3);

        c0 = _mm256_xor_si256(c0, a0);
        c1 = _mm256_xor_si256(c1, a1);
        c0 = _mm256_xor_si256(c0, a2);
        c1 = _mm256_xor_si256(c1, a3);

        a0 = _mm256_permute4x64_epi64(mask2, 0b00000000);
        a1 = _mm256_permute4x64_epi64(mask2, 0b01010101);
        a2 = _mm256_permute4x64_epi64(mask2, 0b10101010);
        a3 = _mm256_permute4x64_epi64(mask2, 0b11111111);

        a0 = _mm256_cmpgt_epi64(zero, a0);
        a1 = _mm256_cmpgt_epi64(zero, a1);
        a2 = _mm256_cmpgt_epi64(zero, a2);
        a3 = _mm256_cmpgt_epi64(zero, a3);

        a0 = _mm256_and_si256(m0, a0);
        a1 = _mm256_and_si256(m1, a1);
        a2 = _mm256_and_si256(m2, a2);
        a3 = _mm256_and_si256(m3, a3);

        c2 = _mm256_xor_si256(c2, a0);
        c3 = _mm256_xor_si256(c3, a1);
        c2 = _mm256_xor_si256(c2, a2);
        c3 = _mm256_xor_si256(c3, a3);

        mask1 = _mm256_slli_epi64(mask1, 4);
        mask2 = _mm256_slli_epi64(mask2, 4);
      }
    }
    c0 = _mm256_xor_si256(c0, c1);
    c2 = _mm256_xor_si256(c2, c3);
    _mm256_store_si256((__m256i*)c[party], c0);
    _mm256_store_si256((__m256i*)c[party + 1], c2);
  }
}

void POLY_mul_add_N(const GF a[AIMER_NUM_MPC_PARTIES], const GF b,
                    GF lo[AIMER_NUM_MPC_PARTIES],
                    GF hi[AIMER_NUM_MPC_PARTIES])
{
  __m128i x[2], y[3], z[4], t[4];

  y[0] = _mm_load_si128((const __m128i*)&b[0]); // b0 b1
  y[1] = _mm_load_si128((const __m128i*)&b[2]); // b2 b3
  y[2] = _mm_xor_si128(y[0], y[1]);

  for (size_t party = 0; party < AIMER_NUM_MPC_PARTIES; party++)
  {
    // polynomial multiplication
    x[0] = _mm_load_si128((const __m128i*)&a[party][0]); // a0 a1
    x[1] = _mm_load_si128((const __m128i*)&a[party][2]); // a2 a3
    z[0] = _mm_load_si128((const __m128i*)&lo[party][0]);
    z[1] = _mm_load_si128((const __m128i*)&lo[party][2]);
    z[2] = _mm_load_si128((const __m128i*)&hi[party][0]);
    z[3] = _mm_load_si128((const __m128i*)&hi[party][2]);

    // [t2 t3] = x[0] * y[0]
    t[0] = _mm_clmulepi64_si128(x[0], y[0], 0x10);
    t[1] = _mm_clmulepi64_si128(x[0], y[0], 0x01);
    t[2] = _mm_clmulepi64_si128(x[0], y[0], 0x00);
    t[3] = _mm_clmulepi64_si128(x[0], y[0], 0x11);

    t[0] = _mm_xor_si128(t[0], t[1]);
    t[1] = _mm_srli_si128(t[0], 8);
    t[0] = _mm_slli_si128(t[0], 8);
    t[2] = _mm_xor_si128(t[2], t[0]);
    t[3] = _mm_xor_si128(t[3], t[1]);

    // [z0 z1 z2 z3] += [t2 t3 0 0] + [0 t2 t3 0]
    z[0] = _mm_xor_si128(z[0], t[2]);
    z[1] = _mm_xor_si128(z[1], t[3]);
    z[1] = _mm_xor_si128(z[1], t[2]);
    z[2] = _mm_xor_si128(z[2], t[3]);

    // [t2 t3] = x[1] * y[1]
    t[0] = _mm_clmulepi64_si128(x[1], y[1], 0x10);
    t[1] = _mm_clmulepi64_si128(x[1], y[1], 0x01);
    t[2] = _mm_clmulepi64_si128(x[1], y[1], 0x00);
    t[3] = _mm_clmulepi64_si128(x[1], y[1], 0x11);

    t[0] = _mm_xor_si128(t[0], t[1]);
    t[1] = _mm_srli_si128(t[0], 8);
    t[0] = _mm_slli_si128(t[0], 8);
    t[2] = _mm_xor_si128(t[2], t[0]);
    t[3] = _mm_xor_si128(t[3], t[1]);

    // [z0 z1 z2 z3] += [0 t2 t3 0] + [0 0 t2 t3]
    z[1] = _mm_xor_si128(z[1], t[2]);
    z[2] = _mm_xor_si128(z[2], t[3]);
    z[2] = _mm_xor_si128(z[2], t[2]);
    z[3] = _mm_xor_si128(z[3], t[3]);

    // [t2 t3] = (x[0] + x[1]) * (y[0] + y[1])
    x[0] = _mm_xor_si128(x[0], x[1]);

    t[0] = _mm_clmulepi64_si128(x[0], y[2], 0x10);
    t[1] = _mm_clmulepi64_si128(x[0], y[2], 0x01);
    t[2] = _mm_clmulepi64_si128(x[0], y[2], 0x00);
    t[3] = _mm_clmulepi64_si128(x[0], y[2], 0x11);

    t[0] = _mm_xor_si128(t[0], t[1]);
    t[1] = _mm_srli_si128(t[0], 8);
    t[0] = _mm_slli_si128(t[0], 8);
    t[2] = _mm_xor_si128(t[2], t[0]);
    t[3] = _mm_xor_si128(t[3], t[1]);

    // [z0 z1 z2 z3] += [0 t2 t3 0]
    z[1] = _mm_xor_si128(z[1], t[2]);
    z[2] = _mm_xor_si128(z[2], t[3]);

    _mm_store_si128((__m128i*)&lo[party][0], z[0]);
    _mm_store_si128((__m128i*)&lo[party][2], z[1]);
    _mm_store_si128((__m128i*)&hi[party][0], z[2]);
    _mm_store_si128((__m128i*)&hi[party][2], z[3]);
  }
}

void POLY_red_N(const GF hi[AIMER_NUM_MPC_PARTIES],
                      GF lo[AIMER_NUM_MPC_PARTIES])
{
  __m128i x[6], z[8];
  __m128i irr = _mm_set_epi64x(0x0, 0x425);
  for (size_t party = 0; party < AIMER_NUM_MPC_PARTIES; party += 2)
  {
    // load
    z[0] = _mm_load_si128((const __m128i*)&lo[party][0]);
    z[1] = _mm_load_si128((const __m128i*)&lo[party][2]);
    z[2] = _mm_load_si128((const __m128i*)&hi[party][0]);
    z[3] = _mm_load_si128((const __m128i*)&hi[party][2]);

    z[4] = _mm_load_si128((const __m128i*)&lo[party + 1][0]);
    z[5] = _mm_load_si128((const __m128i*)&lo[party + 1][2]);
    z[6] = _mm_load_si128((const __m128i*)&hi[party + 1][0]);
    z[7] = _mm_load_si128((const __m128i*)&hi[party + 1][2]);

    // modular reduction
    x[0] = _mm_clmulepi64_si128(z[2], irr, 0x01); // 2 ^ 64
    x[1] = _mm_clmulepi64_si128(z[3], irr, 0x00); // 2 ^ 128
    x[2] = _mm_clmulepi64_si128(z[3], irr, 0x01); // 2 ^ 192

    x[3] = _mm_clmulepi64_si128(z[6], irr, 0x01); // 2 ^ 64
    x[4] = _mm_clmulepi64_si128(z[7], irr, 0x00); // 2 ^ 128
    x[5] = _mm_clmulepi64_si128(z[7], irr, 0x01); // 2 ^ 192

    z[0] = _mm_xor_si128(z[0], _mm_slli_si128(x[0], 8));
    z[1] = _mm_xor_si128(z[1], _mm_srli_si128(x[0], 8));
    z[1] = _mm_xor_si128(z[1], x[1]);
    z[1] = _mm_xor_si128(z[1], _mm_slli_si128(x[2], 8));
    z[2] = _mm_xor_si128(z[2], _mm_srli_si128(x[2], 8));

    z[4] = _mm_xor_si128(z[4], _mm_slli_si128(x[3], 8));
    z[5] = _mm_xor_si128(z[5], _mm_srli_si128(x[3], 8));
    z[5] = _mm_xor_si128(z[5], x[4]);
    z[5] = _mm_xor_si128(z[5], _mm_slli_si128(x[5], 8));
    z[6] = _mm_xor_si128(z[6], _mm_srli_si128(x[5], 8));

    x[0] = _mm_clmulepi64_si128(z[2], irr, 0x00); // 2 ^ 0
    x[1] = _mm_clmulepi64_si128(z[6], irr, 0x00); // 2 ^ 0
    z[0] = _mm_xor_si128(z[0], x[0]);
    z[4] = _mm_xor_si128(z[4], x[1]);

    _mm_store_si128((__m128i*)&lo[party][0], z[0]);
    _mm_store_si128((__m128i*)&lo[party][2], z[1]);
    _mm_store_si128((__m128i*)&lo[party + 1][0], z[4]);
    _mm_store_si128((__m128i*)&lo[party + 1][2], z[5]);
  }
}
