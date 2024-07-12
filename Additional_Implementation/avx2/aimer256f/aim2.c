// SPDX-License-Identifier: MIT

#include <stdint.h>
#include <stdalign.h>
#include "aim2.h"
#include "field.h"
#include "hash.h"
#include "params.h"

// inverse Mersenne S-box with e1 = 11
// (2 ^ 11 - 1) ^ (-1) mod (2 ^ 256 - 1)
// = 0xb6d6dadb5b6b6d6dadb5b6b6d6dadb5b6b6d6dadb5b6b6d6dadb5b6b6d6dadb5
// b6d6dadb5 b6 b6d6dadb5 b6 b6d6dadb5 b6 b6d6dadb5 b6 b6d6dadb5 b6 b6d6dadb5
void GF_exp_invmer_e_1(const GF in, GF out)
{
  unsigned int i;
  alignas(32) GF t1 = {0,};
  alignas(32) GF table_5 = {0,}, table_6 = {0,};
  alignas(32) GF table_a = {0,}, table_b = {0,}, table_d = {0,};

  // t1 = in ^ 4
  GF_sqr(in, table_d);
  GF_sqr(table_d, t1);

  // table_5 = in ^ 5
  GF_mul(t1, in, table_5);
  // table_6 = in ^ 6
  GF_mul(table_5, in, table_6);
  // table_a = in ^ 10 = (in ^ 5) ^ 2
  GF_sqr(table_5, table_a);
  // table_b = in ^ 11
  GF_mul(table_5, table_6, table_b);
  // table_d = in ^ 13
  GF_mul(table_b, table_d, table_d);

  // table_b = in ^ (0xb6), table_5 = in ^ (0xb5)
  GF_sqr(table_b, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, table_5, table_5);
  GF_mul(t1, table_6, table_b);

  // t1 = in ^ (0xb6 d)
  GF_sqr(table_b, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, table_d, t1);

  // t1 = in ^ (0xb6d 6)
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, table_6, t1);

  // t1 = in ^ (0xb6d6 d)
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, table_d, t1);

  // t1 = in ^ (0xb6d6d a)
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, table_a, t1);

  // t1 = in ^ (0xb6d6da d)
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, table_d, t1);

  // table_5 = in ^ (0xb6d6dad b5)
  for (i = 0; i < 8; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_5, table_5);

  // t1 = in ^ (0xb6d6dadb5 b6)
  GF_sqr(table_5, t1);
  for (i = 1; i < 8; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_b, t1);

  // t1 = in ^ (0xb6d6dadb5b6 b6d6dadb5)
  for (i = 0; i < 36; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_5, t1);

  // t1 = in ^ (0xb6d6dadb5b6b6d6dadb5 b6)
  for (i = 0; i < 8; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_b, t1);

  // t1 = in ^ (0xb6d6dadb5b6b6d6dadb5b6 b6d6dadb5)
  for (i = 0; i < 36; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_5, t1);

  // t1 = in ^ (0xb6d6dadb5b6b6d6dadb5b6b6d6dadb5 b6)
  for (i = 0; i < 8; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_b, t1);

  // t1 = in ^ (0xb6d6dadb5b6b6d6dadb5b6b6d6dadb5b6 b6d6dadb5)
  for (i = 0; i < 36; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_5, t1);

  // t1 = in ^ (0xb6d6dadb5b6b6d6dadb5b6b6d6dadb5b6b6d6dadb5 b6)
  for (i = 0; i < 8; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_b, t1);

  // t1 = in ^ (0xb6d6dadb5b6b6d6dadb5b6b6d6dadb5b6b6d6dadb5b6 b6d6dadb5)
  for (i = 0; i < 36; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_5, t1);

  // t1 = in ^ (0xb6d6dadb5b6b6d6dadb5b6b6d6dadb5b6b6d6dadb5b6b6d6dadb5 b6)
  for (i = 0; i < 8; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_b, t1);

  // out = in ^ (0xb6d6dadb5b6b6d6dadb5b6b6d6dadb5b6b6d6dadb5b6b6d6dadb5b6 b6d6dadb5)
  for (i = 0; i < 36; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_5, out);
}

// inverse Mersenne S-box with e2 = 141
// (2 ^ 141 - 1) ^ (-1) mod (2 ^ 256 - 1)
// = 0x2224448889112222444888911222244488911122244448891112224444889111
// 222444 8889112 222444 8889112 222444 889111 222444 4889111 222444 4889111
void GF_exp_invmer_e_2(const GF in, GF out)
{
  unsigned int i;
  alignas(32) GF t1 = {0,}, t2 = {0,}, t3 = {0,}, t4 = {0,}, t5 = {0,};
  alignas(32) GF table_9 = {0,};

  // t2 = in ^ (0x11), table_9 = in ^ 9
  GF_sqr(in, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, in, table_9);
  GF_sqr(t1, t1);
  GF_mul(t1, in, t2);

  // t3 = in ^ (0x111)
  GF_sqr(t2, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, in, t3);

  // t4 = in ^ (0x222444)
  GF_sqr(t3, t1);
  for (i = 0; i < 10; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, t3, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t4);

  // t1 = in ^ (0x222444 8889)
  GF_sqr(t4, t1);
  for (i = 1; i < 9; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, t3, t1);

  for (i = 0; i < 7; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_9, t1);

  // t1 = in ^ (0x2224448889 11)
  for (i = 0; i < 8; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, t2, t1);

  // t5 = in ^ (0x222444888911 2)
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, in, t1);
  GF_sqr(t1, t5);

  // t1 = in ^ (0x2224448889112 2224448889112)
  GF_sqr(t5, t1);
  for (i = 1; i < 52; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, t5, t1);

  // t1 = in ^ (0x22244488891122224448889112 222444)
  for (i = 0; i < 24; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, t4, t1);

  // t1 = in ^ (0x22244488891122224448889112222444 889)
  for (i = 0; i < 5; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, t2, t1);

  for (i = 0; i < 7; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_9, t1);

  // t1 = in ^ (0x22244488891122224448889112222444889 111)
  for (i = 0; i < 12; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, t3, t1);

  // t1 = in ^ (0x22244488891122224448889112222444889111 222444)
  for (i = 0; i < 24; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, t4, t1);

  // t1 = in ^ (0x22244488891122224448889112222444889111222444 4)
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, in, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);

  // t1 = in ^ (0x222444888911222244488891122224448891112224444 889)
  for (i = 0; i < 5; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, t2, t1);

  for (i = 0; i < 7; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_9, t1);

  // t1 = in ^ (0x222444888911222244488891122224448891112224444889 111)
  for (i = 0; i < 12; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, t3, t1);

  // t1 = in ^ (0x222444888911222244488891122224448891112224444889111 222444)
  for (i = 0; i < 24; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, t4, t1);

  // t1 = in ^ (0x222444888911222244488891122224448891112224444889111222444 4)
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, in, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);

  // t1 = in ^ (0x2224448889112222444888911222244488911122244448891112224444 889)
  for (i = 0; i < 5; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, t2, t1);

  for (i = 0; i < 7; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_9, t1);

  // out = in ^ (0x2224448889112222444888911222244488911122244448891112224444889 111)
  for (i = 0; i < 12; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, t3, out);
}

// inverse Mersenne S-box with e3 = 7
// (2 ^ 7 - 1) ^ (-1) mod (2 ^ 256 - 1)
// = 0xddbb76eddbb76eddbb76eddbb76eddbb76eddbb76eddbb76eddbb76eddbb76ed
// ddbb76e ddbb76e ddbb76e ddbb76e ddbb76e ddbb76e ddbb76e ddbb76e ddbb76e d
void GF_exp_invmer_e_3(const GF in, GF out)
{
  unsigned int i;
  alignas(32) GF t1 = {0,};
  alignas(32) GF table_6 = {0,}, table_7 = {0,}, table_b = {0,}, table_d = {0,};

  // t1 = in ^ 3
  GF_sqr(in, table_d);
  GF_mul(table_d, in, t1);

  // table_6 = in ^ 6
  GF_sqr(t1, table_6);
  // table_7 = in ^ 7
  GF_mul(table_6, in, table_7);
  // table_b = in ^ 11
  GF_sqr(table_d, table_b);
  GF_mul(table_7, table_b, table_b);
  // table_d = in ^ 13
  GF_mul(table_b, table_d, table_d);

  // t1 = in ^ 0xdd
  GF_sqr(table_d, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, table_d, t1);

  // t1 = in ^ 0xdd b
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, table_b, t1);

  // t1 = in ^ 0xddb b
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, table_b, t1);

  // t1 = in ^ 0xddbb 7
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, table_7, t1);

  // t1 = in ^ 0xddbb7 6
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, table_6, t1);

  // table_7 = in ^ 0xddbb76 e
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, table_7, t1);
  GF_sqr(t1, table_7);

  // t1 = in ^ 0xddbb76e ddbb76e
  GF_sqr(table_7, t1);
  for (i = 1; i < 28; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_7, t1);

  // t1 = in ^ 0xddbb76eddbb76e ddbb76e
  for (i = 0; i < 28; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_7, t1);

  // t1 = in ^ 0xddbb76eddbb76eddbb76e ddbb76e
  for (i = 0; i < 28; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_7, t1);

  // t1 = in ^ 0xddbb76eddbb76eddbb76eddbb76e ddbb76e
  for (i = 0; i < 28; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_7, t1);

  // t1 = in ^ 0xddbb76eddbb76eddbb76eddbb76eddbb76e ddbb76e
  for (i = 0; i < 28; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_7, t1);

  // t1 = in ^ 0xddbb76eddbb76eddbb76eddbb76eddbb76eddbb76e ddbb76e
  for (i = 0; i < 28; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_7, t1);

  // t1 = in ^ 0xddbb76eddbb76eddbb76eddbb76eddbb76eddbb76eddbb76e ddbb76e
  for (i = 0; i < 28; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_7, t1);

  // t1 = in ^ 0xddbb76eddbb76eddbb76eddbb76eddbb76eddbb76eddbb76eddbb76e ddbb76e
  for (i = 0; i < 28; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_7, t1);

  // out = in ^ 0xddbb76eddbb76eddbb76eddbb76eddbb76eddbb76eddbb76eddbb76eddbb76e d
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, table_d, out);
}

// Mersenne exponentiation with e_star = 3
void GF_exp_mer_e_star(const GF in, GF out)
{
  alignas(32) GF t1 = {0,};

  // t1 = a ^ (2 ^ 2 - 1)
  GF_sqr(in, t1);
  GF_mul(t1, in, t1);

  // out = a ^ (2 ^ 3 - 1)
  GF_sqr(t1, t1);
  GF_mul(t1, in, out);
}

void generate_matrices_L_and_U(const uint8_t iv[AIM2_IV_SIZE],
                               GF matrix_L[AIM2_NUM_INPUT_SBOX][AIM2_NUM_BITS_FIELD],
                               GF matrix_U[AIM2_NUM_INPUT_SBOX][AIM2_NUM_BITS_FIELD],
                               GF vector_b)
{
  uint8_t buf[AIM2_NUM_BYTES_FIELD];
  uint64_t ormask, lmask, umask;
  hash_instance ctx;
  GF g = {0,};

  // initialize hash
  hash_init(&ctx);
  hash_update(&ctx, iv, AIM2_IV_SIZE);
  hash_final(&ctx);

  for (size_t num = 0; num < AIM2_NUM_INPUT_SBOX; num++)
  {
    for (size_t row = 0; row < AIM2_NUM_BITS_FIELD; row++)
    {
      hash_squeeze(&ctx, buf, AIM2_NUM_BYTES_FIELD);
      GF_from_bytes(buf, g);

      ormask = ((uint64_t)1) << (row % 64);
      lmask = ((uint64_t)-1) << (row % 64);
      umask = ~lmask;

      size_t inter = row / 64;
      size_t col_word;
      for (col_word = 0; col_word < inter; col_word++)
      {
        // L is zero, U is full
        matrix_L[num][row][col_word] = 0;
        matrix_U[num][row][col_word] = g[col_word];
      }
      matrix_L[num][row][inter] = (g[inter] & lmask) | ormask;
      matrix_U[num][row][inter] = (g[inter] & umask) | ormask;
      for (col_word = inter + 1; col_word < AIM2_NUM_WORDS_FIELD; col_word++)
      {
        // L is full, U is zero
        matrix_L[num][row][col_word] = g[col_word];
        matrix_U[num][row][col_word] = 0;
      }
    }
  }

  hash_squeeze(&ctx, (uint8_t*)vector_b, AIM2_NUM_BYTES_FIELD);
}

void generate_matrix_LU(const uint8_t iv[AIM2_IV_SIZE],
                        GF matrix_A[AIM2_NUM_INPUT_SBOX][AIM2_NUM_BITS_FIELD],
                        GF vector_b)
{
  alignas(32) GF matrix_L[AIM2_NUM_INPUT_SBOX][AIM2_NUM_BITS_FIELD];
  alignas(32) GF matrix_U[AIM2_NUM_INPUT_SBOX][AIM2_NUM_BITS_FIELD];

  generate_matrices_L_and_U(iv, matrix_L, matrix_U, vector_b);

  for (size_t num = 0; num < AIM2_NUM_INPUT_SBOX; num++)
  {
    for (size_t i = 0; i < AIM2_NUM_BITS_FIELD; i++)
    {
      GF_transposed_matmul(matrix_U[num][i], matrix_L[num], matrix_A[num][i]);
    }
  }
}

void aim2(const uint8_t pt[AIM2_NUM_BYTES_FIELD],
          const uint8_t iv[AIM2_IV_SIZE],
          uint8_t ct[AIM2_NUM_BYTES_FIELD])
{
  alignas(32) GF matrix_L[AIM2_NUM_INPUT_SBOX][AIM2_NUM_BITS_FIELD];
  alignas(32) GF matrix_U[AIM2_NUM_INPUT_SBOX][AIM2_NUM_BITS_FIELD];
  alignas(32) GF vector_b = {0,};

  alignas(32) GF state[AIM2_NUM_INPUT_SBOX];
  alignas(32) GF pt_GF = {0,}, ct_GF = {0,};
  GF_from_bytes(pt, pt_GF);

  // generate random matrix
  generate_matrices_L_and_U(iv, matrix_L, matrix_U, vector_b);

  // linear component: constant addition
  GF_add(pt_GF, aim2_constants[0], state[0]);
  GF_add(pt_GF, aim2_constants[1], state[1]);
  GF_add(pt_GF, aim2_constants[2], state[2]);

  // non-linear component: inverse Mersenne S-box
  GF_exp_invmer_e_1(state[0], state[0]);
  GF_exp_invmer_e_2(state[1], state[1]);
  GF_exp_invmer_e_3(state[2], state[2]);

  // linear component: affine layer
  GF_transposed_matmul(state[0], matrix_U[0], state[0]);
  GF_transposed_matmul(state[0], matrix_L[0], state[0]);

  GF_transposed_matmul(state[1], matrix_U[1], state[1]);
  GF_transposed_matmul(state[1], matrix_L[1], state[1]);

  GF_transposed_matmul(state[2], matrix_U[2], state[2]);
  GF_transposed_matmul(state[2], matrix_L[2], state[2]);

  GF_add(state[0], state[1], state[0]);
  GF_add(state[2], vector_b, state[2]);
  GF_add(state[0], state[2], state[0]);

  // non-linear component: Mersenne S-box
  GF_exp_mer_e_star(state[0], state[0]);

  // linear component: feed-forward
  GF_add(state[0], pt_GF, ct_GF);

  GF_to_bytes(ct_GF, ct);
}

void aim2_sbox_outputs(const GF pt, GF sbox_outputs[AIM2_NUM_INPUT_SBOX])
{
  // linear component: constant addition
  GF_add(pt, aim2_constants[0], sbox_outputs[0]);
  GF_add(pt, aim2_constants[1], sbox_outputs[1]);
  GF_add(pt, aim2_constants[2], sbox_outputs[2]);

  // non-linear component: inverse Mersenne S-box
  GF_exp_invmer_e_1(sbox_outputs[0], sbox_outputs[0]);
  GF_exp_invmer_e_2(sbox_outputs[1], sbox_outputs[1]);
  GF_exp_invmer_e_3(sbox_outputs[2], sbox_outputs[2]);
}
