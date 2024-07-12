// SPDX-License-Identifier: MIT

#include <stdint.h>
#include "aim2.h"
#include "field.h"
#include "hash.h"
#include "params.h"

// inverse Mersenne S-box with e1 = 17
// (2 ^ 17 - 1) ^ (-1) mod (2 ^ 192 - 1)
// = 0xad6b56b5ab5ad5ad6ad6b56b5ab5ad5ad6ad6b56b5ab5ad5
// ad6b56b5ab5ad5 ad6 ad6b56b5ab5ad5 ad6 ad6b56b5ab5ad5
void GF_exp_invmer_e_1(const GF in, GF out)
{
  unsigned int i;
  GF t1 = {0,}, t2 = {0,};
  GF table_5 = {0,}, table_6 = {0,};
  GF table_a = {0,}, table_b = {0,}, table_d = {0,};

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
  GF_mul(table_a, in, table_b);
  // table_d = in ^ 13
  GF_mul(table_b, table_d, table_d);

  // t1 = in ^ (0xad)
  GF_sqr(table_a, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, table_d, t1);

  // t2 = in ^ (0xad 6), table_d = in ^ (0xad5)
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, table_6, t2);
  GF_mul(t1, table_5, table_d);

  // t1 = in ^ (0xad6 b)
  GF_sqr(t2, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, table_b, t1);

  // t1 = in ^ (0xad6b 5)
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, table_5, t1);

  // t1 = in ^ (0xad6b5 6)
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, table_6, t1);

  // t1 = in ^ (0xad6b56 b)
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, table_b, t1);

  // t1 = in ^ (0xad6b56b 5)
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, table_5, t1);

  // t1 = in ^ (0xad6b56b5 a)
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, table_a, t1);

  // t1 = in ^ (0xad6b56b5a b)
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, table_b, t1);

  // t1 = in ^ (0xad6b56b5ab 5)
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, table_5, t1);

  // table_d = in ^ (0xad6b56b5ab5 ad5)
  for (i = 0; i < 12; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_d, table_d);

  // t1 = n ^ (0xad6b56b5ab5ad5 ad6)
  GF_sqr(table_d, t1);
  for (i = 1; i < 12; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, t2, t1);

  // t1 = in ^ (0xad6b56b5ab5ad5ad6 ad6b56b5ab5ad5)
  for (i = 0; i < 56; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_d, t1);

  // t1 = in ^ (0xad6b56b5ab5ad5ad6ad6b56b5ab5ad5 ad6)
  for (i = 0; i < 12; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, t2, t1);

  // t1 = in ^ (0xad6b56b5ab5ad5ad6ad6b56b5ab5ad5ad6 ad6b56b5ab5ad5)
  for (i = 0; i < 56; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_d, out);
}

// inverse Mersenne S-box with e2 = 47
// (2 ^ 47 - 1) ^ (-1) mod (2 ^ 192 - 1)
// = 0xddddddddddddbbbbbbbbbbbb777777777776eeeeeeeeeeed
// dddd dddd dddd bb bb bb bb bb bb 77 77 77 77 77 76 ee ee ee ee ee ed
void GF_exp_invmer_e_2(const GF in, GF out)
{
  unsigned int i;
  GF t1 = {0,}, t2 = {0,};
  GF table_6 = {0,}, table_7 = {0,};
  GF table_b = {0,}, table_d = {0,}, table_e = {0,};

  // t1 = in ^ 3
  GF_sqr(in, table_d);
  GF_mul(table_d, in, t1);

  // table_6 = (in ^ 3) ^ 2
  GF_sqr(t1, table_6);
  // table_7 = in ^ 7
  GF_mul(table_6, in, table_7);
  // table_b = in ^ 11
  GF_sqr(table_d, table_b);
  GF_mul(table_b, table_7, table_b);
  // table_d = in ^ 13
  GF_mul(table_6, table_7, table_d);
  // table_e = in ^ 14
  GF_sqr(table_7, table_e);

  // table_b = in ^ (0xbb)
  GF_sqr(table_b, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, table_b, table_b);

  // table_7 = in ^ (0x77), table_6 = in ^ (0x76)
  GF_sqr(table_7, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, table_6, table_6);
  GF_mul(t1, table_7, table_7);

  // t2 = in ^ (0xdd)
  GF_sqr(table_d, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, table_d, t2);

  // table_e = in ^ (0xee), table_d = in ^ (0xed)
  GF_sqr(table_e, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, table_d, table_d);
  GF_mul(t1, table_e, table_e);

  // t2 = in ^ (0xdd dd)
  GF_sqr(t2, t1);
  for (i = 1; i < 8; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, t2, t2);

  // t1 = in ^ (0xdddd dddd)
  GF_sqr(t2, t1);
  for (i = 1; i < 16; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, t2, t1);

  // t1 = in ^ (0xdddddddd dddd)
  for (i = 0; i < 16; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, t2, t1);

  // t1 = in ^ (0xdddddddddddd bb)
  for (i = 0; i < 8; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_b, t1);

  // t1 = in ^ (0xddddddddddddbb bb)
  for (i = 0; i < 8; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_b, t1);

  // t1 = in ^ (0xddddddddddddbbbb bb)
  for (i = 0; i < 8; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_b, t1);

  // t1 = in ^ (0xddddddddddddbbbbbb bb)
  for (i = 0; i < 8; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_b, t1);

  // t1 = in ^ (0xddddddddddddbbbbbbbb bb)
  for (i = 0; i < 8; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_b, t1);

  // t1 = in ^ (0xddddddddddddbbbbbbbbbb bb)
  for (i = 0; i < 8; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_b, t1);

  // t1 = in ^ (0xddddddddddddbbbbbbbbbbbb 77)
  for (i = 0; i < 8; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_7, t1);

  // t1 = in ^ (0xddddddddddddbbbbbbbbbbbb77 77)
  for (i = 0; i < 8; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_7, t1);

  // t1 = in ^ (0xddddddddddddbbbbbbbbbbbb7777 77)
  for (i = 0; i < 8; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_7, t1);

  // t1 = in ^ (0xddddddddddddbbbbbbbbbbbb777777 77)
  for (i = 0; i < 8; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_7, t1);

  // t1 = in ^ (0xddddddddddddbbbbbbbbbbbb77777777 77)
  for (i = 0; i < 8; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_7, t1);

  // t1 = in ^ (0xddddddddddddbbbbbbbbbbbb7777777777 76)
  for (i = 0; i < 8; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_6, t1);

  // t1 = in ^ (0xddddddddddddbbbbbbbbbbbb777777777776 ee)
  for (i = 0; i < 8; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_e, t1);

  // t1 = in ^ (0xddddddddddddbbbbbbbbbbbb777777777776ee ee)
  for (i = 0; i < 8; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_e, t1);

  // t1 = in ^ (0xddddddddddddbbbbbbbbbbbb777777777776eeee ee)
  for (i = 0; i < 8; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_e, t1);

  // t1 = in ^ (0xddddddddddddbbbbbbbbbbbb777777777776eeeeee ee)
  for (i = 0; i < 8; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_e, t1);

  // t1 = in ^ (0xddddddddddddbbbbbbbbbbbb777777777776eeeeeeee ee)
  for (i = 0; i < 8; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_e, t1);

  // out = in ^ (0xddddddddddddbbbbbbbbbbbb777777777776eeeeeeeeee ed)
  for (i = 0; i < 8; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, table_d, out);
}

// Mersenne exponentiation with e_star = 5
void GF_exp_mer_e_star(const GF in, GF out)
{
  GF t1 = {0,};
  GF t2 = {0,};

  // t2 = a ^ (2 ^ 2 - 1)
  GF_sqr(in, t1);
  GF_mul(t1, in, t2);

  // t1 = a ^ (2 ^ 3 - 1)
  GF_sqr(t2, t1);
  GF_mul(t1, in, t1);

  // out = a ^ (2 ^ 5 - 1)
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, t2, out);
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
  GF matrix_L[AIM2_NUM_INPUT_SBOX][AIM2_NUM_BITS_FIELD];
  GF matrix_U[AIM2_NUM_INPUT_SBOX][AIM2_NUM_BITS_FIELD];

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
  GF matrix_L[AIM2_NUM_INPUT_SBOX][AIM2_NUM_BITS_FIELD];
  GF matrix_U[AIM2_NUM_INPUT_SBOX][AIM2_NUM_BITS_FIELD];
  GF vector_b = {0,};

  GF state[AIM2_NUM_INPUT_SBOX];
  GF pt_GF = {0,}, ct_GF = {0,};
  GF_from_bytes(pt, pt_GF);

  // generate random matrix
  generate_matrices_L_and_U(iv, matrix_L, matrix_U, vector_b);

  // linear component: constant addition
  GF_add(pt_GF, aim2_constants[0], state[0]);
  GF_add(pt_GF, aim2_constants[1], state[1]);

  // non-linear component: inverse Mersenne S-box
  GF_exp_invmer_e_1(state[0], state[0]);
  GF_exp_invmer_e_2(state[1], state[1]);

  // linear component: affine layer
  GF_transposed_matmul(state[0], matrix_U[0], state[0]);
  GF_transposed_matmul(state[0], matrix_L[0], state[0]);

  GF_transposed_matmul(state[1], matrix_U[1], state[1]);
  GF_transposed_matmul(state[1], matrix_L[1], state[1]);

  GF_add(state[0], state[1], state[0]);
  GF_add(state[0], vector_b, state[0]);

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

  // non-linear component: inverse Mersenne S-box
  GF_exp_invmer_e_1(sbox_outputs[0], sbox_outputs[0]);
  GF_exp_invmer_e_2(sbox_outputs[1], sbox_outputs[1]);
}
