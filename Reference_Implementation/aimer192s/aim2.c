// SPDX-License-Identifier: MIT

#include <stdint.h>
#include "aim2.h"
#include "field.h"
#include "hash.h"
#include "params.h"

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
  GF_exp(state[0], state[0], aim2_exponents[0]);
  GF_exp(state[1], state[1], aim2_exponents[1]);

  // linear component: affine layer
  GF_transposed_matmul(state[0], matrix_U[0], state[0]);
  GF_transposed_matmul(state[0], matrix_L[0], state[0]);

  GF_transposed_matmul(state[1], matrix_U[1], state[1]);
  GF_transposed_matmul(state[1], matrix_L[1], state[1]);

  GF_add(state[0], state[1], state[0]);
  GF_add(state[0], vector_b, state[0]);

  // non-linear component: Mersenne S-box
  GF_exp(state[0], state[0], aim2_exponents[2]);

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
  GF_exp(sbox_outputs[0], sbox_outputs[0], aim2_exponents[0]);
  GF_exp(sbox_outputs[1], sbox_outputs[1], aim2_exponents[1]);
}
