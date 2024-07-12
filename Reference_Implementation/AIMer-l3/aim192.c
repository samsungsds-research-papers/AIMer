// -----------------------------------------------------------------------------
// File Name   : aim192.c
// Description : 
// SPDX-License-Identifier: MIT
// -----------------------------------------------------------------------------

#include "aim192.h"
#include "hash.h"

#include "aim_common.c"

void mersenne_exp_5(const GF in, GF out);
void mersenne_exp_7(const GF in, GF out);
void mersenne_exp_29(const GF in, GF out);
void GF_setbit(GF a, int n, uint64_t bit);
uint64_t GF_getbit(const GF a, int n);

void transpose_lower_to_upper(const GF* a, GF* transpose_a);
void transpose_upper_to_lower(const GF* a, GF* transpose_a);
void generate_linear_layer_separate(const uint8_t* initial_vector,
                                    GF** matrix_A, GF vector_B);

void mersenne_exp_5(const GF in, GF out)
{
  GF t1 = {0,};
  GF t2 = {0,};

  // t2 = a ^ (2^2 - 1)
  GF_sqr(in, t1);
  GF_mul(t1, in, t2);

  // t1 = a ^ (2^4 - 1)
  GF_sqr(t2, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, t2, t1);

  // out = a ^ (2^5 - 1)
  GF_sqr(t1, t1);
  GF_mul(t1, in, out);
}

void mersenne_exp_7(const GF in, GF out)
{
  GF t1 = {0,};
  GF t2 = {0,};

  // t1 = a ^ (2^2 - 1)
  GF_sqr(in, t1);
  GF_mul(t1, in, t1);

  // t2 = a ^ (2^3 - 1)
  GF_sqr(t1, t1);
  GF_mul(t1, in, t2);

  // t1 = a ^ (2^6 - 1)
  GF_sqr(t2, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, t2, t1);

  // out = a ^ (2^7 - 1)
  GF_sqr(t1, t1);
  GF_mul(t1, in, out);
}

void mersenne_exp_29(const GF in, GF out)
{
  int i;
  GF t1 = {0,};
  GF t2 = {0,};

  // t1 = a ^ (2^2 - 1)
  GF_sqr(in, t1);
  GF_mul(t1, in, t1);

  // t2 = a ^ (2^3 - 1)
  GF_sqr(t1, t1);
  GF_mul(t1, in, t2);

  // t2 = a ^ (2^6 - 1)
  GF_sqr(t2, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, t2, t2);

  // t2 = a ^ (2^7 - 1)
  GF_sqr(t2, t1);
  GF_mul(t1, in, t2);

  // t2 = a ^ (2^14 - 1)
  GF_sqr(t2, t1);
  for (i = 1; i < 7; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, t2, t2);

  // t2 = a ^ (2^28 - 1)
  GF_sqr(t2, t1);
  for (i = 1; i < 14; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, t2, t2);

  // out = a ^ (2^29 - 1)
  GF_sqr(t2, t1);
  GF_mul(t1, in, out);
}

uint64_t GF_getbit(const GF a, int n)
{
  int i, j;

  i = n / NUMBITS_WORD;
  j = n % NUMBITS_WORD;

  return (a[i] >> j) & 1;
}

void GF_setbit(GF a, int n, uint64_t bit)
{
  int i, j;

  i = n / NUMBITS_WORD;
  j = n % NUMBITS_WORD;

  a[i] |= (bit << j);
}

void transpose_lower_to_upper(const GF* a, GF* transpose_a)
{
  for (size_t i = 1; i <= NUMBITS_FIELD; i++)
  {
    for (size_t j = NUMBITS_FIELD - 1; j > NUMBITS_FIELD - i; j--)
    {
      GF_setbit(transpose_a[j], NUMBITS_FIELD - i, GF_getbit(a[NUMBITS_FIELD - i],
                j));
    }
    GF_setbit(transpose_a[NUMBITS_FIELD - i], NUMBITS_FIELD - i, 1);
  }
}

void transpose_upper_to_lower(const GF* a, GF* transpose_a)
{
  for (size_t i = 1; i <= NUMBITS_FIELD; i++)
  {
    for (size_t j = 0; j < NUMBITS_FIELD - i + 1; j++)
    {
      GF_setbit(transpose_a[j], NUMBITS_FIELD - i, GF_getbit(a[NUMBITS_FIELD - i],
                j));
    }
    GF_setbit(transpose_a[NUMBITS_FIELD - i], NUMBITS_FIELD - i, 1);
  }
}

void generate_linear_layer(const uint8_t* initial_vector, GF** matrix_A,
                           GF vector_B)
{
  hash_instance ctx;
  size_t field_size = sizeof(GF);
  size_t squeeze_len;
  size_t block_index;

  GF* lower = malloc((field_size * field_size) << 3);
  GF* upper = malloc((field_size * field_size) << 3);

  // Initialize Hash
  hash_init(&ctx, (field_size << 1));
  hash_update(&ctx, initial_vector, field_size);
  hash_final(&ctx);

  uint8_t* out = malloc(TAPE_LEN);
  uint8_t* out_ptr;

  for (size_t num = 0; num < 2; num++)
  {
    hash_squeeze(&ctx, out, TAPE_LEN);

    memset(lower, 0, (field_size * field_size) << 3);
    memset(upper, 0, (field_size * field_size) << 3);
    matrix_A[num] = calloc(field_size << 3, field_size);

    squeeze_len = 8;
    block_index = 0;

    out_ptr = out;
    for (size_t row = 0; row < NUMBITS_FIELD; row++)
    {
      uint64_t tmp = 0;
      memcpy((uint8_t *)(&tmp) + (8 - squeeze_len), out_ptr, squeeze_len);
      out_ptr += squeeze_len;
      lower[row][block_index] =
        (le64toh(tmp) & and_mask_lower[row & 0x3f]) | or_mask[row & 0x3f];
      for (size_t i = block_index + 1; i < NUMWORDS_FIELD; i++)
      {
        memcpy((uint8_t *)(&tmp), out_ptr, 8);
        out_ptr += 8;
        lower[row][i] = le64toh(tmp);
      }
      if ((row & 7) == 6)
      {
        squeeze_len--;
      }
      if ((row & 0x3f) == 0x3f)
      {
        block_index++;
        squeeze_len = 8;
      }
    }

    squeeze_len = 0;
    block_index = 0;
    for (size_t row = 0; row < NUMBITS_FIELD; row++)
    {
      uint64_t tmp = 0;
      memcpy((uint8_t *)(&tmp), out_ptr, squeeze_len);
      out_ptr += squeeze_len;
      upper[row][block_index] =
        (le64toh(tmp) & and_mask_upper[row & 0x3f]) | or_mask[row & 0x3f];
      for (size_t i = 0; i < block_index; i++)
      {
        memcpy((uint8_t *)(&tmp), out_ptr, 8);
        out_ptr += 8;
        upper[row][i] = le64toh(tmp);
      }
      if ((row & 7) == 0)
      {
        squeeze_len++;
      }
      if ((row & 0x3f) == 0x3f)
      {
        block_index++;
        squeeze_len = 0;
      }
    }

    // upper_tr -> L'
    GF* upper_tr = calloc(field_size << 3, field_size);
    transpose_upper_to_lower(upper, upper_tr);

    // L' x transpose(L) -> L' x U'
    GF_matmul_mat((const uint64_t*)upper_tr, (const uint64_t*)lower,
                  (uint64_t*)matrix_A[num]);

    free(upper_tr);
  }

  // generate vector B
  hash_squeeze(&ctx, (uint8_t*)vector_B, field_size);

  free(lower);
  free(upper);
  free(out);
}

void generate_linear_layer_separate(const uint8_t* initial_vector,
                                    GF** matrix_A, GF vector_B)
{
  hash_instance ctx;
  size_t field_size = sizeof(GF);
  size_t squeeze_len;
  size_t block_index;

  // Initialize Hash
  hash_init(&ctx, (field_size << 1));
  hash_update(&ctx, initial_vector, field_size);
  hash_final(&ctx);

  uint8_t* out = malloc(TAPE_LEN);
  uint8_t* out_ptr;

  for (size_t num = 0; num < 2; num++)
  {
    hash_squeeze(&ctx, out, TAPE_LEN);

    // Lower
    matrix_A[2 * num] = calloc(field_size << 3, field_size);

    squeeze_len = 8;
    block_index = 0;

    out_ptr = out;
    for (size_t row = 0; row < NUMBITS_FIELD; row++)
    {
      uint64_t tmp = 0;
      memcpy((uint8_t *)(&tmp) + (8 - squeeze_len), out_ptr, squeeze_len);
      out_ptr += squeeze_len;
      matrix_A[2 * num][row][block_index] =
        (le64toh(tmp) & and_mask_lower[row & 0x3f]) | or_mask[row & 0x3f];
      for (size_t i = block_index + 1; i < NUMWORDS_FIELD; i++)
      {
        memcpy((uint8_t *)(&tmp), out_ptr, 8);
        out_ptr += 8;
        matrix_A[2 * num][row][i] = le64toh(tmp);
      }
      if ((row & 7) == 6)
      {
        squeeze_len--;
      }
      if ((row & 0x3f) == 0x3f)
      {
        block_index++;
        squeeze_len = 8;
      }
    }

    // Upper
    matrix_A[2 * num + 1] = calloc(field_size << 3, field_size);

    squeeze_len = 0;
    block_index = 0;
    for (size_t row = 0; row < NUMBITS_FIELD; row++)
    {
      uint64_t tmp = 0;
      memcpy((uint8_t *)(&tmp), out_ptr, squeeze_len);
      out_ptr += squeeze_len;
      matrix_A[2 * num + 1][row][block_index] =
        (le64toh(tmp) & and_mask_upper[row & 0x3f]) | or_mask[row & 0x3f];

      for (size_t i = 0; i < block_index; i++)
      {
        memcpy((uint8_t *)(&tmp), out_ptr, 8);
        out_ptr += 8;
        matrix_A[2 * num + 1][row][i] = le64toh(tmp);
      }
      if ((row & 7) == 0)
      {
        squeeze_len++;
      }
      if ((row & 0x3f) == 0x3f)
      {
        block_index++;
        squeeze_len = 0;
      }
    }
  }

  // generate vector B
  hash_squeeze(&ctx, (uint8_t*)vector_B, field_size);

  free(out);
}

void aim192(const uint8_t* input, const uint8_t* initial_vector,
            uint8_t* output)
{
  GF state[2] = {0,};
  GF input_GF = {0,};
  GF vector_B = {0,};

  size_t field_size = sizeof(GF);

  GF **matrix_A = malloc(4 * sizeof(GF*));

  // input convert to field element
  GF_from_bytes(input, input_GF);

  // Mersenne layer
  mersenne_exp_5(input_GF, state[0]);
  mersenne_exp_29(input_GF, state[1]);

  // generate linear layer
  generate_linear_layer_separate(initial_vector, matrix_A, vector_B);

  GF* state0_lower_tr = calloc(field_size << 3, field_size);
  transpose_lower_to_upper(matrix_A[0], state0_lower_tr);
  GF* state0_upper_tr = calloc(field_size << 3, field_size);
  transpose_upper_to_lower(matrix_A[1], state0_upper_tr);

  GF* state1_lower_tr = calloc(field_size << 3, field_size);
  transpose_lower_to_upper(matrix_A[2], state1_lower_tr);
  GF* state1_upper_tr = calloc(field_size << 3, field_size);
  transpose_upper_to_lower(matrix_A[3], state1_upper_tr);

  // linear layer
  GF_matmul_vec(state[0], (uint64_t*)state0_lower_tr, state[0]);
  GF_matmul_vec(state[0], (uint64_t*)state0_upper_tr, state[0]);

  GF_matmul_vec(state[1], (uint64_t*)state1_lower_tr, state[1]);
  GF_matmul_vec(state[1], (uint64_t*)state1_upper_tr, state[1]);

  GF_add(state[0], state[1], state[0]);
  GF_add(state[0], vector_B, state[0]);

  // Mersenne layer
  mersenne_exp_7(state[0], state[0]);

  // feedback
  GF_add(state[0], input_GF, state[0]);

  // output conversion
  GF_to_bytes(state[0], output);

  for (size_t i = 0; i < 4; i++)
  {
    free(matrix_A[i]);
  }
  free(matrix_A);

  free(state0_lower_tr);
  free(state0_upper_tr);
  free(state1_lower_tr);
  free(state1_upper_tr);
}

void aim192_with_sbox_output(const uint8_t* input, const GF** matrix_A,
                             const GF vector_B, uint8_t* output,
                             GF* sbox_outputs)
{
  GF state[2] = {0,};
  GF input_GF = {0,};

  // input convert to field element
  GF_from_bytes(input, input_GF);

  // Mersenne layer
  mersenne_exp_5(input_GF, state[0]);
  mersenne_exp_29(input_GF, state[1]);
  GF_copy(state[0], sbox_outputs[0]);
  GF_copy(state[1], sbox_outputs[1]);

  // linear layer
  GF_matmul_vec(state[0], (uint64_t*)matrix_A[0], state[0]);
  GF_matmul_vec(state[1], (uint64_t*)matrix_A[1], state[1]);

  GF_add(state[0], state[1], state[0]);
  GF_add(state[0], vector_B, state[0]);

  // Mersenne layer
  mersenne_exp_7(state[0], state[0]);
  GF_copy(state[0], sbox_outputs[2]);

  // feedback
  GF_add(state[0], input_GF, state[0]);

  // output conversion
  GF_to_bytes(state[0], output);
}

void aim192_mpc(const uint8_t* input, const GF** matrix_A,
                const GF vector_B, const uint8_t* output,
                const size_t num_parties, GF *t_shares, GF *s_shares)
{
  size_t party_shares_index;
  size_t field_size = sizeof(GF);

  GF state[2] = {0,};
  GF output_GF = {0,};
  GF temp = {0,};

  GF_from_bytes(output, output_GF);

  for (size_t party = 0; party < num_parties; party++)
  {
    party_shares_index = party * 3;

    GF_from_bytes(input + party * field_size, s_shares[party_shares_index]);

    // Mersenne layer
    GF_copy(t_shares[party_shares_index], state[0]);
    GF_copy(t_shares[party_shares_index + 1], state[1]);

    // linear layer
    GF_matmul_vec(state[0], (uint64_t*)matrix_A[0], state[0]);
    GF_matmul_vec(state[1], (uint64_t*)matrix_A[1], state[1]);

    if (party == 0)
    {
      GF_add(state[1], vector_B, state[1]);
    }
    GF_add(state[0], state[1], s_shares[party_shares_index + 2]);

    // optimization B.1.2 in eprint 2022/588
    // t1 * pt = s1^(2^e1)
    // t2 * pt = s2^(2^e2)
    // s3 * pt = s3^(2^e3) - s3 * output

    // (s1, t1) -> (t1, s1^(2^e1))
    // (s2, t2) -> (t2, s2^(2^e2))
    GF_copy(t_shares[party_shares_index], state[0]);
    GF_copy(t_shares[party_shares_index + 1], state[1]);

    GF_sqr(s_shares[party_shares_index], temp);
    GF_sqr(temp, temp);
    GF_sqr(temp, temp);
    GF_sqr(temp, temp);
    GF_sqr(temp, t_shares[party_shares_index]);

    GF_sqr(t_shares[party_shares_index], temp);
    for (size_t i = 6; i < 28; i++)
    {
      GF_sqr(temp, temp);
    }
    GF_sqr(temp, t_shares[party_shares_index + 1]);

    GF_copy(state[0], s_shares[party_shares_index]);
    GF_copy(state[1], s_shares[party_shares_index + 1]);

    // (s3, t3) -> (s3, s3^(2^e3) - s3 * output)
    GF_sqr(s_shares[party_shares_index + 2], temp);
    GF_sqr(temp, temp);
    GF_sqr(temp, temp);
    GF_sqr(temp, temp);
    GF_sqr(temp, temp);
    GF_sqr(temp, temp);
    GF_sqr(temp, temp);

    GF_mul(output_GF, s_shares[party_shares_index + 2], state[0]);
    GF_add(temp, state[0], t_shares[party_shares_index + 2]);
  }
}
