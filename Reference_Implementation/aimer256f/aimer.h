// SPDX-License-Identifier: MIT

#ifndef AIMER_H
#define AIMER_H

#include <stdint.h>
#include <stddef.h>
#include "field.h"
#include "params.h"

typedef struct tape_t
{
  GF pt_share;
  GF t_shares[AIM2_NUM_INPUT_SBOX];
  GF a_share;
  GF c_share;
} tape_t;

typedef struct mult_chk_t
{
  GF pt_share;
  GF x_shares[AIM2_NUM_INPUT_SBOX + 1];
  GF z_shares[AIM2_NUM_INPUT_SBOX + 1];
  GF a_share;
  GF c_share;
} mult_chk_t;

typedef struct proof_t
{
  uint8_t reveal_path[AIMER_NUM_MPC_PARTIES_LOG2][AIMER_SEED_SIZE];
  uint8_t missing_commitment[AIMER_COMMIT_SIZE];
  uint8_t delta_pt_bytes[AIM2_NUM_BYTES_FIELD];
  uint8_t delta_ts_bytes[AIM2_NUM_INPUT_SBOX][AIM2_NUM_BYTES_FIELD];
  uint8_t delta_c_bytes[AIM2_NUM_BYTES_FIELD];
  uint8_t missing_alpha_share_bytes[AIM2_NUM_BYTES_FIELD];
} proof_t;

typedef struct signature_t
{
  uint8_t salt[AIMER_SALT_SIZE];
  uint8_t h_1[AIMER_COMMIT_SIZE];
  uint8_t h_2[AIMER_COMMIT_SIZE];
  proof_t proofs[AIMER_NUM_REPETITIONS];
} signature_t;

int aimer_sign(const uint8_t* public_key, const uint8_t* private_key,
               const uint8_t* message, const size_t message_len,
               uint8_t* signature);

int aimer_verify(const uint8_t* public_key, const uint8_t* signature,
                 const uint8_t* message, const size_t message_len);

#endif // AIMER_H
