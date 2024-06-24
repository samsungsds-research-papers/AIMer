// SPDX-License-Identifier: MIT

#ifndef SIGN_H
#define SIGN_H

#include "field.h"
#include "hash.h"
#include "params.h"
#include <stddef.h>
#include <stdint.h>

typedef struct tape_t
{
  GF pt_share;
  GF t_shares[AIMER_L];
  GF a_share;
  GF c_share;
} tape_t;

typedef struct mult_chk_N_t
{
  GF pt_share[AIMER_N];
  GF x_shares[AIMER_L + 1][AIMER_N];
  GF z_shares[AIMER_L + 1][AIMER_N];
} mult_chk_N_t;

typedef struct proof_t
{
  uint8_t reveal_path[AIMER_LOGN][AIMER_SEED_SIZE];
  uint8_t missing_commitment[AIMER_COMMIT_SIZE];
  uint8_t delta_pt_bytes[AIM2_NUM_BYTES_FIELD];
  uint8_t delta_ts_bytes[AIMER_L][AIM2_NUM_BYTES_FIELD];
  uint8_t delta_c_bytes[AIM2_NUM_BYTES_FIELD];
  uint8_t missing_alpha_share_bytes[AIM2_NUM_BYTES_FIELD];
} proof_t;

typedef struct signature_t
{
  uint8_t salt[AIMER_SALT_SIZE];
  uint8_t h_1[AIMER_COMMIT_SIZE];
  uint8_t h_2[AIMER_COMMIT_SIZE];
  proof_t proofs[AIMER_T];
} signature_t;

#define aim2_mpc_N AIMER_NAMESPACE(aim2_mpc_N)
void aim2_mpc_N(mult_chk_N_t *mult_chk,
                const GF matrix_A[AIMER_L][AIM2_NUM_BITS_FIELD], const GF ct_GF);

#define commit_and_expand_tape_x4 AIMER_NAMESPACE(commit_and_expand_tape_x4)
void commit_and_expand_tape_x4(tape_t *tapes, uint8_t *commits,
                               const hash_instance_x4 *ctx_precom,
                               const uint8_t *seeds, size_t rep, size_t party);

#define run_phase_1 AIMER_NAMESPACE(run_phase_1)
void run_phase_1(signature_t *sign,
                 uint8_t commits[AIMER_T][AIMER_N][AIMER_COMMIT_SIZE],
                 uint8_t nodes[AIMER_T][2 * AIMER_N - 1][AIMER_SEED_SIZE],
                 mult_chk_N_t mult_chk[AIMER_T],
                 GF alpha_v_shares[AIMER_T][2][AIMER_N],
                 const uint8_t *sk, const uint8_t *m, size_t mlen);

#define run_phase_2_and_3 AIMER_NAMESPACE(run_phase_2_and_3)
void run_phase_2_and_3(signature_t *sign,
                       GF alpha_v_shares[AIMER_T][2][AIMER_N],
                       const mult_chk_N_t mult_chk[AIMER_T]);

#endif // SIGN_H
