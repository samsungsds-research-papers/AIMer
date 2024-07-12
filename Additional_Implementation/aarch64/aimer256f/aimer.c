// SPDX-License-Identifier: MIT

#include <stdlib.h>
#include <stdalign.h>
#include "aimer.h"
#include "aim2.h"
#include "field.h"
#include "hash.h"
#include "params.h"
#include "rng.h"
#include "tree.h"

#define L   AIM2_NUM_INPUT_SBOX
#define TAU AIMER_NUM_REPETITIONS
#define N   AIMER_NUM_MPC_PARTIES

void commit_to_seed_and_expand_tape_x4(const hash_instance_x4* ctx_precom,
                                       const uint8_t* seeds,
                                       size_t rep, size_t party,
                                       uint8_t* commits, tape_t* tapes)
{
  hash_instance_x4 ctx;
  uint8_t bufs[4][AIMER_SEED_SIZE + 2];
  const uint8_t* in_ptrs[4] = {bufs[0], bufs[1], bufs[2], bufs[3]};
  uint8_t* out_ptrs[4];

  memcpy(&ctx, ctx_precom, sizeof(hash_instance_x4));

  bufs[0][0] = (uint8_t)(rep);
  bufs[1][0] = (uint8_t)(rep);
  bufs[2][0] = (uint8_t)(rep);
  bufs[3][0] = (uint8_t)(rep);

  bufs[0][1] = (uint8_t)(party + 0);
  bufs[1][1] = (uint8_t)(party + 1);
  bufs[2][1] = (uint8_t)(party + 2);
  bufs[3][1] = (uint8_t)(party + 3);

  memcpy(&bufs[0][2], seeds + 0 * AIMER_SEED_SIZE, AIMER_SEED_SIZE);
  memcpy(&bufs[1][2], seeds + 1 * AIMER_SEED_SIZE, AIMER_SEED_SIZE);
  memcpy(&bufs[2][2], seeds + 2 * AIMER_SEED_SIZE, AIMER_SEED_SIZE);
  memcpy(&bufs[3][2], seeds + 3 * AIMER_SEED_SIZE, AIMER_SEED_SIZE);

  hash_update_x4(&ctx, in_ptrs, AIMER_SEED_SIZE + 2);
  hash_final_x4(&ctx);

  out_ptrs[0] = commits + 0 * AIMER_COMMIT_SIZE;
  out_ptrs[1] = commits + 1 * AIMER_COMMIT_SIZE;
  out_ptrs[2] = commits + 2 * AIMER_COMMIT_SIZE;
  out_ptrs[3] = commits + 3 * AIMER_COMMIT_SIZE;

  hash_squeeze_x4(&ctx, out_ptrs, AIMER_COMMIT_SIZE);

  out_ptrs[0] = (uint8_t*)(tapes);
  out_ptrs[1] = (uint8_t*)(tapes + 1);
  out_ptrs[2] = (uint8_t*)(tapes + 2);
  out_ptrs[3] = (uint8_t*)(tapes + 3);

  hash_squeeze_x4(&ctx, out_ptrs, sizeof(tape_t));
}

void aim2_mpc_N(const GF matrix_A[][AIM2_NUM_BITS_FIELD],
                const GF ct, mult_chk_N_t* mult_chk)
{
  // pt + c = t ^ {2 ^ e - 1}
  // --> t ^ {2 ^ e} + t * c = t * pt
  // --> z = x * pt
  GF_mul_N(mult_chk->x_shares[0], aim2_constants[0], mult_chk->z_shares[0]);
  GF_transposed_matmul_add_N(mult_chk->x_shares[0], aim2_e1_power_matrix,
                             mult_chk->z_shares[0]);
  GF_transposed_matmul_add_N(mult_chk->x_shares[0], matrix_A[0],
                             mult_chk->x_shares[L]);

  GF_mul_N(mult_chk->x_shares[1], aim2_constants[1], mult_chk->z_shares[1]);
  GF_transposed_matmul_add_N(mult_chk->x_shares[1], aim2_e2_power_matrix,
                             mult_chk->z_shares[1]);
  GF_transposed_matmul_add_N(mult_chk->x_shares[1], matrix_A[1],
                             mult_chk->x_shares[L]);

  GF_sqr_N(mult_chk->x_shares[2], mult_chk->z_shares[2]);
  GF_sqr_N(mult_chk->z_shares[2], mult_chk->z_shares[2]);
  GF_sqr_N(mult_chk->z_shares[2], mult_chk->z_shares[2]);
  GF_sqr_N(mult_chk->z_shares[2], mult_chk->z_shares[2]);
  GF_sqr_N(mult_chk->z_shares[2], mult_chk->z_shares[2]);
  GF_sqr_N(mult_chk->z_shares[2], mult_chk->z_shares[2]);
  GF_sqr_N(mult_chk->z_shares[2], mult_chk->z_shares[2]);
  GF_mul_add_N(mult_chk->x_shares[2], aim2_constants[2], mult_chk->z_shares[2]);
  GF_transposed_matmul_add_N(mult_chk->x_shares[2], matrix_A[2],
                             mult_chk->x_shares[L]);

  // x ^ {2 ^ e - 1} = pt + ct
  // --> x ^ {2 ^ e} + x * ct = x * pt
  // --> z = x * pt
  GF_sqr_N(mult_chk->x_shares[L], mult_chk->z_shares[L]);
  GF_sqr_N(mult_chk->z_shares[L], mult_chk->z_shares[L]);
  GF_sqr_N(mult_chk->z_shares[L], mult_chk->z_shares[L]);
  GF_mul_add_N(mult_chk->x_shares[L], ct, mult_chk->z_shares[L]);
}

// committing to the seeds and the execution views of the parties
void run_phase_1(const uint8_t* public_key, const uint8_t* private_key,
                 const uint8_t* message, const size_t message_len,
                 signature_t* sig, uint8_t commits[TAU][N][AIMER_COMMIT_SIZE],
                 uint8_t nodes[TAU][2 * N - 1][AIMER_SEED_SIZE],
                 mult_chk_N_t mult_chk[TAU],
                 GF alpha_v_shares[TAU][2][N],
                 uint8_t h_1[AIMER_COMMIT_SIZE])
{
  alignas(32) GF pt = {0,}, ct = {0,};
  GF_from_bytes(private_key, pt);
  GF_from_bytes(public_key + AIM2_IV_SIZE, ct);

  // message pre-hashing
  hash_instance ctx;
  hash_init_prefix(&ctx, HASH_PREFIX_0);
  hash_update(&ctx, public_key, AIM2_IV_SIZE + AIM2_NUM_BYTES_FIELD);
  hash_update(&ctx, message, message_len);
  hash_final(&ctx);

  uint8_t mu[AIMER_COMMIT_SIZE];
  hash_squeeze(&ctx, mu, AIMER_COMMIT_SIZE);

  // compute first L sboxes' outputs
  alignas(32) GF sbox_outputs[AIM2_NUM_INPUT_SBOX];
  aim2_sbox_outputs(pt, sbox_outputs);

  // derive the binary matrix and the vector from the initial vector
  alignas(32) GF matrix_A[AIM2_NUM_INPUT_SBOX][AIM2_NUM_BITS_FIELD];
  alignas(32) GF vector_b = {0,};
  generate_matrix_LU(public_key, matrix_A, vector_b); // iv is at public_key

  // generate per-signature randomness
  uint8_t random[SECURITY_BYTES];
  randombytes(random, SECURITY_BYTES);

  // generate salt
  hash_init_prefix(&ctx, HASH_PREFIX_3);
  hash_update(&ctx, private_key, AIM2_NUM_BYTES_FIELD);
  hash_update(&ctx, mu, AIMER_COMMIT_SIZE);
  hash_update(&ctx, random, SECURITY_BYTES);
  hash_final(&ctx);
  hash_squeeze(&ctx, sig->salt, AIMER_SALT_SIZE);

  // generate root seeds and expand seed trees
  for (size_t rep = 0; rep < TAU; rep++)
  {
    hash_squeeze(&ctx, nodes[rep][0], AIMER_SEED_SIZE);
  }
  expand_trees(sig->salt, nodes);

  // hash_instance for h_1
  hash_init_prefix(&ctx, HASH_PREFIX_1);
  hash_update(&ctx, mu, AIMER_COMMIT_SIZE);
  hash_update(&ctx, sig->salt, AIMER_SALT_SIZE);

  hash_instance_x4 ctx_precom;
  hash_init_prefix_x4(&ctx_precom, HASH_PREFIX_5);
  hash_update_x4_1(&ctx_precom, sig->salt, AIMER_SALT_SIZE);

  for (size_t rep = 0; rep < TAU; rep++)
  {
    // initialize adjustment values
    alignas(32) tape_t tapes[4];
    alignas(32) tape_t delta;
    memset(&delta, 0, sizeof(tape_t));

    // initialize x_star
    memset(mult_chk[rep].x_shares[L], 0, sizeof(GF) * (N - 1));
    GF_copy(vector_b, mult_chk[rep].x_shares[L][N - 1]);

    for (size_t party = 0; party < N; party += 4)
    {
      // generate execution views and commitments
      commit_to_seed_and_expand_tape_x4(&ctx_precom, nodes[rep][party + N - 1],
                                        rep, party, commits[rep][party],
                                        tapes);
      hash_update(&ctx, commits[rep][party], 4 * AIMER_COMMIT_SIZE);

      for (size_t j = 0; j < 4; j++)
      {
        // compute offsets
        GF_add(delta.pt_share, tapes[j].pt_share, delta.pt_share);
        GF_add(delta.t_shares[0], tapes[j].t_shares[0], delta.t_shares[0]);
        GF_add(delta.t_shares[1], tapes[j].t_shares[1], delta.t_shares[1]);
        GF_add(delta.t_shares[2], tapes[j].t_shares[2], delta.t_shares[2]);
        GF_add(delta.a_share, tapes[j].a_share, delta.a_share);
        GF_add(delta.c_share, tapes[j].c_share, delta.c_share);

        if (party + j == N - 1)
        {
          // adjust the last share and prepare the proof and h_1
          GF_add(delta.pt_share, pt, delta.pt_share);
          GF_add(delta.t_shares[0], sbox_outputs[0], delta.t_shares[0]);
          GF_add(delta.t_shares[1], sbox_outputs[1], delta.t_shares[1]);
          GF_add(delta.t_shares[2], sbox_outputs[2], delta.t_shares[2]);
          GF_mul_add(pt, delta.a_share, delta.c_share);

          GF_to_bytes(delta.pt_share, sig->proofs[rep].delta_pt_bytes);
          GF_to_bytes(delta.t_shares[0], sig->proofs[rep].delta_ts_bytes[0]);
          GF_to_bytes(delta.t_shares[1], sig->proofs[rep].delta_ts_bytes[1]);
          GF_to_bytes(delta.t_shares[2], sig->proofs[rep].delta_ts_bytes[2]);
          GF_to_bytes(delta.c_share, sig->proofs[rep].delta_c_bytes);

          GF_add(delta.pt_share, tapes[j].pt_share, tapes[j].pt_share);
          GF_add(delta.t_shares[0], tapes[j].t_shares[0], tapes[j].t_shares[0]);
          GF_add(delta.t_shares[1], tapes[j].t_shares[1], tapes[j].t_shares[1]);
          GF_add(delta.t_shares[2], tapes[j].t_shares[2], tapes[j].t_shares[2]);
          GF_add(delta.c_share, tapes[j].c_share, tapes[j].c_share);
        }

        // run the MPC simulation and prepare the mult check inputs
        GF_copy(tapes[j].pt_share, mult_chk[rep].pt_share[party + j]);
        GF_copy(tapes[j].t_shares[0], mult_chk[rep].x_shares[0][party + j]);
        GF_copy(tapes[j].t_shares[1], mult_chk[rep].x_shares[1][party + j]);
        GF_copy(tapes[j].t_shares[2], mult_chk[rep].x_shares[2][party + j]);
        GF_copy(tapes[j].a_share, alpha_v_shares[rep][0][party + j]);
        GF_copy(tapes[j].c_share, alpha_v_shares[rep][1][party + j]);
      }
    }
    aim2_mpc_N(matrix_A, ct, &mult_chk[rep]);

    // NOTE: depend on the order of values in proof_t
    hash_update(&ctx, sig->proofs[rep].delta_pt_bytes,
                AIM2_NUM_BYTES_FIELD * (L + 2));
  }

  // commit to salt, (all commitments of parties' seeds,
  // delta_pt, delta_t, delta_c) for all repetitions
  hash_final(&ctx);
  hash_squeeze(&ctx, h_1, AIMER_COMMIT_SIZE);
}

void run_phase_2_and_3(const uint8_t h_1[AIMER_COMMIT_SIZE],
                       const uint8_t salt[AIMER_SALT_SIZE],
                       const mult_chk_N_t mult_chk[TAU],
                       GF alpha_v_shares[TAU][2][N],
                       uint8_t h_2[AIMER_COMMIT_SIZE])
{
  hash_instance ctx_e;
  hash_init(&ctx_e);
  hash_update(&ctx_e, h_1, AIMER_COMMIT_SIZE);
  hash_final(&ctx_e);

  hash_instance ctx;
  hash_init_prefix(&ctx, HASH_PREFIX_2);
  hash_update(&ctx, h_1, AIMER_COMMIT_SIZE);
  hash_update(&ctx, salt, AIMER_SALT_SIZE);

  for (size_t rep = 0; rep < TAU; rep++)
  {
    alignas(32) GF epsilons[L + 1];
    hash_squeeze(&ctx_e, (uint8_t*)epsilons, AIM2_NUM_BYTES_FIELD * (L + 1));

    alignas(32) GF alpha = {0,};
    alignas(32) GF alpha_v_shares_hi[2][N];
    memset(alpha_v_shares_hi, 0, AIM2_NUM_BYTES_FIELD * 2 * N);

    POLY_mul_add_N(mult_chk[rep].x_shares[0], epsilons[0],
                   alpha_v_shares[rep][0], alpha_v_shares_hi[0]);
    POLY_mul_add_N(mult_chk[rep].z_shares[0], epsilons[0],
                   alpha_v_shares[rep][1], alpha_v_shares_hi[1]);

    POLY_mul_add_N(mult_chk[rep].x_shares[1], epsilons[1],
                   alpha_v_shares[rep][0], alpha_v_shares_hi[0]);
    POLY_mul_add_N(mult_chk[rep].z_shares[1], epsilons[1],
                   alpha_v_shares[rep][1], alpha_v_shares_hi[1]);

    POLY_mul_add_N(mult_chk[rep].x_shares[2], epsilons[2],
                   alpha_v_shares[rep][0], alpha_v_shares_hi[0]);
    POLY_mul_add_N(mult_chk[rep].z_shares[2], epsilons[2],
                   alpha_v_shares[rep][1], alpha_v_shares_hi[1]);

    POLY_mul_add_N(mult_chk[rep].x_shares[3], epsilons[3],
                   alpha_v_shares[rep][0], alpha_v_shares_hi[0]);
    POLY_mul_add_N(mult_chk[rep].z_shares[3], epsilons[3],
                   alpha_v_shares[rep][1], alpha_v_shares_hi[1]);

    POLY_red_N(alpha_v_shares_hi[0], alpha_v_shares[rep][0]);
    POLY_red_N(alpha_v_shares_hi[1], alpha_v_shares[rep][1]);

    for (size_t party = 0; party < N; party++)
    {
      GF_add(alpha, alpha_v_shares[rep][0][party], alpha);
    }

    // alpha is opened, so we can finish calculating v_share
    GF_mul_add_N(mult_chk[rep].pt_share, alpha, alpha_v_shares[rep][1]);
    hash_update(&ctx, (const uint8_t*)alpha_v_shares[rep],
                AIM2_NUM_BYTES_FIELD * 2 * N);
  }
  hash_final(&ctx);
  hash_squeeze(&ctx, h_2, AIMER_COMMIT_SIZE);
}

int aimer_sign(const uint8_t* public_key, const uint8_t* private_key,
               const uint8_t* message, const size_t message_len,
               uint8_t* signature)
{
  int ret = 0;

  hash_instance ctx;

  signature_t* sig = (signature_t*)signature;

  //////////////////////////////////////////////////////////////////////////
  // Phase 1: Committing to the seeds and the execution views of parties. //
  //////////////////////////////////////////////////////////////////////////

  // nodes for seed trees
  uint8_t (*nodes)[2 * N - 1][AIMER_SEED_SIZE] = malloc(sizeof(*nodes) * TAU);
  memset(nodes, 0, TAU * sizeof(*nodes));

  // commitments for seeds
  uint8_t (*commits)[N][AIMER_COMMIT_SIZE] = malloc(sizeof(*commits) * TAU);
  memset(commits, 0, TAU * sizeof(*commits));

  // multiplication check inputs
  mult_chk_N_t* mult_chk = aligned_alloc(32, sizeof(*mult_chk) * TAU);
  memset(mult_chk, 0, TAU * sizeof(*mult_chk));

  // multiplication check outputs
  GF (*alpha_v_shares)[2][N] = aligned_alloc(32, sizeof(*alpha_v_shares) * TAU);
  memset(alpha_v_shares, 0, TAU * sizeof(*alpha_v_shares));

  run_phase_1(public_key, private_key, message, message_len,
              sig, commits, nodes, mult_chk, alpha_v_shares, sig->h_1);

  /////////////////////////////////////////////////////////////////
  // Phase 2, 3: Challenging and committing to the simulation of //
  //             the multiplication checking protocol.           //
  /////////////////////////////////////////////////////////////////

  // compute the commitment of phase 3
  run_phase_2_and_3(sig->h_1, sig->salt, mult_chk, alpha_v_shares, sig->h_2);

  //////////////////////////////////////////////////////
  // Phase 4: Challenging views of the MPC protocols. //
  //////////////////////////////////////////////////////

  hash_init(&ctx);
  hash_update(&ctx, sig->h_2, AIMER_COMMIT_SIZE);
  hash_final(&ctx);

  uint8_t indices[TAU]; // N <= 256
  hash_squeeze(&ctx, indices, TAU);
  for (size_t rep = 0; rep < TAU; rep++)
  {
    indices[rep] &= (1 << AIMER_NUM_MPC_PARTIES_LOG2) - 1;
  }

  //////////////////////////////////////////////////////
  // Phase 5: Opening the views of the MPC protocols. //
  //////////////////////////////////////////////////////

  for (size_t rep = 0; rep < TAU; rep++)
  {
    size_t i_bar = indices[rep];
    reveal_all_but(nodes[rep], i_bar, sig->proofs[rep].reveal_path);
    memcpy(sig->proofs[rep].missing_commitment, commits[rep][i_bar],
           AIMER_COMMIT_SIZE);
    GF_to_bytes(alpha_v_shares[rep][0][i_bar],
                sig->proofs[rep].missing_alpha_share_bytes);
  }

  free(nodes);
  free(commits);
  free(mult_chk);
  free(alpha_v_shares);

  return ret;
}

int aimer_verify(const uint8_t* public_key, const uint8_t* signature,
                 const uint8_t* message, const size_t message_len)
{
  int ret = 0;
  const signature_t* sig = (const signature_t*)signature;

  uint8_t iv[AIM2_IV_SIZE];
  memcpy(iv, public_key, AIM2_IV_SIZE);

  alignas(32) GF ct = {0,};
  GF_from_bytes(public_key + AIM2_IV_SIZE, ct);

  // derive the binary matrix and the vector from the initial vector
  alignas(32) GF matrix_A[AIM2_NUM_INPUT_SBOX][AIM2_NUM_BITS_FIELD];
  alignas(32) GF vector_b = {0,};
  generate_matrix_LU(iv, matrix_A, vector_b);

  hash_instance ctx_e, ctx_h1, ctx_h2;

  // indices = Expand(h_2)
  hash_init(&ctx_e);
  hash_update(&ctx_e, sig->h_2, AIMER_COMMIT_SIZE);
  hash_final(&ctx_e);

  uint8_t indices[TAU]; // N <= 256
  hash_squeeze(&ctx_e, indices, TAU);
  for (size_t rep = 0; rep < TAU; rep++)
  {
    indices[rep] &= (1 << AIMER_NUM_MPC_PARTIES_LOG2) - 1;
  }

  // epsilons = Expand(h_1)
  hash_init(&ctx_e);
  hash_update(&ctx_e, sig->h_1, AIMER_COMMIT_SIZE);
  hash_final(&ctx_e);

  // message pre-hashing
  uint8_t mu[AIMER_COMMIT_SIZE];
  hash_init_prefix(&ctx_h1, HASH_PREFIX_0);
  hash_update(&ctx_h1, public_key, AIM2_IV_SIZE + AIM2_NUM_BYTES_FIELD);
  hash_update(&ctx_h1, message, message_len);
  hash_final(&ctx_h1);
  hash_squeeze(&ctx_h1, mu, AIMER_COMMIT_SIZE);

  // ready for computing h_1' and h_2'
  hash_init_prefix(&ctx_h1, HASH_PREFIX_1);
  hash_update(&ctx_h1, mu, AIMER_COMMIT_SIZE);
  hash_update(&ctx_h1, sig->salt, AIMER_SALT_SIZE);

  hash_init_prefix(&ctx_h2, HASH_PREFIX_2);
  hash_update(&ctx_h2, sig->h_1, AIMER_COMMIT_SIZE);
  hash_update(&ctx_h2, sig->salt, AIMER_SALT_SIZE);

  hash_instance_x4 ctx_precom;
  hash_init_prefix_x4(&ctx_precom, HASH_PREFIX_5);
  hash_update_x4_1(&ctx_precom, sig->salt, AIMER_SALT_SIZE);

  for (size_t rep = 0; rep < TAU; rep++)
  {
    uint8_t nodes[2 * AIMER_NUM_MPC_PARTIES - 2][AIMER_SEED_SIZE] = {0,};
    size_t i_bar = indices[rep];
    reconstruct_seed_tree(sig->proofs[rep].reveal_path, i_bar,
                          sig->salt, rep, nodes);

    alignas(32) mult_chk_N_t mult_chk;
    memset(mult_chk.x_shares[L], 0, sizeof(GF) * (N - 1));
    GF_copy(vector_b, mult_chk.x_shares[L][N - 1]);

    alignas(32) GF epsilons[L + 1];
    hash_squeeze(&ctx_e, (uint8_t*)epsilons, AIM2_NUM_BYTES_FIELD * (L + 1));

    alignas(32) GF alpha_v = {0,};
    alignas(32) GF alpha_v_shares[2][N], alpha_v_shares_hi[2][N];
    memset(alpha_v_shares_hi, 0, AIM2_NUM_BYTES_FIELD * 2 * N);

    for (size_t party = 0; party < N; party += 4)
    {
      tape_t tapes[4];
      uint8_t commits[4][AIMER_COMMIT_SIZE];
      commit_to_seed_and_expand_tape_x4(&ctx_precom, nodes[party + N - 2],
                                        rep, party, commits[0],
                                        tapes);

      if (party / 4 == i_bar / 4)
      {
        memcpy(commits[i_bar % 4], sig->proofs[rep].missing_commitment,
               AIMER_COMMIT_SIZE);
      }
      hash_update(&ctx_h1, commits[0], 4 * AIMER_COMMIT_SIZE);

      if (party == N - 4)
      {
        GF_add(tapes[3].pt_share, (uint64_t*)sig->proofs[rep].delta_pt_bytes,
               tapes[3].pt_share);
        GF_add(tapes[3].t_shares[0], (uint64_t*)sig->proofs[rep].delta_ts_bytes[0],
               tapes[3].t_shares[0]);
        GF_add(tapes[3].t_shares[1], (uint64_t*)sig->proofs[rep].delta_ts_bytes[1],
               tapes[3].t_shares[1]);
        GF_add(tapes[3].t_shares[2], (uint64_t*)sig->proofs[rep].delta_ts_bytes[2],
               tapes[3].t_shares[2]);
        GF_add(tapes[3].c_share, (uint64_t*)sig->proofs[rep].delta_c_bytes,
               tapes[3].c_share);
      }
      for (unsigned j = 0; j < 4; j++)
      {
        GF_copy(tapes[j].pt_share, mult_chk.pt_share[party + j]);
        GF_copy(tapes[j].t_shares[0], mult_chk.x_shares[0][party + j]);
        GF_copy(tapes[j].t_shares[1], mult_chk.x_shares[1][party + j]);
        GF_copy(tapes[j].t_shares[2], mult_chk.x_shares[2][party + j]);
        GF_copy(tapes[j].a_share, alpha_v_shares[0][party + j]);
        GF_copy(tapes[j].c_share, alpha_v_shares[1][party + j]);
      }
    }

    aim2_mpc_N(matrix_A, ct, &mult_chk);
    hash_update(&ctx_h1, sig->proofs[rep].delta_pt_bytes,
                AIM2_NUM_BYTES_FIELD * (L + 2));

    POLY_mul_add_N(mult_chk.x_shares[0], epsilons[0],
                   alpha_v_shares[0], alpha_v_shares_hi[0]);
    POLY_mul_add_N(mult_chk.z_shares[0], epsilons[0],
                   alpha_v_shares[1], alpha_v_shares_hi[1]);

    POLY_mul_add_N(mult_chk.x_shares[1], epsilons[1],
                   alpha_v_shares[0], alpha_v_shares_hi[0]);
    POLY_mul_add_N(mult_chk.z_shares[1], epsilons[1],
                   alpha_v_shares[1], alpha_v_shares_hi[1]);

    POLY_mul_add_N(mult_chk.x_shares[2], epsilons[2],
                   alpha_v_shares[0], alpha_v_shares_hi[0]);
    POLY_mul_add_N(mult_chk.z_shares[2], epsilons[2],
                   alpha_v_shares[1], alpha_v_shares_hi[1]);

    POLY_mul_add_N(mult_chk.x_shares[3], epsilons[3],
                   alpha_v_shares[0], alpha_v_shares_hi[0]);
    POLY_mul_add_N(mult_chk.z_shares[3], epsilons[3],
                   alpha_v_shares[1], alpha_v_shares_hi[1]);

    POLY_red_N(alpha_v_shares_hi[0], alpha_v_shares[0]);
    POLY_red_N(alpha_v_shares_hi[1], alpha_v_shares[1]);

    // open alpha
    GF_from_bytes(sig->proofs[rep].missing_alpha_share_bytes, alpha_v_shares[0][i_bar]);
    GF_set0(alpha_v);
    for (size_t party = 0; party < N; party++)
    {
      GF_add(alpha_v, alpha_v_shares[0][party], alpha_v);
    }

    // alpha is opened, so we can finish calculating v_share
    GF_mul_add_N(mult_chk.pt_share, alpha_v, alpha_v_shares[1]);
    GF_set0(alpha_v);
    for (size_t party = 0; party < N; party++)
    {
      GF_add(alpha_v, alpha_v_shares[1][party], alpha_v);
    }
    GF_add(alpha_v, alpha_v_shares[1][i_bar], alpha_v_shares[1][i_bar]);
    hash_update(&ctx_h2, (const uint8_t*)alpha_v_shares,
                AIM2_NUM_BYTES_FIELD * 2 * N);
  }

  uint8_t h_1_prime[AIMER_COMMIT_SIZE];
  hash_final(&ctx_h1);
  hash_squeeze(&ctx_h1, h_1_prime, AIMER_COMMIT_SIZE);

  uint8_t h_2_prime[AIMER_COMMIT_SIZE];
  hash_final(&ctx_h2);
  hash_squeeze(&ctx_h2, h_2_prime, AIMER_COMMIT_SIZE);

  if (memcmp(h_1_prime, sig->h_1, AIMER_COMMIT_SIZE) != 0
      || memcmp(h_2_prime, sig->h_2, AIMER_COMMIT_SIZE) != 0)
  {
    ret = -1;
  }

  return ret;
}

#undef N
#undef TAU
#undef L
