// SPDX-License-Identifier: MIT

#include <stdlib.h>
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

void commit_to_seed_and_expand_tape(const uint8_t seed[AIMER_SEED_SIZE],
                                    const uint8_t salt[AIMER_SALT_SIZE],
                                    size_t rep, size_t party,
                                    uint8_t commit[AIMER_COMMIT_SIZE],
                                    tape_t* tape)
{
  hash_instance ctx;

  hash_init_prefix(&ctx, HASH_PREFIX_5);
  hash_update(&ctx, salt, AIMER_SALT_SIZE);
  hash_update(&ctx, (const uint8_t*)&rep, sizeof(uint8_t));
  hash_update(&ctx, (const uint8_t*)&party, sizeof(uint8_t));
  hash_update(&ctx, seed, AIMER_SEED_SIZE);
  hash_final(&ctx);

  hash_squeeze(&ctx, commit, AIMER_COMMIT_SIZE);
  hash_squeeze_GF(&ctx, tape->pt_share);
  for (size_t i = 0; i < AIM2_NUM_INPUT_SBOX; i++)
  {
    hash_squeeze_GF(&ctx, tape->t_shares[i]);
  }
  hash_squeeze_GF(&ctx, tape->a_share);
  hash_squeeze_GF(&ctx, tape->c_share);
}

void aim2_mpc(const GF matrix_A[L][AIM2_NUM_BITS_FIELD], const GF vector_b,
              const GF ct, const tape_t* tape, mult_chk_t* mult_chk)
{
  GF_copy(vector_b, mult_chk->x_shares[L]);
  for (size_t ell = 0; ell < L; ell++)
  {
    // pt + c = t ^ {2 ^ e - 1}
    // --> t ^ {2 ^ e} + t * c = t * pt
    // --> z = x * pt
    GF_copy(tape->t_shares[ell], mult_chk->x_shares[ell]);
    GF_exp_power_of_2(tape->t_shares[ell], mult_chk->z_shares[ell],
                      aim2_invmer_exponents[ell]);
    GF_mul_add(tape->t_shares[ell], aim2_constants[ell],
               mult_chk->z_shares[ell]);

    // x_* = sum_i A[i] * t[i] + b
    GF_transposed_matmul_add(tape->t_shares[ell], matrix_A[ell],
                             mult_chk->x_shares[L]);
  }
  // x ^ {2 ^ e - 1} = pt + ct
  // --> x ^ {2 ^ e} + x * ct = x * pt
  // --> z = x * pt
  GF_exp_power_of_2(mult_chk->x_shares[L], mult_chk->z_shares[L],
                    aim2_invmer_exponents[L]);
  GF_mul_add(ct, mult_chk->x_shares[L], mult_chk->z_shares[L]);
}

// committing to the seeds and the execution views of the parties
void run_phase_1(const uint8_t* public_key, const uint8_t* private_key,
                 const uint8_t* message, const size_t message_len,
                 signature_t* sig, uint8_t commits[TAU][N][AIMER_COMMIT_SIZE],
                 uint8_t nodes[TAU][2 * N][AIMER_SEED_SIZE],
                 mult_chk_t mult_chk[TAU][N])
{
  GF pt = {0,}, ct = {0,};
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
  GF sbox_outputs[AIM2_NUM_INPUT_SBOX];
  aim2_sbox_outputs(pt, sbox_outputs);

  // derive the binary matrix and the vector from the initial vector
  GF matrix_A[AIM2_NUM_INPUT_SBOX][AIM2_NUM_BITS_FIELD];
  GF vector_b = {0,};
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

  // generate root seeds of seed trees
  uint8_t root_seeds[TAU][AIMER_SEED_SIZE];
  hash_squeeze(&ctx, root_seeds[0], sizeof(root_seeds));

  // hash_instance for h_1
  hash_init_prefix(&ctx, HASH_PREFIX_1);
  hash_update(&ctx, mu, AIMER_COMMIT_SIZE);
  hash_update(&ctx, sig->salt, AIMER_SALT_SIZE);

  for (size_t rep = 0; rep < TAU; rep++)
  {
    // compute parties' seeds using binary tree
    expand_tree(root_seeds[rep], sig->salt, rep, nodes[rep]);

    // initialize adjustment values
    tape_t delta;
    memset(&delta, 0, sizeof(tape_t));

    GF lin_const = {0,};
    for (size_t party = 0; party < N; party++)
    {
      // generate execution views and commitments
      tape_t tape;
      commit_to_seed_and_expand_tape(nodes[rep][N + party], sig->salt, rep,
                                     party, commits[rep][party], &tape);
      hash_update(&ctx, commits[rep][party], AIMER_COMMIT_SIZE);

      // compute offsets
      GF_add(delta.pt_share, tape.pt_share, delta.pt_share);
      for (size_t i = 0; i < L; i++)
      {
        GF_add(delta.t_shares[i], tape.t_shares[i], delta.t_shares[i]);
      }
      GF_add(delta.a_share, tape.a_share, delta.a_share);
      GF_add(delta.c_share, tape.c_share, delta.c_share);

      // adjust the last share and prepare the proof and h_1
      if (party == N - 1)
      {
        GF_add(delta.pt_share, pt, delta.pt_share);
        for (size_t i = 0; i < L; i++)
        {
          GF_add(delta.t_shares[i], sbox_outputs[i], delta.t_shares[i]);
        }
        GF_mul_add(pt, delta.a_share, delta.c_share);

        GF_to_bytes(delta.pt_share, sig->proofs[rep].delta_pt_bytes);
        for (size_t i = 0; i < L; i++)
        {
          GF_to_bytes(delta.t_shares[i], sig->proofs[rep].delta_ts_bytes[i]);
        }
        GF_to_bytes(delta.c_share, sig->proofs[rep].delta_c_bytes);

        GF_add(delta.pt_share, tape.pt_share, tape.pt_share);
        for (size_t i = 0; i < L; i++)
        {
          GF_add(delta.t_shares[i], tape.t_shares[i], tape.t_shares[i]);
        }
        GF_add(delta.c_share, tape.c_share, tape.c_share);

        GF_copy(vector_b, lin_const);
      }

      // run the MPC simulation and prepare the mult check inputs
      aim2_mpc(matrix_A, lin_const, ct, &tape, &mult_chk[rep][party]);
      GF_copy(tape.pt_share, mult_chk[rep][party].pt_share);
      GF_copy(tape.a_share, mult_chk[rep][party].a_share);
      GF_copy(tape.c_share, mult_chk[rep][party].c_share);
    }

    // NOTE: depend on the order of values in proof_t
    hash_update(&ctx, sig->proofs[rep].delta_pt_bytes,
                AIM2_NUM_BYTES_FIELD * (L + 2));
  }

  // commit to salt, (all commitments of parties' seeds,
  // delta_pt, delta_t, delta_c) for all repetitions
  hash_final(&ctx);
  hash_squeeze(&ctx, sig->h_1, AIMER_COMMIT_SIZE);
}

void run_phase_3(const uint8_t h_1[AIMER_COMMIT_SIZE],
                 const uint8_t salt[AIMER_SALT_SIZE],
                 const GF epsilons[TAU][L + 1],
                 const mult_chk_t mult_chk[TAU][N],
                 GF alpha_v_shares[TAU][N][2],
                 uint8_t h_2[AIMER_COMMIT_SIZE])
{
  hash_instance ctx;
  hash_init_prefix(&ctx, HASH_PREFIX_2);
  hash_update(&ctx, h_1, AIMER_COMMIT_SIZE);
  hash_update(&ctx, salt, AIMER_SALT_SIZE);

  for (size_t rep = 0; rep < TAU; rep++)
  {
    GF alpha = {0,};
    for (size_t party = 0; party < N; party++)
    {
      // alpha_share = a_share + sum x_share[i] * eps[i]
      // v_share = c_share - pt_share * alpha + sum z_share[i] * eps[i]
      GF_copy(mult_chk[rep][party].a_share, alpha_v_shares[rep][party][0]);
      GF_copy(mult_chk[rep][party].c_share, alpha_v_shares[rep][party][1]);
      for (size_t ell = 0; ell < L + 1; ell++)
      {
        GF_mul_add(mult_chk[rep][party].x_shares[ell], epsilons[rep][ell],
                   alpha_v_shares[rep][party][0]);
        GF_mul_add(mult_chk[rep][party].z_shares[ell], epsilons[rep][ell],
                   alpha_v_shares[rep][party][1]);
      }
      GF_add(alpha, alpha_v_shares[rep][party][0], alpha);
      hash_update_GF(&ctx, alpha_v_shares[rep][party][0]);
    }

    // alpha is opened, so we can finish calculating v_share
    for (size_t party = 0; party < N; party++)
    {
      GF_mul_add(mult_chk[rep][party].pt_share, alpha,
                 alpha_v_shares[rep][party][1]);
      hash_update_GF(&ctx, alpha_v_shares[rep][party][1]);
    }
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
  uint8_t (*nodes)[2 * N][AIMER_SEED_SIZE] = malloc(sizeof(*nodes) * TAU);
  memset(nodes, 0, TAU * sizeof(*nodes));

  // commitments for seeds
  uint8_t (*commits)[N][AIMER_COMMIT_SIZE] = malloc(sizeof(*commits) * TAU);
  memset(commits, 0, TAU * sizeof(*commits));

  // multiplication check inputs
  mult_chk_t (*mult_chk)[N] = malloc(sizeof(*mult_chk) * TAU);
  memset(mult_chk, 0, TAU * sizeof(*mult_chk));

  // commitments for phase 1
  run_phase_1(public_key, private_key, message, message_len,
              sig, commits, nodes, mult_chk);

  ////////////////////////////////////////////////////////////////
  // Phase 2: Challenging the multiplication checking protocol. //
  ////////////////////////////////////////////////////////////////

  GF (*epsilons)[L + 1] = malloc(sizeof(*epsilons) * TAU);
  memset(epsilons, 0, TAU * sizeof(*epsilons));

  hash_init(&ctx);
  hash_update(&ctx, sig->h_1, AIMER_COMMIT_SIZE);
  hash_final(&ctx);
  for (size_t rep = 0; rep < TAU; rep++)
  {
    for (size_t ell = 0; ell < L + 1; ell++)
    {
      hash_squeeze_GF(&ctx, epsilons[rep][ell]);
    }
  }

  /////////////////////////////////////////////////////////////////
  // Phase 3: Committing to the simulation of the multiplication //
  //          checking protocol.                                 //
  /////////////////////////////////////////////////////////////////

  // multiplication check outputs
  GF (*alpha_v_shares)[N][2] = malloc(sizeof(*alpha_v_shares) * TAU);
  memset(alpha_v_shares, 0, TAU * sizeof(*alpha_v_shares));

  // compute the commitment of phase 3
  run_phase_3(sig->h_1, sig->salt, epsilons, mult_chk, alpha_v_shares, sig->h_2);

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
    GF_to_bytes(alpha_v_shares[rep][i_bar][0],
                sig->proofs[rep].missing_alpha_share_bytes);
  }

  free(nodes);
  free(commits);
  free(mult_chk);
  free(alpha_v_shares);
  free(epsilons);

  return ret;
}

int aimer_verify(const uint8_t* public_key, const uint8_t* signature,
                 const uint8_t* message, const size_t message_len)
{
  int ret = 0;
  const signature_t* sig = (const signature_t*)signature;

  uint8_t iv[AIM2_IV_SIZE];
  memcpy(iv, public_key, AIM2_IV_SIZE);

  GF ct = {0,};
  GF_from_bytes(public_key + AIM2_IV_SIZE, ct);

  // derive the binary matrix and the vector from the initial vector
  GF matrix_A[AIM2_NUM_INPUT_SBOX][AIM2_NUM_BITS_FIELD];
  GF vector_b = {0,};
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

  for (size_t rep = 0; rep < TAU; rep++)
  {
    uint8_t nodes[2 * AIMER_NUM_MPC_PARTIES][AIMER_SEED_SIZE] = {0,};
    size_t i_bar = indices[rep];
    reconstruct_seed_tree(sig->proofs[rep].reveal_path, i_bar,
                          sig->salt, rep, nodes);

    GF lin_const = {0,}, alpha = {0,}, alpha_share = {0,};
    GF v_shares[N], pt_shares[N];
    GF_set0(v_shares[i_bar]);

    GF epsilons[L + 1];
    for (size_t ell = 0; ell < L + 1; ell++)
    {
      hash_squeeze_GF(&ctx_e, epsilons[ell]);
    }

    for (size_t party = 0; party < N; party++)
    {
      if (party == i_bar)
      {
        hash_update(&ctx_h1, sig->proofs[rep].missing_commitment,
                    AIMER_COMMIT_SIZE);

        hash_update(&ctx_h2, sig->proofs[rep].missing_alpha_share_bytes,
                    AIM2_NUM_BYTES_FIELD);
        GF_from_bytes(sig->proofs[rep].missing_alpha_share_bytes, alpha_share);
        GF_add(alpha, alpha_share, alpha);
        continue;
      }

      tape_t tape;
      uint8_t commit[AIMER_COMMIT_SIZE];
      commit_to_seed_and_expand_tape(nodes[N + party], sig->salt, rep,
                                     party, commit, &tape);
      hash_update(&ctx_h1, commit, AIMER_COMMIT_SIZE);

      // adjust last shares
      if (party == N - 1)
      {
        GF_copy(vector_b, lin_const);

        GF delta;
        GF_from_bytes(sig->proofs[rep].delta_pt_bytes, delta);
        GF_add(tape.pt_share, delta, tape.pt_share);

        GF_from_bytes(sig->proofs[rep].delta_c_bytes, delta);
        GF_add(tape.c_share, delta, tape.c_share);

        for (size_t ell = 0; ell < L; ell++)
        {
          GF_from_bytes(sig->proofs[rep].delta_ts_bytes[ell], delta);
          GF_add(tape.t_shares[ell], delta, tape.t_shares[ell]);
        }
      }

      // run the MPC simulation and prepare the mult check inputs
      mult_chk_t mult_chk;
      memset(&mult_chk, 0, sizeof(mult_chk_t));

      aim2_mpc(matrix_A, lin_const, ct, &tape, &mult_chk);
      GF_copy(tape.pt_share, mult_chk.pt_share);
      GF_copy(tape.a_share, mult_chk.a_share);
      GF_copy(tape.c_share, mult_chk.c_share);

      // mult check
      GF_copy(mult_chk.pt_share, pt_shares[party]);
      GF_copy(mult_chk.a_share, alpha_share);
      GF_copy(mult_chk.c_share, v_shares[party]);
      for (size_t ell = 0; ell < L + 1; ell++)
      {
        GF_mul_add(mult_chk.x_shares[ell], epsilons[ell], alpha_share);
        GF_mul_add(mult_chk.z_shares[ell], epsilons[ell], v_shares[party]);
      }
      GF_add(alpha, alpha_share, alpha);
      hash_update_GF(&ctx_h2, alpha_share);
    }

    // alpha is opened, so we can finish calculating v_share
    for (size_t party = 0; party < N; party++)
    {
      if (party == i_bar)
      {
        continue;
      }

      GF_mul_add(pt_shares[party], alpha, v_shares[party]);
      GF_add(v_shares[i_bar], v_shares[party], v_shares[i_bar]);
    }

    // v is opened
    for (size_t party = 0; party < N; party++)
    {
      hash_update_GF(&ctx_h2, v_shares[party]);
    }

    // NOTE: depend on the order of values in proof_t
    hash_update(&ctx_h1, sig->proofs[rep].delta_pt_bytes,
                AIM2_NUM_BYTES_FIELD * (L + 2));
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
