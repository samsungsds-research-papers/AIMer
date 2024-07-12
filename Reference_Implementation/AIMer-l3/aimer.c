// -----------------------------------------------------------------------------
// File Name   : aimer.c
// Description : 
// SPDX-License-Identifier: MIT
// -----------------------------------------------------------------------------

#include "aim.h"
#include "aimer.h"
#include "aimer_internal.h"
#include "rng.h"
#include "tree.h"

#include <stdlib.h>

int _aimer_sign();
int _aimer_verify();
int serialize_signature();
int deserialize_signature();

int aimer_keygen(aimer_params_t param, aimer_publickey_t* public_key,
                 aimer_privatekey_t* private_key)
{
  if (!public_key || !private_key)
  {
    return -1;
  }

  const aimer_instance_t* instance = aimer_instance_get(param);
  if (!instance)
  {
    return -1;
  }

  randombytes(public_key->data, instance->aim_params.block_size);
  randombytes(private_key->data, instance->aim_params.block_size);

  aim(private_key->data, public_key->data,
      &public_key->data[instance->aim_params.block_size]);

  public_key->params = param;
  return 0;
}

int _aimer_sign(const aimer_instance_t*   instance,
                const aimer_publickey_t*  public_key,
                const aimer_privatekey_t* private_key,
                const uint8_t* message, const size_t message_len,
                signature_t* sig)
{
  const size_t L                = instance->aim_params.num_input_sboxes + 1;
  const size_t block_size       = instance->aim_params.block_size;
  const size_t random_tape_size = instance->aim_params.block_size +
                                  (L - 1) * instance->field_size +
                                  instance->field_size +
                                  instance->field_size;
  const size_t tau              = instance->num_repetitions;
  const size_t N                = instance->num_MPC_parties;

  int ret = 0;
  GF *sbox_pairs;
  GF input_GF = {0,};
  GF vector_B = {0,};
  GF **matrix_A = malloc((L - 1) * sizeof(GF*));

  uint8_t *iv, *output, *output2;
  uint8_t *master_seeds, *party_seed_commitments;

  tree_t        **seed_trees;
  random_tape_t *random_tapes;

  iv  = malloc(block_size);
  output  = malloc(block_size);
  output2 = malloc(block_size);
  memcpy(iv, public_key->data, block_size);
  memcpy(output, public_key->data + block_size, block_size);
  GF_from_bytes(private_key->data, input_GF);

  generate_linear_layer(iv, matrix_A, vector_B);

  sbox_pairs =
    (GF*)malloc(instance->field_size * L);

  aim_with_sbox_output((uint8_t*)private_key->data, (const GF**)matrix_A,
                       vector_B, output2, sbox_pairs);
  if (memcmp(output, output2, block_size))
  {
    for (size_t i = 0; i < L - 1; i++)
    {
      free(matrix_A[i]);
    }
    free(matrix_A);

    free(sbox_pairs);
    free(iv);
    free(output);
    free(output2);

    return -1;
  }

  master_seeds = malloc(instance->seed_size * tau);
  generate_salt_and_seeds(instance, public_key, private_key, message,
                          message_len, sig->salt, master_seeds);

  party_seed_commitments = malloc(tau * N * instance->digest_size);

  random_tapes                   = malloc(sizeof(random_tape_t));
  random_tapes->tape             = malloc(tau * N * random_tape_size);
  random_tapes->random_tape_size = random_tape_size;

  seed_trees = (tree_t**)malloc(sizeof(tree_t*) * tau);
  for (size_t repetition = 0; repetition < tau; repetition++)
  {
    seed_trees[repetition] =
      seed_tree(&master_seeds[repetition * instance->seed_size],
                instance->seed_size, sig->salt, instance->salt_size, N,
                repetition);

    for (size_t party = 0; party < N; party++)
    {
      commit_to_party_seed(instance, get_leaf(seed_trees[repetition], party),
                           sig->salt, repetition, party,
                           party_seed_commitments +
                           (repetition * N + party) * instance->digest_size);

      generate_tape(instance, get_leaf(seed_trees[repetition], party),
                    sig->salt, repetition, party, random_tapes);
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Phase 1: commit to MPC executions and polynomials
  //////////////////////////////////////////////////////////////////////////////
  GF* repetition_shared_inputs  =
    malloc(tau * N     * instance->aim_params.block_size);

  GF* repetition_shared_x     = malloc(tau * N * L * instance->field_size);
  GF* repetition_shared_z     = malloc(tau * N * L * instance->field_size);

  GF* repetition_shared_dot_a = malloc(tau * N     * instance->field_size);
  GF* repetition_shared_dot_c = malloc(tau * N     * instance->field_size);

  GF  input_delta = {0,};
  for (size_t repetition = 0; repetition < tau; repetition++)
  {
    proof_t* proof = &sig->proofs[repetition];
    allocate_proof(instance, proof);

    // Generate sharing of secret key
    GF_copy(input_GF, input_delta);

    for (size_t party = 0; party < N; party++)
    {
      GF* shared_input =
        repetition_shared_inputs + (repetition * N + party);

      uint8_t* random_share =
        random_tapes->tape + (repetition * N + party) * random_tape_size;

      GF_from_bytes(random_share, shared_input[0]);
      GF_add(input_delta, shared_input[0], input_delta);
    }

    // Fix first share
    GF* first_share_input = repetition_shared_inputs + (repetition * N);

    GF_add(first_share_input[0], input_delta, first_share_input[0]);
    GF_to_bytes(input_delta, proof->input_delta);

    for (size_t ell = 0; ell < (L - 1); ell++)
    {
      GF_copy(sbox_pairs[ell], proof->z_delta[ell]);
    }

    for (size_t party = 0; party < N; party++)
    {
      GF* shared_z =
        repetition_shared_z + (repetition * N + party) * L;

      uint8_t* random_z_shares =
        random_tapes->tape + (repetition * N + party) * random_tape_size +
        instance->aim_params.block_size;

      for (size_t ell = 0; ell < (L - 1); ell++)
      {
        GF_from_bytes(random_z_shares + ell * instance->aim_params.block_size,
                      shared_z[ell]);
        GF_add(proof->z_delta[ell], shared_z[ell], proof->z_delta[ell]);
      }
    }

    GF* first_shared_z = repetition_shared_z + (repetition * N * L);

    for (size_t ell = 0; ell < (L - 1); ell++)
    {
      GF_add(first_shared_z[ell], proof->z_delta[ell], first_shared_z[ell]);
    }

    // Get shares of sbox inputs by simulating MPC execution
    GF* shared_x_offset = repetition_shared_x + (repetition * N * L);
    GF* shared_z_offset = repetition_shared_z + (repetition * N * L);

    uint8_t* shared_inputs_offset =
      (uint8_t*)(repetition_shared_inputs + (repetition * N));

    aim_mpc(shared_inputs_offset, (const GF**)matrix_A, vector_B, output,
            N, shared_z_offset, shared_x_offset);
  }

  GF a = {0,};
  GF c = {0,};
  for (size_t repetition = 0; repetition < tau; repetition++)
  {
    GF_set0(a);
    GF_set0(c);

    for (size_t party = 0; party < N; party++)
    {
      uint8_t* random_triple = random_tapes->tape +
        ((repetition * N + party) * random_tape_size +
          instance->aim_params.block_size + (L - 1) * instance->field_size);

      GF* a_share = repetition_shared_dot_a + (repetition * N + party);

      GF_from_bytes(random_triple, a_share[0]);
      GF_add(a, a_share[0], a);

      GF* c_share = repetition_shared_dot_c + repetition * N + party;

      GF_from_bytes(random_triple + instance->field_size, c_share[0]);
      GF_add(c, c_share[0], c);
    }

    // Calculate c_delta = -c + a*b, c is negated above already
    // Calculate c_delta that fixes the dot triple
    GF_mul(a, input_GF, a);
    GF_add(a, c, sig->proofs[repetition].c_delta);

    // Fix party 0's share
    GF_add(repetition_shared_dot_c[repetition * N],
           sig->proofs[repetition].c_delta,
           repetition_shared_dot_c[repetition * N]);
  }

  //////////////////////////////////////////////////////////////////////////////
  // Phase 2: challenge the checking polynomials
  //////////////////////////////////////////////////////////////////////////////
  // Commit to salt, (all commitments of parties seeds,
  // input_delta, z_delta, c_delta) for all repetitions
  GF** epsilons = malloc(tau * sizeof(GF*));

  phase_1_commitment(instance, sig, public_key, message, message_len,
                     party_seed_commitments, sig->h_1);

  // Expand challenge hash to epsilon values
  phase_1_expand(instance, sig->h_1, epsilons);

  //////////////////////////////////////////////////////////////////////////////
  // Phase 3: commit to the views of the checking protocol
  //////////////////////////////////////////////////////////////////////////////
  GF*  repetition_alpha_shares =
         calloc(tau * N * instance->field_size, sizeof(uint8_t));
  GF** v_shares = malloc(tau * sizeof(GF*));

  GF alpha = {0,};
  for (size_t repetition = 0; repetition < tau; repetition++)
  {
    GF* x_shares     = repetition_shared_x     + (repetition * N) * L;
    GF* z_shares     = repetition_shared_z     + (repetition * N) * L;
    GF* dot_a_share  = repetition_shared_dot_a + (repetition * N);
    GF* dot_c_share  = repetition_shared_dot_c + (repetition * N);
    GF* alpha_share  = repetition_alpha_shares + (repetition * N);

    // Execute sacrificing check protocol
    // alpha = eps_i * x_i + a
    GF_set0(alpha);
    for (size_t party = 0; party < N; party++)
    {
      GF temp = {0,};

      GF_copy(dot_a_share[party], alpha_share[party]);
      for (size_t ell = 0; ell < L; ell++)
      {
        GF_mul(x_shares[party * L + ell], epsilons[repetition][ell], temp);
        GF_add(alpha_share[party], temp, alpha_share[party]);
      }
      GF_add(alpha, alpha_share[party], alpha);
    }

    // v^i = - c^i + dot(alpha, y^i) - dot(eps, z^i)
    v_shares[repetition] = malloc(N * instance->field_size);
    for (size_t party = 0; party < N; party++)
    {
      GF_copy(dot_c_share[party], v_shares[repetition][party]);

      GF temp = {0,};
      GF* input_share = repetition_shared_inputs + (repetition * N + party);

      GF_mul(alpha, input_share[0], temp);
      GF_add(v_shares[repetition][party], temp, v_shares[repetition][party]);
      for (size_t ell = 0; ell < L; ell++)
      {
        GF_mul(epsilons[repetition][ell], z_shares[party * L + ell], temp);
        GF_add(v_shares[repetition][party], temp, v_shares[repetition][party]);
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Phase 4: challenge the views of the checking protocol
  //////////////////////////////////////////////////////////////////////////////
  uint16_t* missing_parties = malloc(tau * sizeof(uint16_t));

  phase_2_commitment(instance, sig->salt, sig->h_1, repetition_alpha_shares,
                     (const GF**)v_shares, sig->h_2);

  phase_2_expand(instance, sig->h_2, missing_parties);

  //////////////////////////////////////////////////////////////////////////////
  // Phase 5: Open the views of the checking protocol
  //////////////////////////////////////////////////////////////////////////////
  // Build signature
  for (size_t repetition = 0; repetition < tau; repetition++)
  {
    proof_t* proof = &sig->proofs[repetition];

    size_t missing_party = missing_parties[repetition];
    proof->reveal_list = reveal_all_but(seed_trees[repetition], missing_party);

    memcpy(proof->missing_commitment,
           party_seed_commitments +
           ((repetition * N + missing_party) * instance->digest_size),
           instance->digest_size);

    GF_copy(repetition_alpha_shares[repetition * N + missing_party],
            proof->missing_alpha_share);
  }

  //////////////////////////////////////////////////////////////////////////////
  // Free memory allocations
  //////////////////////////////////////////////////////////////////////////////
  for (size_t i = 0; i < tau; i++)
  {
    free_tree(seed_trees[i]);

    free(epsilons[i]);
    free(v_shares[i]);
  }
  for (size_t i = 0; i < L - 1; i++)
  {
    free(matrix_A[i]);
  }
  free(matrix_A);  

  free(iv);
  free(output);
  free(output2);
  free(sbox_pairs);
  free(master_seeds);
  free(seed_trees);
  free(party_seed_commitments);
  free(random_tapes->tape);
  free(random_tapes);
  free(repetition_shared_inputs);
  free(repetition_shared_x);
  free(repetition_shared_z);
  free(repetition_shared_dot_a);
  free(repetition_shared_dot_c);
  free(epsilons);
  free(v_shares);
  free(missing_parties);
  free(repetition_alpha_shares);

  return ret;
}

int serialize_signature(const aimer_instance_t* instance,
                        const signature_t* sig,
                        uint8_t* signature, size_t* signature_len)
{
  const size_t L   = instance->aim_params.num_input_sboxes;
  const size_t tau = instance->num_repetitions;
  uint8_t* signature_offset = signature;

  memcpy(signature_offset, sig->salt, instance->salt_size);
  signature_offset += instance->salt_size;

  memcpy(signature_offset, sig->h_1, instance->digest_size);
  signature_offset += instance->digest_size;

  memcpy(signature_offset, sig->h_2, instance->digest_size);
  signature_offset += instance->digest_size;

  size_t reveal_list_size =
    ceil_log2(instance->num_MPC_parties) * instance->seed_size;

  for (size_t repetition = 0; repetition < tau; repetition++)
  {
    proof_t* proof = &sig->proofs[repetition];

    memcpy(signature_offset, proof->reveal_list.data, reveal_list_size);
    signature_offset += reveal_list_size;

    memcpy(signature_offset, proof->missing_commitment, instance->digest_size);
    signature_offset += instance->digest_size;

    memcpy(signature_offset, proof->input_delta,
           instance->aim_params.block_size);
    signature_offset += instance->aim_params.block_size;

    for (size_t ell = 0; ell < L; ell++)
    {
      GF_to_bytes(proof->z_delta[ell], signature_offset);
      signature_offset += instance->field_size;
    }

    GF_to_bytes(proof->c_delta, signature_offset);
    signature_offset += instance->field_size;

    GF_to_bytes(proof->missing_alpha_share, signature_offset);
    signature_offset += instance->field_size;
  }

  *signature_len = (int)(signature_offset - signature);
  if (*signature_len > AIMER_MAX_SIGNATURE_SIZE)
  {
    return -1;
  }
  return 0;
}

int aimer_sign(const aimer_publickey_t*  public_key,
               const aimer_privatekey_t* private_key,
               const uint8_t* message, const size_t message_len,
               uint8_t* signature, size_t* signature_len)
{
  int ret = 0;
  signature_t sig;
  aimer_params_t params = public_key->params;

  const aimer_instance_t *instance = aimer_instance_get(params);
  if (!instance)
  {
    return -1;
  }

  allocate_signature(instance, &sig);
  ret =
    _aimer_sign(instance, public_key, private_key, message, message_len, &sig);
  if (ret == -1)
  {
    free_signature(instance, &sig);
    return -1;
  }

  ret = serialize_signature(instance, &sig, signature, signature_len);
  if (ret == -1)
  {
    free_signature(instance, &sig);
    return -1;
  }

  free_signature(instance, &sig);
  return ret;
}

int deserialize_signature(const aimer_instance_t* instance,
                          const uint8_t* signature, const size_t signature_len,
                          signature_t* sig, uint16_t* missing_parties)
{
  const size_t L   = instance->aim_params.num_input_sboxes;
  const size_t tau = instance->num_repetitions;

  size_t offset = 0;

  memcpy(sig->salt, signature, instance->salt_size);
  offset += instance->salt_size;

  memcpy(sig->h_1, signature + offset, instance->digest_size);
  offset += instance->digest_size;
  memcpy(sig->h_2, signature + offset, instance->digest_size);
  offset += instance->digest_size;

  phase_2_expand(instance, sig->h_2, missing_parties);

  size_t reveal_list_size = ceil_log2(instance->num_MPC_parties);

  for (size_t repetition = 0; repetition < tau; repetition++)
  {
    proof_t* proof = &sig->proofs[repetition];

    allocate_proof(instance, proof);

    proof->reveal_list.seed_size = instance->seed_size;
    proof->reveal_list.missing_leaf = missing_parties[repetition];
    proof->reveal_list.data = malloc(reveal_list_size * instance->seed_size);

    memcpy(proof->reveal_list.data,
           signature + offset, reveal_list_size * instance->seed_size);
    offset += reveal_list_size * instance->seed_size;

    memcpy(proof->missing_commitment,
           signature + offset, instance->digest_size);
    offset += instance->digest_size;

    memcpy(proof->input_delta,
           signature + offset, instance->aim_params.block_size);
    offset += instance->aim_params.block_size;

    for (size_t ell = 0; ell < L; ell++)
    {
      GF_from_bytes(signature + offset, proof->z_delta[ell]);
      offset += instance->field_size;
    }

    GF_from_bytes(signature + offset, proof->c_delta);
    offset += instance->field_size;

    GF_from_bytes(signature + offset, proof->missing_alpha_share);
    offset += instance->field_size;
  }

  if (offset != signature_len)
  {
    return -1;
  }

  return 0;
}

int _aimer_verify(const aimer_instance_t*  instance,
                  const aimer_publickey_t* public_key,
                  const signature_t* sig,
                  const uint16_t* missing_parties,
                  const uint8_t* message, const size_t message_len)
{
  const size_t L                = instance->aim_params.num_input_sboxes + 1;
  const size_t block_size       = instance->aim_params.block_size;
  const size_t random_tape_size = instance->aim_params.block_size +
                                  (L - 1) * instance->field_size +
                                  instance->field_size +
                                  instance->field_size;
  const size_t tau              = instance->num_repetitions;
  const size_t N                = instance->num_MPC_parties;

  int ret = 0;

  // h_1 expansion,
  // h_2 expansion already happened in deserialize
  GF** epsilons = malloc(tau * sizeof(GF*));

  phase_1_expand(instance, sig->h_1, epsilons);

  // Rebuild seed trees
  tree_t** seed_trees = malloc(tau * sizeof(tree_t*));
  uint8_t* party_seed_commitments = malloc(tau * N * instance->digest_size);

  random_tape_t* random_tapes = malloc(sizeof(random_tape_t));
  random_tapes->tape = malloc(tau * N * random_tape_size);
  random_tapes->random_tape_size = random_tape_size;

  uint8_t* dummy = calloc(instance->seed_size, sizeof(uint8_t));
  for (size_t repetition = 0; repetition < tau; repetition++)
  {
    proof_t* proof = &sig->proofs[repetition];
    seed_trees[repetition] =
      seed_tree_with_reveal_list(&proof->reveal_list, sig->salt,
                                 instance->salt_size, N, repetition);
    for (size_t party = 0; party < N; party++)
    {
      if (party != missing_parties[repetition])
      {
        commit_to_party_seed(instance, get_leaf(seed_trees[repetition], party),
                             sig->salt, repetition, party, 
                             party_seed_commitments +
                             (repetition * N + party) * instance->digest_size);

        generate_tape(instance, get_leaf(seed_trees[repetition], party),
                      sig->salt, repetition, party, random_tapes);
      }
      else
      {
        generate_tape(instance, dummy, sig->salt, repetition, party,
                      random_tapes);
      }
    }
    memcpy
    (party_seed_commitments +
     (repetition * N + missing_parties[repetition]) * instance->digest_size,
     proof->missing_commitment, instance->digest_size);
  }
  free(dummy);

  // Recompute commitments to executions of block cipher
  uint8_t* repetition_shared_inputs =
    malloc(tau * N * instance->aim_params.block_size);
  GF* repetition_shared_x     = malloc(tau * N * L * instance->field_size);
  GF* repetition_shared_z     = malloc(tau * N * L * instance->field_size);
  GF* repetition_shared_dot_a = malloc(tau * N * instance->field_size);
  GF* repetition_shared_dot_c = malloc(tau * N * instance->field_size);

  uint8_t *iv = malloc(block_size);
  uint8_t *output = malloc(block_size);

  memcpy(iv, public_key->data, block_size);
  memcpy(output, public_key->data + block_size, block_size);

  // generate linear layer for AIM mpc
  GF** matrix_A = malloc((L - 1) * sizeof(GF*));
  GF vector_B = {0,};
  generate_linear_layer(iv, matrix_A, vector_B);

  for (size_t repetition = 0; repetition < tau; repetition++)
  {
    proof_t* proof = &sig->proofs[repetition];
    // Generate sharing of secret key
    for (size_t party = 0; party < N; party++)
    {
      uint8_t* shared_input =
        repetition_shared_inputs +
        (repetition * N + party) * instance->aim_params.block_size;
      uint8_t* random_share =
        random_tapes->tape + (repetition * N + party) * random_tape_size;

      memcpy(shared_input, random_share, instance->aim_params.block_size);
    }

    // Fix first share
    uint8_t* first_input_share =
      repetition_shared_inputs +
      (repetition * N * instance->aim_params.block_size);
    for (size_t j = 0; j < instance->aim_params.block_size; j++)
    {
      first_input_share[j] ^= proof->input_delta[j];
    }

    // Generate sharing of z values
    for (size_t party = 0; party < N; party++)
    {
      GF* shared_z = repetition_shared_z + (repetition * N + party) * L;
      uint8_t* random_z_shares =
        random_tapes->tape + (repetition * N + party) * random_tape_size +
        instance->aim_params.block_size;

      for (size_t ell = 0; ell < (L - 1); ell++)
      {
        GF_from_bytes(random_z_shares + ell * instance->aim_params.block_size,
                      shared_z[ell]);
      }
    }

    // Fix first share
    GF* first_shared_z = repetition_shared_z + (repetition * N * L);
    for (size_t ell = 0; ell < (L - 1); ell++)
    {
      GF_add(first_shared_z[ell], proof->z_delta[ell], first_shared_z[ell]);
    }

    // Get shares of sbox inputs by executing MPC AIM
    GF* shared_x_offset = repetition_shared_x + (repetition * N * L);
    GF* shared_z_offset = repetition_shared_z + (repetition * N * L);

    uint8_t* shared_inputs_offset =
      repetition_shared_inputs +
      (repetition * N * instance->aim_params.block_size);

    aim_mpc(shared_inputs_offset, (const GF**)matrix_A, vector_B, output,
            N, shared_z_offset, shared_x_offset);
  }

  // Recompute shares of dot product
  for (size_t repetition = 0; repetition < tau; repetition++)
  {
    proof_t* proof = &sig->proofs[repetition];
    // Also generate valid dot triple a,z,c and save c_delta
    for (size_t party = 0; party < N; party++)
    {
      if (party != missing_parties[repetition])
      {
         uint8_t* random_triple =
           random_tapes->tape +
           ((repetition * N + party) * random_tape_size +
            instance->aim_params.block_size + (L - 1) * instance->field_size);

        GF* a_share = repetition_shared_dot_a + (repetition * N + party);
        GF_from_bytes(random_triple, a_share[0]);

        GF* c_share = repetition_shared_dot_c + (repetition * N + party);
        GF_from_bytes(random_triple + instance->field_size, c_share[0]);
      }
    }

    // Fix party 0's share
    if (0 != missing_parties[repetition])
    {
      // Adjust first share with delta from signature
      GF_add(repetition_shared_dot_c[repetition * N], proof->c_delta,
             repetition_shared_dot_c[repetition * N]);
    }
  }

  // Recompute views of sacrificing checks
  GF* repetition_alpha_shares =
    calloc(tau * N * instance->field_size, sizeof(uint8_t));
  GF** v_shares = malloc(tau * sizeof(GF*));
  GF alpha = {0,};
  for (size_t repetition = 0; repetition < tau; repetition++)
  {
    proof_t* proof = &sig->proofs[repetition];
    size_t missing_party = missing_parties[repetition];

    GF_set0(alpha);
    for (size_t party = 0; party < N; party++)
    {
      GF* alpha_share = repetition_alpha_shares + (repetition * N + party);
      GF* a_share     = repetition_shared_dot_a + (repetition * N + party);
      GF* x_shares    = repetition_shared_x + (repetition * N + party) * L;

      if (party != missing_party)
      {
        GF temp = {0,};

        GF_copy(a_share[0], alpha_share[0]);
        for (size_t ell = 0; ell < L; ell++)
        {
          GF_mul(x_shares[ell], epsilons[repetition][ell], temp);
          GF_add(alpha_share[0], temp, alpha_share[0]);
        }
        GF_add(alpha, alpha_share[0], alpha);
      }
    }

    // Fill missing shares
    GF* missing_alpha =
      repetition_alpha_shares + (repetition * N + missing_party);

    GF_copy(proof->missing_alpha_share, missing_alpha[0]);
    GF_add(alpha, missing_alpha[0], alpha);

    // v^i = - c^i + dot(alpha, y) - dot(eps, z^i)
    v_shares[repetition] = calloc(N * instance->field_size, 1);
    for (size_t party = 0; party < N; party++)
    {
      if (party != missing_party)
      {
        GF temp = {0,};
        GF* z_shares = repetition_shared_z + (repetition * N + party) * L;
        GF* c_share  = repetition_shared_dot_c + (repetition * N + party);
        GF input_share;

        GF_from_bytes(repetition_shared_inputs +
                      (repetition * N + party) * instance->field_size,
                      input_share);

        GF_mul(alpha, input_share, v_shares[repetition][party]);
        GF_add(v_shares[repetition][party], c_share[0],
               v_shares[repetition][party]);
        for (size_t ell = 0; ell < L; ell++)
        {
          GF_mul(epsilons[repetition][ell], z_shares[ell], temp);
          GF_add(v_shares[repetition][party], temp,
                 v_shares[repetition][party]);
        }
      }
    }
    // Calculate missing shares as 0 - sum_{i!=missing} v^i
    for (size_t party = 0; party < N; party++)
    {
      if (party != missing_party)
      {
        GF_add(v_shares[repetition][missing_party], v_shares[repetition][party],
               v_shares[repetition][missing_party]);
      }
    }
  }

  // Recompute h_1 and h_2
  uint8_t* h_1 = malloc(instance->digest_size);
  uint8_t* h_2 = malloc(instance->digest_size);

  phase_1_commitment(instance, sig, public_key, message, message_len,
                     party_seed_commitments, h_1);

  phase_2_commitment(instance, sig->salt, h_1, repetition_alpha_shares,
                     (const GF**)v_shares, h_2);

  if (memcmp(h_1, sig->h_1, instance->digest_size) != 0)
  {
    ret = -1;
  }

  if (memcmp(h_2, sig->h_2, instance->digest_size) != 0)
  {
    ret = -1;
  }

  for (size_t i = 0; i < tau; i++)
  {
    free_tree(seed_trees[i]);
    free(epsilons[i]);
    free(v_shares[i]);
  }
  free(seed_trees);
  free(epsilons);
  free(v_shares);

  for (size_t i = 0; i < L - 1; i++)
  {
    free(matrix_A[i]);
  }
  free(matrix_A);

  free(repetition_alpha_shares);
  free(repetition_shared_inputs);
  free(repetition_shared_x);
  free(repetition_shared_z);
  free(repetition_shared_dot_a);
  free(repetition_shared_dot_c);

  free(party_seed_commitments);
  free(random_tapes->tape);
  free(random_tapes);

  free(h_1);
  free(h_2);
  free(iv);
  free(output);

  return ret;
}

int aimer_verify(const aimer_publickey_t* public_key,
                 const uint8_t* signature, const size_t signature_len,
                 const uint8_t* message, const size_t message_len)
{
  int ret = 0;
  signature_t sig;
  aimer_params_t params = public_key->params;

  const aimer_instance_t *instance = aimer_instance_get(params);
  if (!instance)
  {
    return -1;
  }

  uint16_t* missing_parties =
    malloc(instance->num_repetitions * sizeof(uint16_t));

  allocate_signature(instance, &sig);
  ret = deserialize_signature(instance, signature, signature_len, &sig,
                              missing_parties);
  if (ret == -1)
  {
    free_signature(instance, &sig);
    free(missing_parties);
    return -1;
  }

  ret = _aimer_verify(instance, public_key, &sig, missing_parties, message,
                        message_len);
  if (ret == -1)
  {
    free_signature(instance, &sig);
    free(missing_parties);
    return -1;
  }

  free_signature(instance, &sig);
  free(missing_parties);
  return ret;
}
