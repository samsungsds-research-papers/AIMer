// -----------------------------------------------------------------------------
// File Name   : aimer_internal.c
// Description : 
// SPDX-License-Identifier: MIT
// -----------------------------------------------------------------------------

#include "aimer_internal.h"

#include <string.h>

void allocate_proof(const aimer_instance_t* instance, proof_t* proof)
{
  proof->missing_commitment =
    malloc(instance->digest_size);
  proof->input_delta =
    malloc(instance->aim_params.block_size);
  proof->z_delta =
    malloc((instance->aim_params.num_input_sboxes) * instance->field_size);
}

void allocate_signature(const aimer_instance_t* instance, signature_t* sig)
{
  sig->salt   = malloc(instance->salt_size);
  sig->h_1    = malloc(instance->digest_size);
  sig->h_2    = malloc(instance->digest_size);
  sig->proofs = malloc(instance->num_repetitions * sizeof(proof_t));
}

void free_proof(proof_t* proof)
{
  free(proof->reveal_list.data);
  free(proof->missing_commitment);
  free(proof->input_delta);
  free(proof->z_delta);
}

void free_signature(const aimer_instance_t* instance, signature_t* sig)
{
  for (size_t i = 0; i < instance->num_repetitions; i++)
  {
    free_proof(&sig->proofs[i]);
  }
  free(sig->proofs);
  free(sig->h_1);
  free(sig->h_2);
  free(sig->salt);
}

void generate_tape(const aimer_instance_t* instance, const uint8_t* seed,
                   const uint8_t* salt, size_t repetition, size_t party,
                   random_tape_t* tapes)
{
    hash_instance ctx;

    hash_init(&ctx, instance->digest_size);
    hash_update(&ctx, seed, instance->seed_size);
    hash_update(&ctx, salt, instance->salt_size);
    hash_update_uint16(&ctx, (uint16_t)repetition);
    hash_update_uint16(&ctx, (uint16_t)party);
    hash_final(&ctx);

    hash_squeeze(&ctx,
                tapes->tape +
                (repetition * instance->num_MPC_parties + party) *
                tapes->random_tape_size,
                tapes->random_tape_size);
}

void generate_salt_and_seeds(const aimer_instance_t* instance,
                             const aimer_publickey_t* public_key,
                             const aimer_privatekey_t* private_key,
                             const uint8_t* message, const size_t message_len,
                             uint8_t* salt, uint8_t* seeds)
{
  size_t repetition, seed_size;
  uint8_t buffer[instance->seed_size];
  hash_instance ctx;

  memset(buffer, 0, sizeof(buffer));

  hash_init(&ctx, instance->digest_size);
  hash_update_uint16(&ctx, (uint16_t)instance->params);
  hash_update(&ctx, private_key->data, sizeof(private_key->data));
  hash_update(&ctx, public_key->data, sizeof(public_key->data));
  hash_update(&ctx, message, message_len);
  hash_final(&ctx);

  hash_squeeze(&ctx, salt, instance->salt_size);
  seed_size = instance->seed_size;
  for (repetition = 0; repetition < instance->num_repetitions; repetition++)
  {
    hash_squeeze(&ctx, buffer, seed_size);
    memcpy(seeds + seed_size * repetition, buffer, seed_size);
  }
}

void commit_to_party_seed(const aimer_instance_t* instance,
                          const uint8_t* seed, const uint8_t* salt,
                          size_t repetition_index, size_t party_index,
                          uint8_t* commitment)
{
  hash_instance ctx;

  hash_init(&ctx, instance->digest_size);
  hash_update(&ctx, salt, instance->salt_size);
  hash_update_uint16(&ctx, (uint16_t)repetition_index);
  hash_update_uint16(&ctx, (uint16_t)party_index);
  hash_update(&ctx, seed, instance->seed_size);
  hash_final(&ctx);

  hash_squeeze(&ctx, commitment, instance->digest_size);
}

void phase_1_commitment(const aimer_instance_t* instance,
                        const signature_t* sig,
                        const aimer_publickey_t* public_key,
                        const uint8_t* message, size_t message_len,
                        const uint8_t* party_seed_commitments,
                        uint8_t* h_1)
{
  hash_instance ctx;
  uint8_t buffer[instance->field_size];

  memset(buffer, 0, sizeof(buffer));

  hash_init_prefix(&ctx, instance->digest_size, HASH_PREFIX_1);
  hash_update(&ctx, message, message_len);
  hash_update(&ctx, public_key->data, sizeof(public_key->data));
  hash_update(&ctx, sig->salt, instance->salt_size);

  for (size_t repetition = 0; repetition < instance->num_repetitions;
       repetition++)
  {
    proof_t* proof = &sig->proofs[repetition];
    for (size_t party = 0; party < instance->num_MPC_parties; party++)
    {
      hash_update(&ctx, party_seed_commitments +
                  ((repetition * instance->num_MPC_parties + party) *
                   instance->digest_size), instance->digest_size);
    }

    hash_update(&ctx, proof->input_delta, instance->aim_params.block_size);

    for (size_t ell = 0; ell < instance->aim_params.num_input_sboxes; ell++)
    {
      GF_to_bytes(proof->z_delta[ell], buffer);
      hash_update(&ctx, buffer, instance->field_size);
    }

    GF_to_bytes(proof->c_delta, buffer);
    hash_update(&ctx, (const uint8_t*)proof->c_delta, instance->field_size);
  }
  hash_final(&ctx);

  hash_squeeze(&ctx, h_1, instance->digest_size);
}

void phase_1_expand(const aimer_instance_t* instance, const uint8_t* h_1,
                    GF** epsilons)
{
  hash_instance ctx;
  uint8_t buffer[instance->field_size];

  memset(buffer, 0, sizeof(buffer));

  hash_init(&ctx, instance->digest_size);
  hash_update(&ctx, h_1, instance->digest_size);
  hash_final(&ctx);

  for (size_t e = 0; e < instance->num_repetitions; e++)
  {
    epsilons[e] =
      malloc((instance->aim_params.num_input_sboxes + 1) * instance->field_size);
    for (size_t ell = 0; ell < (instance->aim_params.num_input_sboxes + 1); ell++)
    {
      hash_squeeze(&ctx, buffer, instance->field_size);
      GF_from_bytes(buffer, epsilons[e][ell]);
    }
  }
}

void phase_2_commitment(const aimer_instance_t* instance, const uint8_t* salt,
                        const uint8_t* h_1, const GF* repetition_alpha_shares,
                        const GF** v_shares, uint8_t* h_2)
{
  hash_instance ctx;
  uint8_t buffer[instance->field_size];

  memset(buffer, 0, sizeof(buffer));

  hash_init_prefix(&ctx, instance->digest_size, HASH_PREFIX_2);
  hash_update(&ctx, salt, instance->salt_size);
  hash_update(&ctx, h_1 , instance->digest_size);

  for (size_t repetition = 0; repetition < instance->num_repetitions;
       repetition++)
  {
    for (size_t party = 0; party < instance->num_MPC_parties; party++)
    {
      GF_to_bytes(repetition_alpha_shares
                  [repetition * instance->num_MPC_parties + party],
                  buffer);
      hash_update(&ctx, buffer, instance->field_size);

      GF_to_bytes(v_shares[repetition][party], buffer);
      hash_update(&ctx, buffer, instance->field_size);
    }
  }
  hash_final(&ctx);

  hash_squeeze(&ctx, h_2, instance->digest_size);
}

void phase_2_expand(const aimer_instance_t* instance, const uint8_t* h_2,
                    uint16_t* opened_parties)
{
  size_t num_squeeze_bytes = instance->num_MPC_parties > 256 ? 2 : 1;
  uint16_t mask = (1ULL << ceil_log2(instance->num_MPC_parties)) - 1;
  uint16_t party;

  hash_instance ctx;
  hash_init(&ctx, instance->digest_size);
  hash_update(&ctx, h_2, instance->digest_size);
  hash_final(&ctx);

  for (size_t e = 0; e < instance->num_repetitions; e++)
  {
    do
    {
      hash_squeeze(&ctx, (uint8_t *)&party, num_squeeze_bytes);
      party = le16toh(party) & mask;
    }
    while (party >= instance->num_MPC_parties);

    opened_parties[e] = party;
  }
}
