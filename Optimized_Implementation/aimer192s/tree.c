// SPDX-License-Identifier: MIT

#include <stdint.h>
#include <string.h>
#include "tree.h"
#include "hash.h"

void expand_seed(const uint8_t seed[AIMER_SEED_SIZE],
                 const uint8_t salt[AIMER_SALT_SIZE],
                 const size_t rep_index,
                 const size_t node_index,
                 uint8_t out[2 * AIMER_SEED_SIZE])
{
  hash_instance ctx;

  hash_init_prefix(&ctx, HASH_PREFIX_4);
  hash_update(&ctx, salt, AIMER_SALT_SIZE);
  hash_update(&ctx, (const uint8_t*)&rep_index, sizeof(uint8_t));
  hash_update(&ctx, (const uint8_t*)&node_index, sizeof(uint8_t));
  hash_update(&ctx, seed, AIMER_SEED_SIZE);
  hash_final(&ctx);
  hash_squeeze(&ctx, out, 2 * AIMER_SEED_SIZE);
}

//  Example of tree for [N = 8]
//  x
//  d = 0: 1
//  d = 1: 2         3
//  d = 2: 4   5     6     7
//  d = 3: 8 9 10 11 12 13 14 15

void expand_trees(const uint8_t salt[AIMER_SALT_SIZE],
                  uint8_t nodes[][2 * AIMER_NUM_MPC_PARTIES - 1][AIMER_SEED_SIZE])
{
  uint8_t bufs[4][AIMER_SEED_SIZE + 2];
  const uint8_t* in_ptrs[4] = {bufs[0], bufs[1], bufs[2], bufs[3]};
  uint8_t* out_ptrs[4];

  hash_instance_x4 ctx, ctx_;
  hash_init_prefix_x4(&ctx_, HASH_PREFIX_4);
  hash_update_x4_1(&ctx_, salt, AIMER_SALT_SIZE);

  size_t rep, index, depth;
  for (rep = 0; rep + 4 <= AIMER_NUM_REPETITIONS; rep += 4)
  {
    for (index = 1; index < AIMER_NUM_MPC_PARTIES; index++)
    {
      memcpy(&ctx, &ctx_, sizeof(hash_instance_x4));

      bufs[0][0] = (uint8_t)(rep + 0);
      bufs[1][0] = (uint8_t)(rep + 1);
      bufs[2][0] = (uint8_t)(rep + 2);
      bufs[3][0] = (uint8_t)(rep + 3);

      bufs[0][1] = (uint8_t)(index);
      bufs[1][1] = (uint8_t)(index);
      bufs[2][1] = (uint8_t)(index);
      bufs[3][1] = (uint8_t)(index);

      memcpy(bufs[0] + 2, nodes[rep + 0][index - 1], AIMER_SEED_SIZE);
      memcpy(bufs[1] + 2, nodes[rep + 1][index - 1], AIMER_SEED_SIZE);
      memcpy(bufs[2] + 2, nodes[rep + 2][index - 1], AIMER_SEED_SIZE);
      memcpy(bufs[3] + 2, nodes[rep + 3][index - 1], AIMER_SEED_SIZE);

      hash_update_x4(&ctx, in_ptrs, AIMER_SEED_SIZE + 2);
      hash_final_x4(&ctx);

      out_ptrs[0] = nodes[rep + 0][2 * index - 1];
      out_ptrs[1] = nodes[rep + 1][2 * index - 1];
      out_ptrs[2] = nodes[rep + 2][2 * index - 1];
      out_ptrs[3] = nodes[rep + 3][2 * index - 1];

      hash_squeeze_x4(&ctx, out_ptrs, AIMER_SEED_SIZE << 1);
    }
  }

  hash_instance ctx1, ctx1_;
  hash_init_prefix(&ctx1_, HASH_PREFIX_4);
  hash_update(&ctx1_, salt, AIMER_SALT_SIZE);

  for (; rep < AIMER_NUM_REPETITIONS; rep++)
  {
    for (depth = 0; depth < 2; depth++)
    {
      for (index = (1U << depth); index < (2U << depth); index++) 
      {
        memcpy(&ctx1, &ctx1_, sizeof(hash_instance));
        bufs[0][0] = (uint8_t)(rep);
        bufs[0][1] = (uint8_t)(index);
        memcpy(bufs[0] + 2, nodes[rep][index - 1], AIMER_SEED_SIZE);

        hash_update(&ctx1, bufs[0], AIMER_SEED_SIZE + 2);
        hash_final(&ctx1);
        hash_squeeze(&ctx1, nodes[rep][2 * index - 1], AIMER_SEED_SIZE << 1);
      }
    }
    // depth = 2;
    for (; depth < AIMER_NUM_MPC_PARTIES_LOG2; depth++)
    {
      for (index = (1U << depth); index < (2U << depth); index += 4)
      {
        memcpy(&ctx, &ctx_, sizeof(hash_instance_x4));

        bufs[0][0] = (uint8_t)(rep);
        bufs[1][0] = (uint8_t)(rep);
        bufs[2][0] = (uint8_t)(rep);
        bufs[3][0] = (uint8_t)(rep);

        bufs[0][1] = (uint8_t)(index + 0);
        bufs[1][1] = (uint8_t)(index + 1);
        bufs[2][1] = (uint8_t)(index + 2);
        bufs[3][1] = (uint8_t)(index + 3);

        memcpy(bufs[0] + 2, nodes[rep][index - 1], AIMER_SEED_SIZE);
        memcpy(bufs[1] + 2, nodes[rep][index + 0], AIMER_SEED_SIZE);
        memcpy(bufs[2] + 2, nodes[rep][index + 1], AIMER_SEED_SIZE);
        memcpy(bufs[3] + 2, nodes[rep][index + 2], AIMER_SEED_SIZE);

        hash_update_x4(&ctx, in_ptrs, AIMER_SEED_SIZE + 2);
        hash_final_x4(&ctx);

        out_ptrs[0] = nodes[rep][2 * index - 1];
        out_ptrs[1] = nodes[rep][2 * index + 1];
        out_ptrs[2] = nodes[rep][2 * index + 3];
        out_ptrs[3] = nodes[rep][2 * index + 5];

        hash_squeeze_x4(&ctx, out_ptrs, AIMER_SEED_SIZE << 1);
      }
    }
  }
}

void reveal_all_but(const uint8_t nodes[2 * AIMER_NUM_MPC_PARTIES - 1][AIMER_SEED_SIZE],
                    const size_t cover_index,
                    uint8_t reveal_path[AIMER_NUM_MPC_PARTIES_LOG2][AIMER_SEED_SIZE])
{
  size_t index = cover_index + AIMER_NUM_MPC_PARTIES;
  for (size_t depth = 0; depth < AIMER_NUM_MPC_PARTIES_LOG2; depth++)
  {
    // index ^ 1 is sibling index
    memcpy(reveal_path[depth], nodes[(index ^ 1) - 1], AIMER_SEED_SIZE);

    // go to parent node
    index >>= 1;
  }
}

void reconstruct_seed_tree(const uint8_t reveal_path[AIMER_NUM_MPC_PARTIES_LOG2][AIMER_SEED_SIZE],
                           const size_t cover_index,
                           const uint8_t salt[AIMER_SALT_SIZE],
                           const size_t rep,
                           uint8_t nodes[2 * AIMER_NUM_MPC_PARTIES - 2][AIMER_SEED_SIZE])
{
  uint8_t bufs[4][AIMER_SEED_SIZE + 2];
  const uint8_t* in_ptrs[4] = {bufs[0], bufs[1], bufs[2], bufs[3]};
  uint8_t* out_ptrs[4];

  hash_instance_x4 ctx, ctx_;
  hash_init_prefix_x4(&ctx_, HASH_PREFIX_4);
  hash_update_x4_1(&ctx_, salt, AIMER_SALT_SIZE);

  size_t index, depth, path;

  // d = 1
  path = ((cover_index + AIMER_NUM_MPC_PARTIES)
           >> (AIMER_NUM_MPC_PARTIES_LOG2 - 1)) ^ 1;

  expand_seed(reveal_path[AIMER_NUM_MPC_PARTIES_LOG2 - 1],
              salt, rep, path, nodes[(path << 1) - 2]);

  for (depth = 2; depth < AIMER_NUM_MPC_PARTIES_LOG2; depth++)
  {
    path = ((cover_index + AIMER_NUM_MPC_PARTIES)
             >> (AIMER_NUM_MPC_PARTIES_LOG2 - depth)) ^ 1;

    memcpy(nodes[path - 2], reveal_path[AIMER_NUM_MPC_PARTIES_LOG2 - depth],
           AIMER_SEED_SIZE);

    for (index = (1U << depth); index < (2U << depth); index += 4)
    {
      memcpy(&ctx, &ctx_, sizeof(hash_instance_x4));

      bufs[0][0] = (uint8_t)(rep);
      bufs[1][0] = (uint8_t)(rep);
      bufs[2][0] = (uint8_t)(rep);
      bufs[3][0] = (uint8_t)(rep);

      bufs[0][1] = (uint8_t)(index + 0);
      bufs[1][1] = (uint8_t)(index + 1);
      bufs[2][1] = (uint8_t)(index + 2);
      bufs[3][1] = (uint8_t)(index + 3);

      memcpy(bufs[0] + 2, nodes[index - 2], AIMER_SEED_SIZE);
      memcpy(bufs[1] + 2, nodes[index - 1], AIMER_SEED_SIZE);
      memcpy(bufs[2] + 2, nodes[index + 0], AIMER_SEED_SIZE);
      memcpy(bufs[3] + 2, nodes[index + 1], AIMER_SEED_SIZE);

      hash_update_x4(&ctx, in_ptrs, AIMER_SEED_SIZE + 2);
      hash_final_x4(&ctx);

      out_ptrs[0] = nodes[2 * index - 2];
      out_ptrs[1] = nodes[2 * index + 0];
      out_ptrs[2] = nodes[2 * index + 2];
      out_ptrs[3] = nodes[2 * index + 4];

      hash_squeeze_x4(&ctx, out_ptrs, AIMER_SEED_SIZE << 1);
    }
  }

  path = ((cover_index + AIMER_NUM_MPC_PARTIES)
          >> (AIMER_NUM_MPC_PARTIES_LOG2 - depth)) ^ 1;

  memcpy(nodes[path - 2], reveal_path[AIMER_NUM_MPC_PARTIES_LOG2 - depth],
          AIMER_SEED_SIZE);
}
