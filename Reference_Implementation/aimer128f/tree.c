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

void expand_tree(const uint8_t root_seed[AIMER_SEED_SIZE],
                 const uint8_t salt[AIMER_SALT_SIZE],
                 const size_t rep_index,
                 uint8_t nodes[2 * AIMER_NUM_MPC_PARTIES][AIMER_SEED_SIZE])
{
  memcpy(nodes[1], root_seed, AIMER_SEED_SIZE);

  for (size_t index = 1; index < AIMER_NUM_MPC_PARTIES; index++)
  {
    // parent: index | child: (2index, 2index + 1)
    expand_seed(nodes[index], salt, rep_index, index, nodes[index << 1]);
  }
}

void reveal_all_but(const uint8_t nodes[2 * AIMER_NUM_MPC_PARTIES][AIMER_SEED_SIZE],
                    const size_t cover_index,
                    uint8_t reveal_path[AIMER_NUM_MPC_PARTIES_LOG2][AIMER_SEED_SIZE])
{
  size_t index = cover_index + AIMER_NUM_MPC_PARTIES;
  for (size_t depth = 0; depth < AIMER_NUM_MPC_PARTIES_LOG2; depth++)
  {
    // index ^ 1 is sibling index
    memcpy(reveal_path[depth], nodes[index ^ 1], AIMER_SEED_SIZE);

    // go to parent node
    index >>= 1;
  }
}

void reconstruct_seed_tree(const uint8_t reveal_path[AIMER_NUM_MPC_PARTIES_LOG2][AIMER_SEED_SIZE],
                           const size_t cover_index,
                           const uint8_t salt[AIMER_SALT_SIZE],
                           const size_t rep_index,
                           uint8_t nodes[2 * AIMER_NUM_MPC_PARTIES][AIMER_SEED_SIZE])
{
  size_t index, depth;

  index = cover_index + AIMER_NUM_MPC_PARTIES;
  for (depth = 0; depth < AIMER_NUM_MPC_PARTIES_LOG2; depth++)
  {
    memcpy(nodes[index ^ 1], reveal_path[depth], AIMER_SEED_SIZE);
    index >>= 1;
  }

  for (depth = 1; depth < AIMER_NUM_MPC_PARTIES_LOG2; depth++)
  {
    // skip the parent nodes
    size_t skip_index = (cover_index + AIMER_NUM_MPC_PARTIES) >>
                        (AIMER_NUM_MPC_PARTIES_LOG2 - depth);
    for (index = (1U << depth); index < (2U << depth); index++)
    {
      if (index == skip_index)
      {
        continue;
      }

      expand_seed(nodes[index], salt, rep_index, index, nodes[index << 1]);
    }
  }
}
