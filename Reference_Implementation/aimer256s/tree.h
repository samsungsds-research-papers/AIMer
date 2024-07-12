// SPDX-License-Identifier: MIT

#ifndef TREE_H
#define TREE_H

#include <stdint.h>
#include "params.h"

void expand_tree(const uint8_t root_seed[AIMER_SEED_SIZE],
                 const uint8_t salt[AIMER_SALT_SIZE],
                 const size_t rep_index,
                 uint8_t nodes[2 * AIMER_NUM_MPC_PARTIES][AIMER_SEED_SIZE]);

void reveal_all_but(const uint8_t nodes[2 * AIMER_NUM_MPC_PARTIES][AIMER_SEED_SIZE],
                    const size_t cover_index,
                    uint8_t reveal_path[AIMER_NUM_MPC_PARTIES_LOG2][AIMER_SEED_SIZE]);

void reconstruct_seed_tree(const uint8_t reveal_path[AIMER_NUM_MPC_PARTIES_LOG2][AIMER_SEED_SIZE],
                           const size_t cover_index,
                           const uint8_t salt[AIMER_SALT_SIZE],
                           const size_t rep_index,
                           uint8_t nodes[2 * AIMER_NUM_MPC_PARTIES][AIMER_SEED_SIZE]);

#endif // TREE_H
