// SPDX-License-Identifier: MIT

#ifndef HASH_H
#define HASH_H

#include <stdint.h>
#include "field.h"
#include "shake/KeccakHash.h"
#include "shake/KeccakHashtimes4.h"

static const uint8_t HASH_PREFIX_0    =  0;
static const uint8_t HASH_PREFIX_1    =  1;
static const uint8_t HASH_PREFIX_2    =  2;
static const uint8_t HASH_PREFIX_3    =  3;
static const uint8_t HASH_PREFIX_4    =  4;
static const uint8_t HASH_PREFIX_5    =  5;

typedef Keccak_HashInstance hash_instance;

void hash_init(hash_instance* ctx);
void hash_init_prefix(hash_instance* ctx, const uint8_t prefix);
void hash_update(hash_instance* ctx, const uint8_t* data, size_t data_byte_len);
void hash_update_GF(hash_instance* ctx, const GF g);
void hash_final(hash_instance* ctx);
void hash_squeeze(hash_instance* ctx, uint8_t* buffer, size_t buffer_len);

// x4 parallel hashing
typedef Keccak_HashInstancetimes4 hash_instance_x4;

void hash_init_x4(hash_instance_x4* ctx);
void hash_init_prefix_x4(hash_instance_x4* ctx, const uint8_t prefix);
void hash_update_x4(hash_instance_x4* ctx, const uint8_t** data, size_t data_byte_len);
void hash_update_x4_1(hash_instance_x4* ctx, const uint8_t *data,
                      size_t data_byte_len);
void hash_final_x4(hash_instance_x4* ctx);
void hash_squeeze_x4(hash_instance_x4* ctx, uint8_t** buffer, size_t buffer_len);

#endif // HASH_H
