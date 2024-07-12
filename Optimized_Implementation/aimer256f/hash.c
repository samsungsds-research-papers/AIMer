// SPDX-License-Identifier: MIT

#include <stdint.h>
#include "hash.h"
#include "portable_endian.h"

void hash_init(hash_instance* ctx)
{
  Keccak_HashInitialize_SHAKE256(ctx);
}

void hash_update(hash_instance* ctx, const uint8_t* data, size_t data_byte_len)
{
  Keccak_HashUpdate(ctx, data, data_byte_len << 3);
}

void hash_final(hash_instance* ctx)
{
  Keccak_HashFinal(ctx, NULL);
}

void hash_squeeze(hash_instance* ctx, uint8_t* buffer, size_t buffer_len)
{
  Keccak_HashSqueeze(ctx, buffer, buffer_len << 3);
}

void hash_update_GF(hash_instance* ctx, const GF g)
{
  Keccak_HashUpdate(ctx, (const uint8_t*)g, AIM2_NUM_BYTES_FIELD << 3);
}

void hash_init_prefix(hash_instance* ctx, const uint8_t prefix)
{
  hash_init(ctx);
  hash_update(ctx, &prefix, sizeof(prefix));
}

// x4 parallel hashing
void hash_init_x4(hash_instance_x4* ctx)
{
  Keccak_HashInitializetimes4_SHAKE256(ctx);
}

void hash_update_x4(hash_instance_x4* ctx, const uint8_t** data, size_t data_byte_len)
{
  Keccak_HashUpdatetimes4(ctx, data, data_byte_len << 3);
}

void hash_update_x4_1(hash_instance_x4* ctx, const uint8_t *data,
                      size_t data_byte_len)
{
  const uint8_t* temp[4] = {data, data, data, data};
  Keccak_HashUpdatetimes4(ctx, temp, data_byte_len << 3);
}

void hash_init_prefix_x4(hash_instance_x4* ctx, const uint8_t prefix)
{
  hash_init_x4(ctx);
  hash_update_x4_1(ctx, &prefix, sizeof(prefix));
}

void hash_final_x4(hash_instance_x4* ctx)
{
  Keccak_HashFinaltimes4(ctx, NULL);
}

void hash_squeeze_x4(hash_instance_x4* ctx, uint8_t** buffer, size_t buffer_len)
{
  Keccak_HashSqueezetimes4(ctx, buffer, buffer_len << 3);
}
