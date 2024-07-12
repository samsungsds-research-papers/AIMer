// -----------------------------------------------------------------------------
// File Name   : hash.c
// Description : 
// SPDX-License-Identifier: MIT
// -----------------------------------------------------------------------------

#include "hash.h"
#include "portable_endian.h"

void hash_init(hash_instance* ctx, size_t digest_size)
{
    if (digest_size <= 32)
    {
        Keccak_HashInitialize_SHAKE128(ctx);
    }
    else
    {
        Keccak_HashInitialize_SHAKE256(ctx);
    }
}

void hash_update(hash_instance* ctx, const uint8_t* data, size_t data_byte_len)
{
    Keccak_HashUpdate(ctx, data, data_byte_len << 3);
}

void hash_final(hash_instance* ctx)
{
    Keccak_HashFinal(ctx, NULL);
}

void hash_squeeze(hash_instance* ctx, uint8_t* data, size_t data_byte_len)
{
    Keccak_HashSqueeze(ctx, data, data_byte_len << 3);
}

void hash_update_uint16(hash_instance* ctx, uint16_t data)
{
    const uint16_t data_little_endian = htole16(data);
    hash_update(ctx, (const uint8_t*)&data_little_endian,
                sizeof(data_little_endian));
}

void hash_init_prefix(hash_instance* ctx, size_t digest_size,
                      const uint8_t prefix)
{
    hash_init(ctx, digest_size);
    hash_update(ctx, &prefix, sizeof(prefix));
}
