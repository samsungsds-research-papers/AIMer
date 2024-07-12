// SPDX-License-Identifier: MIT

#ifndef FIELD_H
#define FIELD_H

#include <stdint.h>
#include "params.h"

typedef uint64_t GF[2];

void GF_set0(GF a);
void GF_copy(const GF in, GF out);
void GF_to_bytes(const GF in, uint8_t* out);
void GF_from_bytes(const uint8_t* in, GF out);

void GF_add(const GF a, const GF b, GF c);
void GF_mul(const GF a, const GF b, GF c);
void GF_mul_N(const GF a[AIMER_NUM_MPC_PARTIES], const GF b,
                    GF c[AIMER_NUM_MPC_PARTIES]);
void GF_mul_add(const GF a, const GF b, GF c);
void GF_mul_add_N(const GF a[AIMER_NUM_MPC_PARTIES], const GF b,
                        GF c[AIMER_NUM_MPC_PARTIES]);
void GF_sqr(const GF a, GF c);
void GF_sqr_N(const GF a[AIMER_NUM_MPC_PARTIES], GF c[AIMER_NUM_MPC_PARTIES]);
void GF_transposed_matmul(const GF a, const GF b[AIM2_NUM_BITS_FIELD], GF c);
void GF_transposed_matmul_add_N(const GF a[AIMER_NUM_MPC_PARTIES],
                                const GF b[AIM2_NUM_BITS_FIELD],
                                GF c[AIMER_NUM_MPC_PARTIES]);

void POLY_mul_add_N(const GF a[AIMER_NUM_MPC_PARTIES], const GF b,
                    GF lo[AIMER_NUM_MPC_PARTIES],
                    GF hi[AIMER_NUM_MPC_PARTIES]);
void POLY_red_N(const GF hi[AIMER_NUM_MPC_PARTIES],
                      GF lo[AIMER_NUM_MPC_PARTIES]);

#endif // FIELD_H
