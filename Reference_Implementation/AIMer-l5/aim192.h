// -----------------------------------------------------------------------------
// File Name   : aim192.h
// Description : 
// SPDX-License-Identifier: MIT
// -----------------------------------------------------------------------------

#ifndef AIM192_H
#define AIM192_H

#include "field/field.h"

#define TAPE_LEN 4752

void generate_linear_layer(const uint8_t* initial_vector, GF** matrix_A, GF vector_B);
void aim192(const uint8_t* input, const uint8_t* initial_vector,
            uint8_t* output);
void aim192_with_sbox_output(const uint8_t* input, const GF** matrix_A,
                             const GF vector_B, uint8_t* output,
                             GF* sbox_outputs);
void aim192_mpc(const uint8_t* input, const GF** matrix_A,
                const GF vector_B, const uint8_t* output,
                const size_t num_parties, GF *t_shares, GF *s_shares);

#endif // AIM192_H
