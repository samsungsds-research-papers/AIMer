// -----------------------------------------------------------------------------
// File Name   : field.h
// Description : 
// SPDX-License-Identifier: MIT
// -----------------------------------------------------------------------------

#ifndef FIELD_H
#define FIELD_H

#include "../aimer.h"

#if       _AIMER_L == 1
  #include "field128.h"

  #define GF                GF2_128
  #define GF_set0           GF2_128_set0
  #define GF_set1           GF2_128_set1
  #define GF_copy           GF2_128_copy
  #define GF_to_bytes       GF2_128_to_bytes
  #define GF_from_bytes     GF2_128_from_bytes
  #define GF_add            GF2_128_add
  #define GF_mul            GF2_128_mul
  #define GF_sqr            GF2_128_sqr
  #define GF_matmul_vec     GF2_128_matmul_vec
  #define GF_matmul_mat     GF2_128_matmul_mat

#elif     _AIMER_L == 3
  #include "field192.h"

  #define GF                GF2_192
  #define GF_set0           GF2_192_set0
  #define GF_set1           GF2_192_set1
  #define GF_copy           GF2_192_copy
  #define GF_to_bytes       GF2_192_to_bytes
  #define GF_from_bytes     GF2_192_from_bytes
  #define GF_add            GF2_192_add
  #define GF_mul            GF2_192_mul
  #define GF_sqr            GF2_192_sqr
  #define GF_matmul_vec     GF2_192_matmul_vec
  #define GF_matmul_mat     GF2_192_matmul_mat

#elif     _AIMER_L == 5
  #include "field256.h"
  
  #define GF                GF2_256
  #define GF_set0           GF2_256_set0
  #define GF_set1           GF2_256_set1
  #define GF_copy           GF2_256_copy
  #define GF_to_bytes       GF2_256_to_bytes
  #define GF_from_bytes     GF2_256_from_bytes
  #define GF_add            GF2_256_add
  #define GF_mul            GF2_256_mul
  #define GF_sqr            GF2_256_sqr
  #define GF_matmul_vec     GF2_256_matmul_vec
  #define GF_matmul_mat     GF2_256_matmul_mat

#endif

#endif // FIELD_H
