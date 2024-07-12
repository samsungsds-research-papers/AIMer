// -----------------------------------------------------------------------------
// File Name   : aim.h
// Description : 
// SPDX-License-Identifier: MIT
// -----------------------------------------------------------------------------

#ifndef AIM_H
#define AIM_H

#include "aimer.h"

#if       _AIMER_L == 1
  #include "aim128.h"

  #define aim                   aim128
  #define aim_with_sbox_output  aim128_with_sbox_output
  #define aim_mpc               aim128_mpc

#elif     _AIMER_L == 3
  #include "aim192.h"

  #define aim                   aim192
  #define aim_with_sbox_output  aim192_with_sbox_output
  #define aim_mpc               aim192_mpc

#elif     _AIMER_L == 5
  #include "aim256.h"

  #define aim                   aim256
  #define aim_with_sbox_output  aim256_with_sbox_output
  #define aim_mpc               aim256_mpc

#endif

typedef struct aim_params_t
{
  uint32_t block_size;
  uint32_t num_input_sboxes;
} aim_params_t;

#endif // AIM_H
