// SPDX-License-Identifier: MIT

#include <stdio.h>
#include <string.h>
#include "../field.h"
#include "../aim2.h"

int main()
{
  int succ = 1, i;
  GF a, b, c;

  for (i = 0; i < AIM2_NUM_BYTES_FIELD; i++)
  {
    ((uint8_t*)a)[i] = i;
  }

  GF_exp_power_of_2(a, b, aim2_invmer_exponents[0]);
  GF_transposed_matmul(a, aim2_e1_power_matrix, c);
  if (memcmp(b, c, sizeof(GF)))
  {
    succ = 0;
    printf("[-] GF_exp_power_of_2 != GF_exp_by_mat for e1\n");
  }

  GF_exp_power_of_2(a, b, aim2_invmer_exponents[1]);
  GF_transposed_matmul(a, aim2_e2_power_matrix, c);
  if (memcmp(b, c, sizeof(GF)))
  {
    succ = 0;
    printf("[-] GF_exp_power_of_2 != GF_exp_by_mat for e2\n");
  }

  GF_exp(a, b, aim2_exponents[0]);                   // b ^ {2 ^ e - 1} = a
  GF_mul(a, b, c);                                   // c = a * b
  GF_exp_power_of_2(b, b, aim2_invmer_exponents[0]); // b = b ^ {2 ^ e}

  if (memcmp(b, c, sizeof(GF)))
  {
    succ = 0;
    printf("[-] aim2_exps[0] does not fit with aim2_invmer_exps[0]\n");
  }

  GF_exp(a, b, aim2_exponents[1]);                   // b ^ {2 ^ e - 1} = a
  GF_mul(a, b, c);                                   // c = a * b
  GF_exp_power_of_2(b, b, aim2_invmer_exponents[1]); // b = b ^ {2 ^ e}

  if (memcmp(b, c, sizeof(GF)))
  {
    succ = 0;
    printf("[-] aim2_exps[1] does not fit with aim2_invmer_exps[1]\n");
  }

  // Inverted order compared to Sage!
  uint8_t pt[16] =
    {0xff, 0xee, 0xdd, 0xcc, 0xbb, 0xaa, 0x99, 0x88,
     0x77, 0x66, 0x55, 0x44, 0x33, 0x22, 0x11, 0x00};
  const uint8_t ct_expected[16] =
    {0x18, 0xeb, 0x19, 0xe4, 0xfa, 0x63, 0x54, 0x58,
     0xe8, 0x74, 0x66, 0x7a, 0xf1, 0x44, 0xc1, 0x5a};
  uint8_t iv[16] =
    {0xff, 0xee, 0xdd, 0xcc, 0xbb, 0xaa, 0x99, 0x88,
     0x77, 0x66, 0x55, 0x44, 0x33, 0x22, 0x11, 0x00};

  uint8_t ct[16] = {0,};

  aim2(pt, iv, ct);

  printf("PLAINTEXT                 : ");
  for (i = (int)sizeof(pt) - 1; i >= 0; i--)
  {
    printf("%02x", pt[i]);
  }
  printf("\n");

  printf("IV                        : ");
  for (i = (int)sizeof(iv) - 1; i >= 0; i--)
  {
    printf("%02x", iv[i]);
  }
  printf("\n");

  printf("CIPHERTEXT                : ");
  for (i = (int)sizeof(ct) - 1; i >= 0; i--)
  {
    printf("%02x", ct[i]);
  }
  printf("\n");

  if (memcmp(ct_expected, ct, sizeof(GF)))
  {
    succ = 0;
    printf("[-] ct != ct_expected\n");
  }

  if (succ)
    printf("Passed!\n");
  else
    printf("Error!\n");

  // For  extracting test data -------------------------------------------------
  // printf("CIPHERTEXT(little endian) : ");
  // for (i = 0; i < (int)sizeof(ct) - 1; i++)
  // {
  //   printf("0x%02x, ", ct[i]);
  // }
  // printf("0x%02x\n", ct[i]);
  // For  extracting test data -------------------------------------------------

  return 0;
}
