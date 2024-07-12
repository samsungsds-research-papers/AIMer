// -----------------------------------------------------------------------------
// File Name   : test_aim.c
// Description : 
// SPDX-License-Identifier: MIT
// -----------------------------------------------------------------------------

#include "aim.h"
#include <stdio.h>

int main()
{
  int i;

  // Inverted order compared to Sage!
  uint8_t pt[32] =
    {0xff, 0xee, 0xdd, 0xcc, 0xbb, 0xaa, 0x99, 0x88,
     0x77, 0x66, 0x55, 0x44, 0x33, 0x22, 0x11, 0x00,
     0xff, 0xee, 0xdd, 0xcc, 0xbb, 0xaa, 0x99, 0x88,
     0x77, 0x66, 0x55, 0x44, 0x33, 0x22, 0x11, 0x00};
  uint8_t ct_expected[32] =
    {0xcc, 0x2b, 0xd0, 0x80, 0xf1, 0xde, 0x58, 0x5a,
     0xd8, 0x45, 0xe6, 0x96, 0xd9, 0xcd, 0x83, 0x48,
     0x0a, 0xa7, 0xb4, 0xcb, 0xa4, 0x1a, 0xfd, 0x97,
     0x90, 0x37, 0x88, 0xf0, 0xbc, 0x9f, 0xb5, 0xd6};
  uint8_t iv[32] =
    {0xff, 0xee, 0xdd, 0xcc, 0xbb, 0xaa, 0x99, 0x88,
     0x77, 0x66, 0x55, 0x44, 0x33, 0x22, 0x11, 0x00,
     0xff, 0xee, 0xdd, 0xcc, 0xbb, 0xaa, 0x99, 0x88,
     0x77, 0x66, 0x55, 0x44, 0x33, 0x22, 0x11, 0x00};

  uint8_t ct[32] = {0,};

  aim(pt, iv, ct);

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

    if (ct[i] != ct_expected[i])
    {
      printf("\nError!\n");
      return 0;
    }
  }  
  printf("\nPassed!\n");

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
