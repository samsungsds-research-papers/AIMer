// SPDX-License-Identifier: MIT

#include <stdio.h>
#include <string.h>
#include "../field.h"
#include "../aim2.h"

int main()
{
  int succ = 1, i;

  // Inverted order compared to Sage!
  uint8_t pt[24] =
    {0x77, 0x66, 0x55, 0x44, 0x33, 0x22, 0x11, 0x00,
     0xff, 0xee, 0xdd, 0xcc, 0xbb, 0xaa, 0x99, 0x88,
     0x77, 0x66, 0x55, 0x44, 0x33, 0x22, 0x11, 0x00};
  const uint8_t ct_expected[24] =
    {0x0a, 0x68, 0x3f, 0x52, 0x45, 0xc5, 0x68, 0xb4,
     0xc6, 0x4a, 0xfc, 0x73, 0x4a, 0x51, 0x94, 0xc0,
     0x77, 0xc8, 0x41, 0x35, 0xdf, 0xa5, 0xce, 0xad};
  uint8_t iv[24] =
    {0x77, 0x66, 0x55, 0x44, 0x33, 0x22, 0x11, 0x00,
     0xff, 0xee, 0xdd, 0xcc, 0xbb, 0xaa, 0x99, 0x88,
     0x77, 0x66, 0x55, 0x44, 0x33, 0x22, 0x11, 0x00};

  uint8_t ct[24] = {0,};

  aim2(ct, pt, iv);

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
