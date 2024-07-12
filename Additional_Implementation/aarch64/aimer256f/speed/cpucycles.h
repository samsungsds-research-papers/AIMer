#ifndef CPUCYCLES_H
#define CPUCYCLES_H

#include <stdint.h>
#include <time.h>

static inline uint64_t cpucycles(void)
{

  struct timespec t;

  clock_gettime(CLOCK_REALTIME, &t);

  return (t.tv_sec * 1000000000 + t.tv_nsec);
}

uint64_t cpucycles_overhead(void);

#endif
