#include "parallel.h"

#ifndef USE_OMP
double omp_get_wtime()
{
  struct timeval time;
  gettimeofday(&time, NULL);
  return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

int omp_get_thread_num() {
  return 0;
}

#endif
