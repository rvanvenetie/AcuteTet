#pragma once

#ifdef USE_OMP
#include <omp.h>
#else
#include <time.h>
#include <sys/time.h>
double omp_get_wtime();
int omp_get_thread_num();
#endif
