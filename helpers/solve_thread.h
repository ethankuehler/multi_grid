

#ifndef SOLVE_MUTI_LIBRARY_H
#define SOLVE_MUTI_LIBRARY_H

#include "N_len.h"

void solve_threaded(const double* f, double* u, N_len Nlen, int iters, double w, double dx);

#endif //SOLVE_MUTI_LIBRARY_H

