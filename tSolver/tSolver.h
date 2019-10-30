#include "../multi_grid/multi_grid.h"

#ifndef C_SOR_3D_TSOLVER_H
#define C_SOR_3D_TSOLVER_H

void tSolve(const double* f, double* u, N_len Nlen, int iter, double w, double dx);

#endif //C_SOR_3D_TSOLVER_H
