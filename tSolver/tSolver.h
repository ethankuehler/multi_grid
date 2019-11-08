#include "../multi_grid/multi_grid.h"

#ifndef C_SOR_3D_TSOLVER_H
#define C_SOR_3D_TSOLVER_H

void tSolve(const float* f, float* u, N_len Nlen, int iter, float w, float dx);

#endif //C_SOR_3D_TSOLVER_H
