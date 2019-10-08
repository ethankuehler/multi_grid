//
// Created by dcrush on 10/3/19.
//

#ifndef C_SOR_3D_MULTI_GRID_H
#define C_SOR_3D_MULTI_GRID_H

typedef struct _N_len {
    int i;
    int j;
    int k;
} N_len;

N_len coarsen(N_len Nlen);
N_len refine(N_len Nlen);
int length(N_len Nlen);
char can_coarsen(N_len Nlen);

void multi(double* f, double* u, int m, double dx, double w, int iters);


#endif //C_SOR_3D_MULTI_GRID_H


