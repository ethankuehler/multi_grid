#include "../multi_grid/multi_grid.h"
#include "../multi_grid/operators.h"
#include <stdlib.h>
#include <string.h>


void solver_red(const double* f, double* u, N_len Nlen, double w, double dxs) {
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    int i;
#pragma omp parallel for num_threads(8) schedule(auto) shared(u, f, Nlen, w, dxs) private(i)
    for (i = 1; i < Ni - 1; i++) {
        for (int j = 1; j < Nj - 1; j++) {
            for (int k = 1; k < Nk - 1; k++) {
                int n = loc(i, j, k, Nlen);
                u[n] = (1 - w)*u[n] +
                       (w/-6)*
                       (f[n]*dxs - (u[n - Ni*Nj] + u[n - Nj] + u[n + 1] + u[n - 1] + u[n + Nj] + u[n + Ni*Nj]));
            }
        }
    }

}

void solver_copy_boundary(const double* u, double* u_out, N_len Nlen) {
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    int i;
    for (i = 1; i < Ni - 1; i++) {
        int n = 0;
        for(int j = 0; j < Nj; j++) {
            n = loc(i, j, 0, Nlen);
            u_out[n] = u[n];
            n = loc(i, j, Nk, Nlen);
            u_out[n] = u[n];
        }
        for(int k = 0; k < Nk; k++) {
            n = loc(i, 0, k, Nlen);
            u_out[n] = u[n];
            n = loc(i, Nj, k, Nlen);
            u_out[n] = u[n];
        }
    }
}

void tSolve(const double* f, double* u, N_len Nlen, int iter, double w, double dx) {
    for(int _= 0; _ < iter; _++) {
        solver_red(f, u, Nlen, w, dx*dx);
    }
}
