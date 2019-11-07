#include "../multi_grid/multi_grid.h"
#include "../multi_grid/operators.h"
#include <stdlib.h>
#include <string.h>


void point_sol(double* u, const double* f, int i, int j, int k, double dxs, double w, N_len Nlen) {
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    int n = loc(i, j, k, Nlen);
    u[n] = (1 - w)*u[n] +
           (w/-6)*(f[n]*dxs - (u[n - Ni*Nj] + u[n - Nj] + u[n + 1] + u[n - 1] + u[n + Nj] + u[n + Ni*Nj]));
}

void solver_red(const double* f, double* u, N_len Nlen, double w, double dxs) {
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    int i, j, k;
#pragma omp parallel num_threads(4) shared(u, f, Nlen, w, dxs)
    {
#pragma omp for  schedule(static) private(i, j, k)
        for (i = 1; i < Ni - 1; i += 2) {
            for (j = 1; j < Nj - 1; j += 2) {
                for (k = 1; k < Nk - 1; k += 2) {
                    point_sol(u, f, i, j, k, dxs, w, Nlen);
                }
            }
            for (j = 2; j < Nj - 1; j += 2) {
                for (k = 2; k < Nk - 1; k += 2) {
                    point_sol(u, f, i, j, k, dxs, w, Nlen);
                }
            }
        }
#pragma omp for  schedule(static) private(i, j, k)
        for (i = 2; i < Ni - 1; i += 2) {
            for (j = 1; j < Nj - 1; j += 2) {
                for (k = 2; k < Nk - 1; k += 2) {
                    point_sol(u, f, i, j, k, dxs, w, Nlen);
                }
            }
            for (j = 2; j < Nj - 1; j += 2) {
                for (k = 1; k < Nk - 1; k += 2) {
                    point_sol(u, f, i, j, k, dxs, w, Nlen);
                }
            }
        }
    }
}

void solver_black(const double* f, double* u, N_len Nlen, double w, double dxs) {
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    int i, j, k;
#pragma omp parallel num_threads(4) shared(u, f, Nlen, w, dxs)
    {
#pragma omp for  schedule(static) private(i, j, k)
        for (i = 1; i < Ni - 1; i += 2) {
            for (j = 1; j < Nj - 1; j += 2) {
                for (k = 2; k < Nk - 1; k += 2) {
                    point_sol(u, f, i, j, k, dxs, w, Nlen);
                }
            }
            for (j = 2; j < Nj - 1; j += 2) {
                for (k = 1; k < Nk - 1; k += 2) {
                    point_sol(u, f, i, j, k, dxs, w, Nlen);
                }
            }
        }
#pragma omp for  schedule(static) private(i, j, k)
        for (i = 2; i < Ni - 1; i += 2) {
            for (j = 1; j < Nj - 1; j += 2) {
                for (k = 1; k < Nk - 1; k += 2) {
                    point_sol(u, f, i, j, k, dxs, w, Nlen);
                }
            }
            for (j = 2; j < Nj - 1; j += 2) {
                for (k = 2; k < Nk - 1; k += 2) {
                    point_sol(u, f, i, j, k, dxs, w, Nlen);
                }
            }
        }
    }
}


void tSolve(const double* f, double* u, N_len Nlen, int iter, double w, double dx) {
    for(int _= 0; _ < iter; _++) {
        solver_red(f, u, Nlen, w, dx*dx);
        solver_black(f, u, Nlen, w, dx*dx);
    }
}

