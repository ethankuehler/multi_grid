#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "N_len.h"


int loc(int i, int j, int k, N_len Nlen) {
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    assert(i*Ni*Nj + j*Nj + k < Ni*Nj*Nk);
    return i*Ni*Nj + j*Nj + k;
}


//N is the dim of the input
void reduce(const double* f_in, double* f_out, N_len Nlen) {
    N_len Nclen = coarsen(Nlen);
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    int Nci = Nclen.i;
    int Ncj = Nclen.j;
    //iterating though all interior points
    for (int i = 1; i < Ni - 1; i++) {
        for (int j = 1; j < Nj - 1; j++) {
            for (int k = 1; k < Nk - 1; k++) {
                //taking the
                f_out[loc(i/2, j/2, k/2, Nclen)] =
                        (f_in[loc(i, j, k, Nlen)] +
                         f_in[loc(i + 1, j, k, Nlen)] +
                         f_in[loc(i - 1, j, k, Nlen)] +
                         f_in[loc(i, j + 1, k, Nlen)] +
                         f_in[loc(i, j - 1, k, Nlen)] +
                         f_in[loc(i, j, k + 1, Nlen)] +
                         f_in[loc(i, j, k - 1, Nlen)])/7;
            }
        }
    }
}

void residual(const double* f, const double* u, double* r, N_len Nlen, double dxs) {
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    //because the outer points are fixed by the boundary conditions, their residual is zero.
    for (int i = 1; i < Ni - 1; i++) {
        for (int j = 1; j < Nj - 1; j++) {
            for (int k = 1; k < Nk - 1; k++) {
                int n = loc(i, j, k, Nlen);
                r[n] = f[n] -
                       (u[n - Ni*Nj] + u[n - Nj] + u[n + 1] + -6*u[n] + u[n - 1] + u[n + Nj] + u[n + Ni*Nj])/dxs;
            }
        }
    }
}

void restriction(const double* f, const double* u, double* f_out, N_len Nlen, double dxs) {
    //compute residual
    double* r = calloc(sizeof(double), length(Nlen)); //freed at end of function

    //compute residual.
    residual(f, u, r, Nlen, dxs);
    reduce(r, f_out, Nlen);
    free(r);
}


void interpolate(const double* f, double* f_out, N_len Nlen, double dx) {
    N_len Nclen = coarsen(Nlen);
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    int Nci = Nclen.i;
    int Ncj = Nclen.j;
    //middle area and other walls.
    //we only need to do the interior areas as their should be zero
    for (int i = 0; i < Ni - 2; i += 2) {
        for (int j = 0; j < Nj - 2; j += 2) {
            for (int k = 0; k < Nk - 2; k += 2) {
                int n = loc(i, j, k, Nlen);
                int nc = loc(i/2, j/2, k/2, Nclen);

                //last corner works
                f_out[loc(i + 2, j + 2, k + 2, Nlen)] = f[loc(i/2 + 1, j/2 + 1, k/2 + 1, Nclen)];

                //line parts
                //top right  working
                f_out[loc(i + 2, j + 1, k + 2, Nlen)] =
                        f_out[loc(i + 2, j, k + 2, Nlen)]*(1 - dx) + f_out[loc(i + 2, j + 2, k + 2, Nlen)]*dx;

                //top forward working
                f_out[loc(i + 2, j + 2, k + 1, Nlen)] =
                        f_out[loc(i + 2, j + 2, k, Nlen)]*(1 - dx) + f_out[loc(i + 2, j + 2, k + 2, Nlen)]*dx;

                //forward right working
                f_out[loc(i + 1, j + 2, k + 2, Nlen)] =
                        f_out[loc(i, j + 2, k + 2, Nlen)]*(1 - dx) + f_out[loc(i + 2, j + 2, k + 2, Nlen)]*dx;

                //center of walls
                //right working
                f_out[loc(i + 1, j + 1, k + 2, Nlen)] =
                        f_out[loc(i + 1, j, k + 2, Nlen)]*(1 - dx) + f_out[loc(i + 1, j + 2, k + 2, Nlen)]*dx;

                // forward working
                f_out[loc(i + 1, j + 2, k + 1, Nlen)] =
                        f_out[loc(i + 1, j + 2, k, Nlen)]*(1 - dx) + f_out[loc(i + 1, j + 2, k + 2, Nlen)]*dx;

                // top
                f_out[loc(i + 2, j + 1, k + 1, Nlen)] =
                        f_out[loc(i + 2, j + 1, k, Nlen)]*(1 - dx) + f_out[loc(i + 2, j + 1, k + 2, Nlen)]*dx;

                //middle point working
                f_out[loc(i + 1, j + 1, k + 1, Nlen)] =
                        f_out[loc(i + 1, j + 1, k, Nlen)]*(1 - dx) + f_out[loc(i + 1, j + 1, k + 2, Nlen)]*dx;
            }
        }
    }
}
