#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "solve_block.h"
#include "operators.h"
#include "multi_grid.h"


N_len coarsen(N_len Nlen) {
    int i = (Nlen.i+1)/2;
    int j = (Nlen.k+1)/2;
    int k = (Nlen.k+1)/2;
    assert(!(i < 3 || j < 3 || k < 3));
    return (N_len){i, j, k};
}

N_len refine(N_len Nlen){
    return (N_len){Nlen.i*2 - 1, Nlen.j*2 - 1, Nlen.k*2 - 1};
}

int length(N_len Nlen) {
    return Nlen.i * Nlen.j * Nlen.k;
}

char can_coarsen(N_len Nlen) {
    int i = (Nlen.i-1)/2;
    int j = (Nlen.k-1)/2;
    int k = (Nlen.k-1)/2;
    if (!(i < 3 || j < 3 || k < 3)) {
        return 1;
    }
    return 0;
}

/* f is the right hand side of the equation
 * u is the output of the function and is the potential
 * m is the int that determents the dimensions of the input
 * dx is the step size
 */
void multi(double* f, double* u, int m, double dx, double w, int iters) {
    int N = (int) pow(2, m) + 1;
    int Nc = (int) pow(2, m - 1) + 1;
    N_len Nlen = (N_len){N, N, N};
    N_len Nclen = coarsen(Nlen);
    if (m == 1) {
        int n = 1 + N + N*N;
        u[n] = (1.0/-6.0)*(f[n]*dx - (u[n - N*N] + u[n - N] + u[n + 1] + u[n - 1] + u[n + N] + u[n + N*N]));
    } else {
        double* fc = calloc(sizeof(double), length(Nclen));
        double* uc = calloc(sizeof(double), length(Nclen));
        double* d = calloc(sizeof(double), length(Nlen));
        if (m != 8) {
            solve(f, u, Nlen, iters, 1, dx);
        } else {
            solve(f, u, Nlen, iters, w, dx);
        }

        restriction(f, u, fc, Nlen, dx*dx);
        multi(fc, uc, m - 1, dx*2, w, iters);
        interpolate(uc, d, Nlen, dx);
        for (int i = 0; i < N*N*N; i++) {
            u[i] = u[i] + d[i];
        }
        if (m != 8) {
            solve(f, u, Nlen, iters, 1, dx);
        } else {
            solve(f, u, Nlen, iters, w, dx);
        }
        free(fc);
        free(uc);
        free(d);
    }
}
