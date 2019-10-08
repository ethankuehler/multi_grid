#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
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

int length(N_len Nlen) {
    return Nlen.i * Nlen.j * Nlen.k;
}

char can_coarsen(N_len Nlen) {
    int i = (Nlen.i+1)/2;
    int j = (Nlen.k+1)/2;
    int k = (Nlen.k+1)/2;
    if (!(i < 3 || j < 3 || k < 3)) {
        return true;
    }
    return false;
}

/* f is the right hand side of the equation
 * u is the output of the function and is the potential
 * m is the int that determents the dimensions of the input
 * dx is the step size
 */
void multi(double* f, double* u, N_len Nlen, double dx, double w, int iters, char top){
    
    if (!can_coarsen(Nlen)) {
        solve(f, u, Nlen, iters, 1, dx);
    } else {
        N_len Nclen = coarsen(Nlen);
        double* fc = calloc(sizeof(double), length(Nclen));
        double* uc = calloc(sizeof(double), length(Nclen));
        double* d = calloc(sizeof(double), length(Nlen));
        if (top != true) {
            solve(f, u, Nlen, iters, 1, dx);
        } else {
            solve(f, u, Nlen, iters, w, dx);
        }

        restriction(f, u, fc, Nlen, dx*dx);
        multi(fc, uc, Nclen, dx*2, w, iters, false);
        interpolate(uc, d, Nlen, dx);
        for (int i = 0; i < length(Nlen); i++) {
            u[i] = u[i] + d[i];
        }
        if (top != true) {
            solve(f, u, Nlen, iters, 1, dx);
        } else {
            solve(f, u, Nlen, iters, w, dx);
        }
        free(fc);
        free(uc);
        free(d);
    }
}
