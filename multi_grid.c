#include <math.h>
#include <stdlib.h>
#include "solve_block.h"
#include "operators.h"

/* f is the right hand side of the equation
 * u is the output of the function and is the potential
 * m is the int that determents the dimensions of the input
 * dx is the step size
 */
void multi(double* f, double* u, int m, double dx, double w, int iters) {
    int N = (int) pow(2, m) + 1;
    int Nc = (int) pow(2, m - 1) + 1;
    if (m == 1) {
        int n = 1 + N + N*N;
        u[n] = (1.0/-6.0)*(f[n]*dx - (u[n - N*N] + u[n - N] + u[n + 1] + u[n - 1] + u[n + N] + u[n + N*N]));
    } else {
        double* fc = calloc(sizeof(double), Nc*Nc*Nc);
        double* uc = calloc(sizeof(double), Nc*Nc*Nc);
        double* d = calloc(sizeof(double), N*N*N);
        if (m != 8) {
            solve(f, u, N, N, N, iters, 1, dx);
        } else {
            solve(f, u, N, N, N, iters, w, dx);
        }
        restriction(f, u, fc, m, dx*dx);
        multi(fc, uc, m - 1, dx*2, w, iters);
        interpolate(uc, d, m - 1, dx);
        for (int i = 0; i < N*N*N; i++) {
            u[i] = u[i] + d[i];
        }
        if (m != 8) {
            solve(f, u, N, N, N, iters, 1, dx);
        } else {
            solve(f, u, N, N, N, iters, w, dx);
        }
        free(fc);
        free(uc);
        free(d);
    }
}
