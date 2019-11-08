
#include "multi_grid/operators.h"
#include <omp.h>


void method1(float* u, const float* f, N_len Nlen, float w, float dx) {
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    int i;
//#pragma omp parallel for num_threads(2) schedule(auto) shared(u, f, Nlen, w, dx) private(i)
    for (i = 1; i < Ni - 1; i++) {
        for (int j = 1; j < Nj - 1; j++) {
            for (int k = 1; k < Nk - 1; k++) {
                int n = loc(i, j, k, Nlen);
                u[n] = (1 - w)*u[n] +
                       (w/-6)*
                       (f[n]*dx - (u[n - Ni*Nj] + u[n - Nj] + u[n + 1] + u[n - 1] + u[n + Nj] + u[n + Ni*Nj]));
            }
        }
    }
}

void method2(float* u, const float* f, N_len Nlen, float w, float dx) {
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    for (int i = Ni - 2; i > 0; i--) {
        for (int j = 1; j < Nj - 1; j++) {
            for (int k = 1; k < Nk - 1; k++) {
                int n = loc(i, j, k, Nlen);
                u[n] = (1 - w)*u[n] +
                       (w/-6)*
                       (f[n]*dx - (u[n - Ni*Nj] + u[n - Nj] + u[n + 1] + u[n - 1] + u[n + Nj] + u[n + Ni*Nj]));
            }
        }
    }
}

void method3(float* u, const float* f, N_len Nlen, float w, float dx) {
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    for (int i = 1; i < Ni - 1; i++) {
        for (int j = 1; j < Nj - 1; j++) {
            for (int k = Nk - 2; k > 0; k--) {
                int n = loc(i, j, k, Nlen);
                u[n] = (1 - w)*u[n] +
                       (w/-6)*
                       (f[n]*dx - (u[n - Ni*Nj] + u[n - Nj] + u[n + 1] + u[n - 1] + u[n + Nj] + u[n + Ni*Nj]));
            }
        }
    }
}

void method4(float* u, const float* f, N_len Nlen, float w, float dx) {
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    for (int i = 1; i < Ni - 1; i++) {
        for (int j = Nj - 2; j > 0; j--) {
            for (int k = 1; k < Nk - 1; k++) {
                int n = loc(i, j, k, Nlen);
                u[n] = (1 - w)*u[n] +
                       (w/-6)*
                       (f[n]*dx - (u[n - Ni*Nj] + u[n - Nj] + u[n + 1] + u[n - 1] + u[n + Nj] + u[n + Ni*Nj]));
            }
        }
    }
}

void method5(float* u, const float* f, N_len Nlen, float w, float dx) {
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    for (int i = 1; i < Ni - 1; i++) {
        for (int j = Nj - 2; j > 0; j--) {
            for (int k = Nk - 2; k > 0; k--) {
                int n = loc(i, j, k, Nlen);
                u[n] = (1 - w)*u[n] +
                       (w/-6)*
                       (f[n]*dx - (u[n - Ni*Nj] + u[n - Nj] + u[n + 1] + u[n - 1] + u[n + Nj] + u[n + Ni*Nj]));
            }
        }
    }
}


void method6(float* u, const float* f, N_len Nlen, float w, float dx) {
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    for (int i = Ni - 2; i > 0; i--) {
        for (int j = Nj - 2; j > 0; j--) {
            for (int k = Nk - 2; k > 0; k--) {
                int n = loc(i, j, k, Nlen);
                u[n] = (1 - w)*u[n] +
                       (w/-6)*
                       (f[n]*dx - (u[n - Ni*Nj] + u[n - Nj] + u[n + 1] + u[n - 1] + u[n + Nj] + u[n + Ni*Nj]));
            }
        }
    }
}

void solve(const float* f, float* u, N_len Nlen, int iters, float w, float dx) {
    dx = dx*dx;
    for (int i = 0; i < iters; i++) {
        method1(u, f, Nlen, w, dx);
        //method2(u, f, Nlen, w, dx);
        //method3(u, f, Nlen, w, dx);
        //method4(u, f, Nlen, w, dx);
        //method5(u, f, Nlen, w, dx);
        //method6(u, f, Nlen, w, dx);
    }
}
