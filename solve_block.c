//
// Created by dcrush on 9/26/19.
//
#include "operators.h"


void method1(double* u, const double* f, N_len Nlen, double w, double dx) {
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    for (int i = 1; i < Ni - 1; i++) {
        for (int j = 1; j < Nj - 1; j++) {
            for (int k = 1; k < Nk - 1; k++) {
                int n = loc(i, j, k, Ni, Nj);
                u[n] = (1 - w)*u[n] +
                       (w/-6)*
                       (f[n]*dx - (u[n - Ni*Nj] + u[n - Nj] + u[n + 1] + u[n - 1] + u[n + Nj] + u[n + Ni*Nj]));
            }
        }
    }
}

void method2(double* u, const double* f, N_len Nlen, double w, double dx) {
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    for (int i = Ni - 2; i > 0; i--) {
        for (int j = 1; j < Nj - 1; j++) {
            for (int k = 1; k < Nk - 1; k++) {
                int n = loc(i, j, k, Ni, Nj);
                u[n] = (1 - w)*u[n] +
                       (w/-6)*
                       (f[n]*dx - (u[n - Ni*Nj] + u[n - Nj] + u[n + 1] + u[n - 1] + u[n + Nj] + u[n + Ni*Nj]));
            }
        }
    }
}

void method3(double* u, const double* f, N_len Nlen, double w, double dx) {
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    for (int i = 1; i < Ni - 1; i++) {
        for (int j = 1; j < Nj - 1; j++) {
            for (int k = Nk - 2; k > 0; k--) {
                int n = loc(i, j, k, Ni, Nj);
                u[n] = (1 - w)*u[n] +
                       (w/-6)*
                       (f[n]*dx - (u[n - Ni*Nj] + u[n - Nj] + u[n + 1] + u[n - 1] + u[n + Nj] + u[n + Ni*Nj]));
            }
        }
    }
}

void method4(double* u, const double* f, N_len Nlen, double w, double dx) {
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    for (int i = 1; i < Ni - 1; i++) {
        for (int j = Nj - 2; j > 0; j--) {
            for (int k = 1; k < Nk - 1; k++) {
                int n = loc(i, j, k, Ni, Nj);
                u[n] = (1 - w)*u[n] +
                       (w/-6)*
                       (f[n]*dx - (u[n - Ni*Nj] + u[n - Nj] + u[n + 1] + u[n - 1] + u[n + Nj] + u[n + Ni*Nj]));
            }
        }
    }
}

void method5(double* u, const double* f, N_len Nlen, double w, double dx) {
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    for (int i = 1; i < Ni - 1; i++) {
        for (int j = Nj - 2; j > 0; j--) {
            for (int k = Nk - 2; k > 0; k--) {
                int n = loc(i, j, k, Ni, Nj);
                u[n] = (1 - w)*u[n] +
                       (w/-6)*
                       (f[n]*dx - (u[n - Ni*Nj] + u[n - Nj] + u[n + 1] + u[n - 1] + u[n + Nj] + u[n + Ni*Nj]));
            }
        }
    }
}


void method6(double* u, const double* f, N_len Nlen, double w, double dx) {
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    for (int i = Ni - 2; i > 0; i--) {
        for (int j = Nj - 2; j > 0; j--) {
            for (int k = Nk - 2; k > 0; k--) {
                int n = loc(i, j, k, Ni, Nj);
                u[n] = (1 - w)*u[n] +
                       (w/-6)*
                       (f[n]*dx - (u[n - Ni*Nj] + u[n - Nj] + u[n + 1] + u[n - 1] + u[n + Nj] + u[n + Ni*Nj]));
            }
        }
    }
}

void solve(const double* f, double* u, N_len Nlen, int iters, double w, double dx) {
    dx = dx*dx;
    for (int i = 0; i < iters; i++) {
        method1(u, f, Nlen, w, dx);
        method2(u, f, Nlen, w, dx);
        method3(u, f, Nlen, w, dx);
        method4(u, f, Nlen, w, dx);
        method5(u, f, Nlen, w, dx);
        method6(u, f, Nlen, w, dx);
    }
}
