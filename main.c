#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include "solve_block.h"
#include "operators.h"
#include "multi_grid.h"


void func(double* u, int i, int j, int k, double L, int N, double m, double dx) {
    double x = pow((i*dx - L/2), 2);
    double y = pow((j*dx - L/2), 2);
    double z = pow((k*dx - L/2), 2);
    double rsqrd = x + y + z;
    int n = loc(i, j, k, N, N);
    u[n] = -m/sqrt(rsqrd);
}

void save_params(const char* str, int count, ...) {
    FILE* file = fopen(str, "w");
    va_list list;
    va_start(list, count);
    for (int i = 0; i < count; i++) {
        fprintf(file, "%f\n", va_arg(list, double));
    }
    va_end(list);
    fclose(file);
}

void save_gird(const char* str, double* vec, int length) {
    FILE* file = fopen(str, "w");

    for (int i = 0; i < length; i++) {
        fprintf(file, "%f ", vec[i]);
    }
    fprintf(file, "\n");
    fclose(file);
}


int main() {
    int M = 8;
    int N = (int) pow(2, M) + 1;
    int N2 = (int) pow(2, M - 1) + 1;
    int iters = 5;
    double w = 1.98;
    double L = 20;
    double dx = L/N;
    double dens = 1;
    double R = 1;
    double a = 1;
    double b = 1;
    double c = 1;
    double m = (4.0/3.0)*M_PI*a*b*c*dens;
    double* u = calloc(sizeof(double), N*N*N);
    double* u2 = calloc(sizeof(double), N*N*N);
    double* f = calloc(sizeof(double), N*N*N);

    double x_mid, y_mid, z_mid;
    x_mid = y_mid = z_mid = L/2;


    //setting up source function
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                double x = (i*dx - (x_mid));
                double y = (j*dx - (y_mid));
                double z = (k*dx - (z_mid - 3));
                double rsqrd = (1/(a*a))*x*x + (1/(b*b))*y*y + (1/(c*c))*z*z;
                int n = loc(i, j, k, N, N);
                if (rsqrd < 1) {
                    f[n] = M_PI*4*dens;
                } else {
                    f[n] = 0;
                }
            }
        }
    }

    //setting up source function
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                double x = (i*dx - (x_mid));
                double y = (j*dx - (y_mid));
                double z = (k*dx - (z_mid + 3));
                double rsqrd = (1/(a*a))*x*x + (1/(b*b))*y*y + (1/(c*c))*z*z;
                int n = loc(i, j, k, N, N);
                if (rsqrd < 1) {
                    f[n] = M_PI*4*dens;
                }
            }
        }
    }


    //setting up boundary conditions
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            func(u, i, j, 0, L, N, m, dx);
            func(u, i, j, N - 1, L, N, m, dx);
            func(u, i, 0, j, L, N, m, dx);
            func(u, i, N - 1, j, L, N, m, dx);
            func(u, 0, i, j, L, N, m, dx);
            func(u, N - 1, i, j, L, N, m, dx);
            func(u2, i, j, 0, L, N, m, dx);
            func(u2, i, j, N - 1, L, N, m, dx);
            func(u2, i, 0, j, L, N, m, dx);
            func(u2, i, N - 1, j, L, N, m, dx);
            func(u2, 0, i, j, L, N, m, dx);
            func(u2, N - 1, i, j, L, N, m, dx);
        }
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                double x = (i*dx - L/2);
                double y = (j*dx - L/2);
                double z = (k*dx - L/2);
                double rsqrd = x*x + y*y + z*z;
                int n = loc(i, j, k, N, N);
                if (rsqrd < R*R) {
                    u2[n] = u[n] = -(m/(2*R*R*R))*(3*R*R - rsqrd);
                } else {
                    u2[n] = u[n] = -m/sqrt(rsqrd);
                }
            }
        }
    }



    //printf("number of iters: %d\n", iters*6);
    solve(f, u, N, N, N, iters*2, w, dx);
    save_gird("data2.txt", u, N*N*N);

    multi(f, u2, M, dx, w, iters);
    save_gird("data3.txt", u2, N*N*N);

    //restriction(f, u, u2, M, dx*dx);
    //save_gird("data2.txt", u2, N2*N2*N2);
    //printf("%f", u2[0]);


    save_gird("f.txt", f, N*N*N);
    save_params("params.txt", 9, (double) N, L, dx, R, dens, a, b, c, (double) N);

    return 0;


}