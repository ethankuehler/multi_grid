#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <string.h>
#include "solve_block.h"
#include "multi_grid/operators.h"
#include "multi_grid/multi_grid.h"
#include <time.h>

void func(double* u, int i, int j, int k, double L, int N, double m, double dx) {
    double x = pow((i*dx - L/2), 2);
    double y = pow((j*dx - L/2), 2);
    double z = pow((k*dx - L/2), 2);
    double rsqrd = x + y + z;
    int n = loc(i, j, k, (N_len){N, N, N});
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

void save_line(const char* str, double* vec,int i, int j, N_len Nlen) {
    FILE* file = fopen(str, "w");
    int n = loc(i, j, 0, Nlen);
    for (int i = n; i < n + Nlen.k; i++) {
        fprintf(file, "%f ", vec[i]);
    }
    fprintf(file, "\n");
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

double L2(double* f, double* u, N_len Nlen, double dx) {
    double sum = 0;
    double* r = calloc(sizeof(double*), length(Nlen));
    residual(f, u, r, Nlen, dx*dx);

    for(int i = 0; i < length(Nlen); i++) {
        sum += pow(r[i],2);
    }
    free(r);
    return pow(sum, 1.0/2.0)/length(Nlen);
}

double avg_diff(const double* f1, const double* f2, N_len Nlen) {
    double sum = 0;
    for(int i = 0; i < length(Nlen); i++) {
        sum += fabs(f1[i] - f2[i]);
    }
    return sum/length(Nlen);
}

typedef struct _location_value {
    int i;
    int j;
    int k;
    double f;
}location_value;

location_value max_diff(const double* f1, const double* f2, N_len Nlen) {
    location_value max = (location_value){0, 0, 0, 0};
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    for (int i = 0; i < Ni; i++) {
        for (int j = 0; j < Nj; j++) {
            for (int k = 0; k < Nk; k++) {
                int n = loc(i, j, k, Nlen);
                double one = f1[n];
                double two = f2[n];
                double dif = fabs(one - two);
                if (max.f < dif){
                    max.f = dif;
                    max.i = i;
                    max.j = j;
                    max.k = k;
                }
            }
        }
    }
    return max;
}

void solve_top(const double* f, double* u, N_len Nlen, double dx) {
    solve(f, u, Nlen, 5, 1.9, dx);
}

void solve_coarse(const double* f, double* u, N_len Nlen, double dx) {
    solve(f, u, Nlen, 1, 1, dx);
}

void solve_base(const double* f, double* u, N_len Nlen, double dx) {
    solve(f, u, Nlen, 1, 1, dx);
}



int main() {
    int N = 701;
    int iters = 5;
    double w = 1.9;
    double L = 20;
    double dx = L/N;
    double dens = 1;
    double R = 1;
    double a = 1;
    double b = 1;
    double c = 1;
    //double* u = calloc(sizeof(double), N*N*N);
    double* u2 = calloc(sizeof(double), N*N*N);
    double* f = calloc(sizeof(double), N*N*N);

    /*
    double x_mid, y_mid, z_mid;
    x_mid = y_mid = z_mid = L/2;


    //setting up source function
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                double x = (i*dx - (x_mid));
                double y = (j*dx - (y_mid));
                double z = (k*dx - (z_mid));
                double rsqrd = (1/(a*a))*x*x + (1/(b*b))*y*y + (1/(c*c))*z*z;
                int n = loc(i, j, k, (N_len){N, N, N});
                if (rsqrd < 1) {
                    f[n] = M_PI*4*dens;
                } else {
                    f[n] = 0;
                }
            }
        }
    }

    double m = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                int n = loc(i,j,k, (N_len){N, N, N});
                if (f[n] != 0.0) {
                    m += dens*dx*dx*dx;
                }
            }
        }
    }

    printf("The mass is %f\n", m);

    //setting up boundary conditions
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            func(u, i, j, 0, L, N, m, dx);
            func(u, i, j, N - 1, L, N, m, dx);
            func(u, i, 0, j, L, N, m, dx);
            func(u, i, N - 1, j, L, N, m, dx);
            func(u, 0, i, j, L, N, m, dx);
            func(u, N - 1, i, j, L, N, m, dx);
        }
    }

    for (int i = 1; i < N-1; i++) {
        for (int j = 1; j < N-1; j++) {
            for (int k = 1; k < N-1; k++) {
                double x = (i*dx - L/2);
                double y = (j*dx - L/2);
                double z = (k*dx - L/2 + dx);
                double rsqrd = x*x + y*y + z*z;
                int n = loc(i, j, k, (N_len){N, N, N});
                if (rsqrd < R*R) {
                    u[n] = -(m/(2*R*R*R))*(3*R*R - rsqrd);
                } else {
                    u[n] = -m/sqrt(rsqrd);
                }
            }
        }
    }


    printf("dx = %f\n", dx);
    //save_gird("data.txt", u, N*N*N);
    memcpy(u2, u, length((N_len){N, N, N})*sizeof(double));


    for (int i = 1; i < N-1; i++) {
        for (int j = 1; j < N-1; j++) {
            for (int k = 1; k < N-1; k++) {
                double x = (i*dx - L/2);
                double y = (j*dx - L/2);
                double z = (k*dx - L/2);
                double rsqrd = x*x + y*y + z*z;
                int n = loc(i, j, k, (N_len){N, N, N});
                if (rsqrd < R*R) {
                    u[n] = -(m/(2*R*R*R))*(3*R*R - rsqrd);
                } else {
                    u[n] = -m/sqrt(rsqrd);
                }
            }
        }
    }
    //solve(f, u, (N_len){N, N, N}, 100, w, dx);

    //printf("number of iters: %d\n", iters*6);
    //solve(f, u, (N_len){N, N, N}, 100, w, dx);
    //save_gird("data2.txt", u, N*N*N);
    */
    N_len Nlen = (N_len){N, N, N};
    funcs_args arg = (funcs_args){solve_top, solve_coarse, solve_base};
    multi(f, u2, Nlen, dx, arg, true);
    printf("%f\n", u2[1000]);
    /*
    //save_gird("data3.txt", u2, N*N*N);


    printf("error for multi :%lf\n", L2(f, u2, Nlen,dx));
    printf("error for SOR :%lf\n", L2(f, u, Nlen,dx));
    printf("avg diff :%f\n", avg_diff(u, u2, Nlen));

    location_value bruh = max_diff(u, u2, Nlen);
    save_line("line.txt", u, N/2, N/2, Nlen);
    save_line("line3.txt", u2, N/2, N/2, Nlen);

    printf("biggest diff was at {%i, %i, %i} and was %f\n", bruh.i, bruh.j, bruh.k, bruh.f);
    printf("%i\n", loc( bruh.i, bruh.j, bruh.k, Nlen));


    //save_gird("f.txt", f, N*N*N);
    save_params("params.txt", 9, (double) N, L, dx, R, dens, a, b, c, (double) N);

    */
    return 0;


}