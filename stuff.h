//
// Created by dcrush on 10/29/19.
//

#ifndef C_SOR_3D_STUFF_H
#define C_SOR_3D_STUFF_H

#include "multi_grid/operators.h"
#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>

#define NUM 4

void func(float* u, int i, int j, int k, float L, int N, float m, float dx) {
    float x = pow((i*dx - L/2), 2);
    float y = pow((j*dx - L/2), 2);
    float z = pow((k*dx - L/2), 2);
    float rsqrd = x + y + z;
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

void save_line(const char* str, float* vec,int i, int j, N_len Nlen) {
    FILE* file = fopen(str, "w");
    int n = loc(i, j, 0, Nlen);
    for (int i = n; i < n + Nlen.k; i++) {
        fprintf(file, "%f ", vec[i]);
    }
    fprintf(file, "\n");
    fclose(file);
}

void save_gird(const char* str, float* vec, int length) {
    FILE* file = fopen(str, "w");

    for (int i = 0; i < length; i++) {
        fprintf(file, "%f ", vec[i]);
    }
    fprintf(file, "\n");
    fclose(file);
}

float L2(float* f, float* u, N_len Nlen, float dx) {
    float sum = 0;
    float* r = calloc(sizeof(float*), length(Nlen));
    residual(f, u, r, Nlen, dx*dx);
    int i;
#pragma omp parallel for num_threads(NUM) private(i)
    for(i = 0; i < length(Nlen); i++) {
        sum += pow(r[i],2);
    }
    free(r);
    return pow(sum, 1.0/2.0)/length(Nlen);
}

float avg_diff(const float* f1, const float* f2, N_len Nlen) {
    float sum = 0;
    int i;
#pragma omp parallel for num_threads(NUM) private(i)
    for(i = 0; i < length(Nlen); i++) {
        sum += fabs(f1[i] - f2[i]);
    }
    return sum/length(Nlen);
}

typedef struct _location_value {
    int i;
    int j;
    int k;
    float f;
}location_value;

location_value max_diff(const float* f1, const float* f2, N_len Nlen) {
    location_value max = (location_value){0, 0, 0, 0};
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    int i, j, k;
#pragma omp parallel for num_threads(NUM) private(i, j, k) collapse(3)
    for ( i = 0; i < Ni; i++) {
        for ( j = 0; j < Nj; j++) {
            for ( k = 0; k < Nk; k++) {
                int n = loc(i, j, k, Nlen);
                float one = f1[n];
                float two = f2[n];
                float dif = fabs(one - two);
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

void data(float* u, float* u2, float* f, N_len Nlen, float dx, float L, float R, float dens, float a,float  b, float c) {
    int N = Nlen.i;

    printf("error for multi :%lf\n", L2(f, u2, Nlen,dx));
    printf("error for SOR :%lf\n", L2(f, u, Nlen,dx));
    printf("avg diff :%f\n", avg_diff(u, u2, Nlen));

    location_value bruh = max_diff(u, u2, Nlen);
    //save_line("line.txt", u, N/2, N/2, Nlen);
    //save_line("line3.txt", u2, N/2, N/2, Nlen);

    printf("biggest diff was at {%i, %i, %i} and was %f\n", bruh.i, bruh.j, bruh.k, bruh.f);
    printf("%i\n", loc( bruh.i, bruh.j, bruh.k, Nlen));


    //save_gird("f.txt", f, N*N*N);
    save_params("params.txt", 9, (float) N, L, dx, R, dens, a, b, c, (float) N);

}

void inital(float* u, float* u2, float* f, float dens, float R, N_len Nlen, float L, float dx, float shift) {
    float x_mid, y_mid, z_mid;
    x_mid = y_mid = z_mid = L/2;

    int N = Nlen.i;
    //setting up source function
    int i, j, k;
#pragma omp parallel for num_threads(NUM) private(i, j, k) collapse(3)
    for ( i = 0; i < N; i++) {
        for ( j = 0; j < N; j++) {
            for ( k = 0; k < N; k++) {
                float x = (i*dx - (x_mid));
                float y = (j*dx - (y_mid));
                float z = (k*dx - (z_mid));
                float rsqrd = (x*x + y*y + z*z);
                int n = loc(i, j, k, (N_len){N, N, N});
                if (rsqrd < 1) {
                    f[n] = M_PI*4*dens;
                } else {
                    f[n] = 0;
                }
            }
        }
    }

    float m = 0;
#pragma omp parallel for num_threads(NUM) private(i, j, k) collapse(3)
    for ( i = 0; i < N; i++) {
        for ( j = 0; j < N; j++) {
            for ( k = 0; k < N; k++) {
                int n = loc(i,j,k, (N_len){N, N, N});
                if (f[n] != 0.0) {
                    m += dens*dx*dx*dx;
                }
            }
        }
    }

    printf("The mass is %f\n", m);
#pragma omp parallel for num_threads(NUM) private(i, j, k) collapse(3)
    for ( i = 0; i < N; i++) {
        for ( j = 0; j < N; j++) {
            for ( k = 0; k < N; k++) {
                float x = (i*dx - L/2);
                float y = (j*dx - L/2);
                float z = (k*dx - L/2 + shift);
                float rsqrd = x*x + y*y + z*z;
                int n = loc(i, j, k, (N_len){N, N, N});
                if (rsqrd < R*R) {
                    u2[n] = u[n] = -(m/(2*R*R*R))*(3*R*R - rsqrd);
                } else {
                    u2[n] = u[n] = -m/sqrt(rsqrd);
                }
            }
        }
    }
    printf("dx = %f\n", dx);
    //save_gird("data.txt", u, N*N*N);
    memcpy(u2, u, length((N_len){N, N, N})*sizeof(float));

#pragma omp parallel for num_threads(NUM) private(i, j, k) collapse(3)
    for ( i = 0; i < N; i++) {
        for ( j = 0; j < N; j++) {
            for ( k = 0; k < N; k++) {
                float x = (i*dx - L/2);
                float y = (j*dx - L/2);
                float z = (k*dx - L/2);
                float rsqrd = x*x + y*y + z*z;
                int n = loc(i, j, k, (N_len){N, N, N});
                if (rsqrd < R*R) {
                    u[n] = -(m/(2*R*R*R))*(3*R*R - rsqrd);
                } else {
                    u[n] = -m/sqrt(rsqrd);
                }
            }
        }
    }
}

#endif //C_SOR_3D_STUFF_H
