#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "multi_grid.h"
#include "solve_block.h"
#include "operators.h"

void save_gird2(const char* str, double* vec, int length) {
    FILE* file = fopen(str, "w");

    for (int i = 0; i < length; i++) {
        fprintf(file, "%f ", vec[i]);
    }
    fprintf(file, "\n");
    fclose(file);
}
//f is the right hand side of the equation
//u is the output of the function and is the potential
//m is the int that determents the dimensions of the input
//dx is the step size
void multi_old(double* f, double* u, int m, int N, double dx) {
    if (m == 1) { // if m is one then exact solution is calculated
        int n = 1 + N + N*N;
        u[n] = (1.0/-6.0)*(f[n]*dx - (u[n - N*N] + u[n - N] + u[n + 1] + u[n - 1] + u[n + N] + u[n + N*N]));
    } else { //other wise approximation is done.
        double w = 1.98;//omega for finest grid
        double w2 = 1.0;//omega for all others
        int iter = 1;//number of iterations for finest gird
        int iter2 = 20;// number of iterations for all others.
        if( m != 8) { //doing SOR with inputs depending on if finest grid or not
            solve(f, u, N, N, N, iter2, w2, dx);
        } else {
            solve(f, u, N, N, N, iter, w, dx);
        }
        //setting up data needed for m-1 grid
        int N2 = (int)pow(2, m-1) + 1;
        double* u2 = calloc(sizeof(double), N2*N2*N2);
        double* f2 = calloc(sizeof(double), N2*N2*N2);
        double* d = calloc(sizeof(double), N*N*N);

        restriction(f, u, f2, m, dx*dx);//restricting the current solution to the coarser problem
        multi_old(f2, u2, m - 1, N2, dx * 2);//recursively calling multi function to solve the coarser problem.
        interpolate(u2, d, m - 1, dx);//interpolating solution of coarser problem back to finer gird.

        for(int i = 0; i < N*N*N; i++) { // subtracting coarser solution from fine solution
            u[i] -= d[i];
        }
        //freeing unneeded memory
        free(u2);
        free(f2);
        free(d);

        if(m != 8) { //doing SOR with inputs depending on if finest grid or not
            solve(f, u, N, N, N, iter2, w2, dx);
        } else {
            solve(f, u, N, N, N, iter, w, dx);
        }
    }
}

void multi(double* f, double* u, int m, double dx, double w, int iters) {
    int N = pow(2, m) + 1;
    int Nc = pow(2, m - 1) + 1;
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
        multi(fc, uc, m - 1, dx * 2, w, iters);
        interpolate(uc, d, m - 1, dx);
        for(int i = 0; i < N*N*N; i++) {
            u[i] = u[i] + d[i];
        }
        if( m == 8){
            //save_gird2("data3.txt", d, N *N *N);
        }

        if (m != 8) {
            solve(f, u, N, N, N, iters, 1, dx);
        } else {
            //solve(f, u, N, N, N, iters, w, dx);
        }
        free(fc);
        free(uc);
        free(d);
    }
}
