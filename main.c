#include <stdlib.h>
#include <stdbool.h>
#include "solve_block.h"
#include "multi_grid/multi_grid.h"
#include "tSolver/tSolver.h"
#include <omp.h>
#include "stuff.h"

void solve_top(const double* f, double* u, N_len Nlen, double dx) {
    tSolve(f, u, Nlen, 1, 1.9, dx);
}

void solve_coarse(const double* f, double* u, N_len Nlen, double dx) {
    tSolve(f, u, Nlen, 5, 1, dx);
}

void solve_base(const double* f, double* u, N_len Nlen, double dx) {
    solve(f, u, Nlen, 1, 1, dx);
}


int main() {
    int N = 513;
    N_len Nlen = (N_len){N, N, N};
    int iters = 5;
    double w = 1.2;
    double L = 20;
    double dx = L/N;
    double dens = 1;
    double R = 1;
    double a = 1;
    double b = 1;
    double c = 1;
    double* u = calloc(sizeof(double), N*N*N);
    double* u2 = calloc(sizeof(double), N*N*N);
    double* f = calloc(sizeof(double), N*N*N);



    //dose the initial conditions and stuff
    //its not really important what it dose
    inital(u, u2, f, dens, R, Nlen, L, dx, dx);


    //solve(f, u, (N_len){N, N, N}, 100, w, dx);

    //printf("number of iters: %d\n", iters*6);
    //solve(f, u, (N_len){N, N, N}, 100, w, dx);
    //save_gird("data2.txt", u, N*N*N);


    funcs_args arg = (funcs_args){solve_top, solve_coarse, solve_base};
    double start = omp_get_wtime();
    multi(f, u2, Nlen, dx, arg, true); //muti call
    double end = omp_get_wtime();
    //tSolve(f, u2, (N_len){N, N, N}, 50, w, dx);


    //dose all the output
    //save_gird("data3.txt", u2, N*N*N);
    data(u, u2, f, Nlen, start, end, dx, L, R, dens, a, b, c);

    return 0;
}