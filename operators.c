#include "solve_block.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

//N is the dim of the input
void reduce(const double* f_in, double* f_out, int M) {
    int N = pow(2, M) + 1; // N of the fine grid
    int Nc = pow(2, M - 1) + 1; //N of the coarse grid
    //iterating though all interior points
    for(int i = 1; i < N-1; i++) {
        for(int j = 1; j < N-1; j ++) {
            for(int k = 1; k < N-1; k++) {
                //taking the
                f_out[loc(i/2,j/2,k/2, Nc, Nc)] =
                        (f_in[loc(i, j, k, N, N)] +
                        f_in[loc(i + 1, j, k, N, N)] +
                        f_in[loc(i - 1, j, k, N, N)] +
                        f_in[loc(i, j + 1, k, N, N)] +
                        f_in[loc(i, j - 1, k, N, N)] +
                        f_in[loc(i, j, k + 1, N, N)] +
                        f_in[loc(i, j, k - 1, N, N)])/7;
            }
        }
    }
}

void residual(const double* f, const double* u, double* r, int N, double dxs) {
    //because the outer points are fixed by the boundary conditions, their residual is zero.
    for (int i = 1; i < N - 1; i++) {
        for (int j = 1; j < N - 1; j++) {
            for (int k = 1; k < N - 1; k++) {
                int n = loc(i, j, k, N, N);
                r[n] =  f[n] - (u[n - N*N] + u[n - N] + u[n + 1] + -6*u[n] + u[n - 1] + u[n + N] + u[n + N*N])/dxs;
            }
        }
    }
}

void restriction(const double* f, const double* u, double* f_out, int m, double dxs) {
    //compute residual
    int N = (int) pow(2, m) + 1;//dimensions of fine grid
    double* r = calloc(sizeof(double), N*N*N); //freed at end of function

    //compute residual.
    residual(f, u, r, N, dxs);
    reduce(r, f_out, m);
    free(r);
}


//the m is that of the coarse grid
void interpolate(const double* f, double* f_out, int m, double dx) {
    int N = (int) pow(2, m + 1) + 1;//dimensions of fine grid
    int Nc = (int) pow(2, m) + 1;//dimensions of coarse grid

    //setting up the bottom back left corner.
    f_out[0] = f[0];

    //doing back and left edge of bottom and the left edge of the back wall.
    for (int i = 0; i < N - 2; i += 2) {
        int n = loc(0, 0, i, N, N);
        int nc = loc(0, 0, i/2, Nc, Nc);
        //back works
        f_out[n + 2] = f[nc + 1];
        f_out[n + 1] = f_out[n]*(1 - dx) + f_out[n + 2]*dx;

        //left works
        n = loc(0, i, 0, N, N);
        nc = loc(0, i/2, 0, Nc, Nc);
        f_out[n + 2*N] = f[nc + N];
        f_out[n + N] = f_out[n + 2*N]*(1 - dx) + f_out[n]*dx;

        //up works
        n = loc(i, 0, 0, N, N);
        nc = loc(i/2, 0, 0, Nc, Nc);
        f_out[n + 2*N*N] = f[nc + Nc*Nc];
        f_out[n + N*N] = f_out[n]*(1 - dx) + f_out[n + 2*N*N]*dx;
    }
    //doing bottom left and back wall
    for (int i = 0; i < N - 2; i += 2) {
        for (int j = 0; j < N - 2; j += 2) {

            //bottom tested to work
            int n = loc(0, i, j, N, N);
            int nc = loc(0, i/2, j/2, Nc, Nc);
            f_out[loc(0 , i + 2, j + 2, N, N)] =
                    f[loc(0, i/2 + 1, j/2 + 1, Nc, Nc)];

            f_out[loc(0 , i + 2, j + 1, N, N)] =
                    f_out[loc(0 , i + 2, j, N, N)]*(1 - dx) + f_out[loc(0 , i + 2, j + 2, N, N)]*dx;

            f_out[loc(0 , i + 1, j + 2, N, N)] =
                    f_out[loc(0 , i, j + 2, N, N)]*(1 - dx) + f_out[loc(0 , i + 2, j + 2, N, N)]*dx;

            f_out[loc(0 , i + 1, j + 1, N, N)] =
                    f_out[loc(0 , i + 1, j, N, N)]*(1 - dx) + f_out[loc(0 , i + 1, j + 2, N, N)]*dx;

            //back wall working
            f_out[loc(i + 2, 0, j + 2, N, N)] =
                    f[loc(i/2 + 1, 0, j/2 + 1, Nc, Nc)]; //top right

            f_out[loc(i + 2, 0, j + 1, N, N)] =
                    f_out[loc(i + 2, 0, j, N, N)]*(1 - dx) + f_out[loc(i + 2, 0, j + 2, N, N)]*dx; //top

            f_out[loc(i + 1, 0, j + 2, N, N)] =
                    f_out[loc(i , 0, j + 2, N, N)]*(1 - dx) + f_out[loc(i + 2, 0, j + 2, N, N)]*dx; //right

            f_out[loc(i + 1, 0, j + 1, N, N)] =
                    f_out[loc(i, 0, j + 2, N, N)]*(1 - dx) + f_out[loc(i + 2, 0, j + 2, N, N)]*dx;//middle

            //left wall //ending issue
            n = loc(i, j, 0, N, N);
            nc = loc(i/2, j/2, 0, Nc, Nc);
            f_out[loc(i + 2, j + 2, 0, N, N)] =
                    f[loc(i/2 + 1, j/2 + 1, 0, Nc, Nc)]; //foward top

            f_out[loc(i + 2, j + 1, 0, N, N)] =
                    f_out[loc(i + 2, j, 0, N, N)]*(1 - dx) + f_out[loc(i + 2, j + 2, 0, N, N)]*dx;//top

            f_out[loc(i + 1, j + 2, 0, N, N)] =
                    f_out[loc(i, j + 2, 0, N, N)]*(1 - dx) + f_out[loc(i + 2, j + 2, 0, N, N)]*dx; //foward

            f_out[loc(i + 1, j + 1, 0, N, N)] =
                    f_out[loc(i, j + 2, 0, N, N)]*(1 - dx) + f_out[loc(i + 2, j + 2, 0, N, N)]*dx;//middle
        }
    }

    //middle area and other walls.
    for (int i = 0; i < N - 2; i += 2) {
        for (int j = 0; j < N - 2; j += 2) {
            for (int k = 0; k < N - 2; k += 2) {
                int n = loc(i, j, k, N, N);
                int nc = loc(i/2, j/2, k/2, Nc, Nc);

                //last corner works
                f_out[loc(i + 2, j + 2, k + 2, N, N)] = f[loc(i/2 + 1, j/2 + 1, k/2 + 1, Nc, Nc)];

                //line parts
                //top right  working
                f_out[loc(i + 2, j + 1, k + 2, N, N)] =
                        f_out[loc(i + 2, j , k + 2, N, N)]*(1 - dx) + f_out[loc(i + 2, j + 2, k + 2, N, N)]*dx;

                //top forward working
                f_out[loc(i + 2, j + 2, k + 1, N, N)] =
                        f_out[loc(i + 2, j + 2, k, N, N)]*(1 - dx) + f_out[loc(i + 2, j + 2, k + 2, N, N)]*dx;

                //forward right working
                f_out[loc(i + 1, j + 2, k + 2, N, N)] =
                        f_out[loc(i, j + 2, k + 2, N, N)]*(1 - dx) + f_out[loc(i + 2, j + 2, k + 2, N, N)]*dx;

                //center of walls
                //right working
                f_out[loc(i + 1, j + 1, k + 2, N, N)] =
                        f_out[loc(i + 1, j, k + 2, N, N)]*(1 - dx) + f_out[loc(i + 1, j + 2, k + 2, N, N)]*dx;

                // forward working
                f_out[loc(i + 1, j + 2, k + 1, N, N)] =
                        f_out[loc(i + 1, j + 2, k, N, N)]*(1 - dx) + f_out[loc(i + 1, j + 2, k + 2, N, N)]*dx;

                // top
                f_out[loc(i + 2, j + 1, k + 1, N, N)] =
                        f_out[loc(i + 2, j + 1, k, N, N)]*(1 - dx) + f_out[loc(i + 2, j + 1, k + 2, N, N)]*dx;

                //middle point working
                f_out[loc(i + 1, j + 1, k + 1, N, N)] =
                        f_out[loc(i + 1, j + 1, k, N, N)]*(1 - dx) + f_out[loc(i + 1, j + 1, k + 2, N, N)]*dx;
            }
        }
    }
}
