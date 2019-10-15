#include <pthread.h>
#include "solve_thread.h"
#include "operators.h"
#include <stdbool.h>
#include <stdio.h>

typedef struct _arg_struct {
    const double* f;
    double* u;
    N_len Nlen;
    int i;
    double w;
    double dxs;
} arg_struct;

void solve_plane(const double* f, double* u, N_len Nlen, int si, double w, double dx) {
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    for(int i = si; i < Ni; i += 8) {
        for (int j = 1; j < Nj - 1; j++) {
            for (int k = 1; k < Nk - 1; k++) {
                int n = loc(i, j, k, Nlen);
                u[n] = (1 - w)*u[n] +
                       (w/-6)*(f[n]*dx - (u[n - Ni*Nj] + u[n - Nj] + u[n + 1] + u[n - 1] + u[n + Nj] + u[n + Ni*Nj]));
            }
        }
    }
}

void* call_solve_plane(void* args_in) {
    arg_struct* arg = (arg_struct*)args_in;
    double* u = arg->u;
    const double* f = arg->f;
    N_len Nlen = arg->Nlen;
    int i = arg->i;
    double w = arg->w;
    double dxs = arg->dxs;
    solve_plane(f, u, Nlen, i, w, dxs);
    return NULL;
}

void solve_threaded(const double* f, double* u, N_len Nlen, int iters, double w, double dx) {
    double dxs = dx*dx;
    pthread_t threads[4];
    for(int _ = 0; _ < iters; _++) {

        //red
        //populate threads
        int j = 0;
        for(int i = 0; i < 8; i+=2){
            arg_struct arg = (arg_struct){f, u, Nlen, i, w, dxs};
            arg_struct* ptr = &arg;
            pthread_create(&threads[j], NULL, call_solve_plane, (void *)ptr);
            j++;
        }
        for(int i = 0; i < 4; i++) {
            pthread_join(threads[i], NULL);
        }
        //black
        j = 0;
        for(int i = 1; i < 9; i+=2){
            arg_struct arg = (arg_struct){f, u, Nlen, i, w, dxs};
            arg_struct* ptr = &arg;
            pthread_create(&threads[j], NULL, call_solve_plane, (void *)ptr);
            j++;
        }

        for(int i = 0; i < 4; i++) {
            pthread_join(threads[i], NULL);
        }

    }
}

/*
void solve_plane2(double* u, const double* f, N_len Nlen, int ii, int ie, double w, double dxs) {
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    for (int i = ii; i < ie - 1; i++) {
        for (int j = 1; j < Nj - 1; j++) {
            for (int k = 1; k < Nk - 1; k++) {
                int n = loc(i, j, k, Nlen);
                u[n] = (1 - w)*u[n] +
                       (w/-6)*
                       (f[n]*dxs - (u[n - Ni*Nj] + u[n - Nj] + u[n + 1] + u[n - 1] + u[n + Nj] + u[n + Ni*Nj]));
            }
        }
    }
}


typedef struct _arg_struct2 {
    const double* f;
    double* u;
    N_len Nlen;
    int ii;
    int ie;
    double w;
    double dxs;
} arg_struct2;

void* call_solve_plane2(void* args_in) {
    arg_struct2* arg = (arg_struct2*)args_in;
    double* u = arg->u;
    const double* f = arg->f;
    N_len Nlen = arg->Nlen;
    int ii = arg->ii;
    int ie = arg->ie;
    double w = arg->w;
    double dxs = arg->dxs;
    solve_plane2(u, f, Nlen, ii, ie, w, dxs);
    return NULL;
}

void solve_threaded(const double* f, double* u, N_len Nlen, int iters, double w, double dx) {
    double dxs = dx*dx;
    pthread_t threads[4];
    int quarter = Nlen.i/4;
    for(int _ = 0; _ < iters; _++) {
        for(int i = 0; i < 4; i++) {
            arg_struct2 arg = (arg_struct2){f, u, Nlen, quarter*i,quarter*i  +  quarter, w, dxs};
            arg_struct2* ptr = &arg;
            pthread_create(&threads[i], NULL, call_solve_plane2, (void *)ptr);
        }
    }
    for(int i = 0; i < 4; i++) {
        pthread_join(threads[i], NULL);
    }
    for(int i = 1; i < 4; i++){
        solve_plane(u, f, Nlen, i*quarter - 1, w, dxs);
    }
}
*/
