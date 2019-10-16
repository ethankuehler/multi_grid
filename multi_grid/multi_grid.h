//
// Created by dcrush on 10/3/19.
//

#ifndef C_SOR_3D_MULTI_GRID_H
#define C_SOR_3D_MULTI_GRID_H

#include <stdbool.h>

typedef struct _N_len {
    int i;
    int j;
    int k;
} N_len;

typedef struct _funcs_args {
    void (*solve_top)(const double*, double*, N_len, double);
    void (*solve_coarse)(const double*, double*, N_len, double);
    void (*solve_base)(const double*, double*, N_len, double);
} funcs_args;

/*
 * This is an implantation of multi_grid method.
 * f is the right hand side of the equation.
 * u is the initial guess and output.
 * Nlen is the diemsions of both f and u.
 * dx is the step size
 * args are the solver functions.
 * top should always be set to true.
 */
void multi(const double* f, double* u, N_len Nlen, double dx, funcs_args arg, bool top);


#endif //C_SOR_3D_MULTI_GRID_H


