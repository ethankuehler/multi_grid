
#ifndef C_SOR_3D_MULTI_GRID_H
#define C_SOR_3D_MULTI_GRID_H

#include <stdbool.h>


//This struct is used to store the number of girds in each dimensions for a gird.
typedef struct _N_len {
    int i;
    int j;
    int k;
} N_len;

/*
 * This struct is used to store solver used in multi grid.
 * solve_top is used on the finest grid.
 * solve_coarse is used on all the coarser gird excepted the coarsest.
 * solve_base is used on the coarsest gird.
 */
typedef struct _funcs_args {
    void (*solve_top)(const float*, float*, N_len, float);
    void (*solve_coarse)(const float*, float*, N_len, float);
    void (*solve_base)(const float*, float*, N_len, float);
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
void multi(const float* f, float* u, N_len Nlen, float dx, funcs_args arg, bool top);


#endif //C_SOR_3D_MULTI_GRID_H


