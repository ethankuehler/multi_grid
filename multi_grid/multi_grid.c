#include <stdlib.h>
#include <stdbool.h>
#include "operators.h"
#include "multi_grid.h"

/*
 * This fractionation implements multi grid v cycle recursively. First is checks to see if the gird can be reduced.
 * If it cant it will use the solve_base solver to get the solution to the coarsest gird and return. If it can it will
 * allocate needed memory. Then check to see which solver to use and then use it, restrict the current solution to the
 * coarser gird to get approximate problem. Afterwards it will call the multi recursively to solve approximate problem.
 * It will interpolate solution to approximate problem and add that to the current solution. Finally it will
 * again choose a solver, free memory and return.
 */
void multi(const float* f, float* u, N_len Nlen, float dx, funcs_args arg, bool top) {
    //checking to see if on coarsest
    if (!can_coarsen(Nlen)) {
        arg.solve_base(f, u, Nlen, dx);//getting solution to coarsest grid
    } else {
        N_len Nclen = coarsen(Nlen); //Nlen for coarser gird, hence the c.
        //Arrays to store coarser gird.
        float* fc = calloc(sizeof(float), length(Nclen));
        float* uc = calloc(sizeof(float), length(Nclen));
        //Array to sore interpolated solution.
        float* d = calloc(sizeof(float), length(Nlen));

        //checking to see if finest gird
        if (top != true) {
            arg.solve_coarse(f, u, Nlen, dx); //solver for not finest gird
        } else {
            arg.solve_top(f, u, Nlen, dx); //solver for the finest gird
        }

        restriction(f, u, fc, Nlen, dx*dx); //restricting current solution to get new problem.
        multi(fc, uc, Nclen, dx*2, arg, false); //recursively calling multi grid to get solution to new problem.
        interpolate(uc, d, Nlen, dx); //interpolating solution to new problem to get correction.

        //adding correction too solution
        for (int i = 0; i < length(Nlen); i++) {
            u[i] = u[i] + d[i];
        }

        //checking which solver to use again
        if (top != true) {
            arg.solve_coarse(f, u, Nlen, dx);
        } else {
            arg.solve_top(f, u, Nlen, dx);
        }

        free(fc);
        free(uc);
        free(d);
    }
}
