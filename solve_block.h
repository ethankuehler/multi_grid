//
// Created by ethan kueher on 9/26/19.
//

#ifndef C_SOR_3D_SOLVE_BLOCK_H
#define C_SOR_3D_SOLVE_BLOCK_H

#include "multi_grid/multi_grid.h"

/*
 * function to preform SOR
 * u is where the intial guess would go and also acts as the output.
 * f is the right hand side of the equation.
 * Ni, Nj, Nk are the dimensions of the problem.
 * iters determent how many iterations are done, the amount done is iters*6
 * w is the omega for sor
 * dx is the step size
 */
void solve(const double* f, double* u, N_len Nlen, int iters, double w, double dx) ;

#endif //C_SOR_3D_SOLVE_BLOCK_H
