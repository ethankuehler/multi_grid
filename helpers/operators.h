#ifndef C_SOR_3D_OPERATORS_H
#define C_SOR_3D_OPERATORS_H

#include "N_len.h"

/*
 * function to find index of a vector from three coordinates.
 * i, j, k are the coordinates
 * Ni and Nj are the number of girds in the i'tj and j'th dimensions
 * returns the index of the array that would correspond to the given coordinates
 */
int loc(int i, int j, int k, N_len Nlen);

//reduce a m sized grid into a m-1 grid. m is that of the fine grid
void reduce(const double* f_in, double* f_out, N_len Nlen);

/*
 * calculates the residual of u
 * f is the right hand side of the dif eq
 * u is the approximate solution of the dif eq
 * r is the output
 * n is the diemsions
 * dxs is the step size squared
 */
void residual(const double* f, const double* u, double* r, N_len Nlen, double dxs);

/*
 * This function calculates the residual of u and then reducing the residual to an m-1 grid
 * f is the right hand side of the dif eq
 * u is the approximate solution
 * f_out is the output
 * m is that of the fine grid
 * dxs is the step size squared
 */
void restriction(const double* f, const double* u, double* f_out, N_len Nlen, double dxs);

/*
 * interpolates the function f into f_out where m is that of the coarse grid.
 * Nlen is that of the finer grid.
 * dx is the step size.
 */
void interpolate(const double* f, double* f_out, N_len Nlen, double dx);

#endif //C_SOR_3D_OPERATORS_H
