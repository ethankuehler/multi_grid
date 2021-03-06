//
// Created by dcrush on 9/27/19.
//

#ifndef C_SOR_3D_OPERATORS_H
#define C_SOR_3D_OPERATORS_H

//reduce a m sized grid into a m-1 grid. m is that of the fine grid
void reduce(const double* f_in, double* f_out, int m);

/*
 * calculates the residual of u
 * f is the right hand side of the dif eq
 * u is the approximate solution of the dif eq
 * r is the output
 * n is the diemsions
 * dxs is the step size squared
 */
void residual(const double* f, const double* u, double* r, int N, double dxs);

/*
 * This function calculates the residual of u and then reducing the residual to an m-1 grid
 * f is the right hand side of the dif eq
 * u is the approximate solution
 * f_out is the output
 * m is that of the fine grid
 * dxs is the step size squared
 */
void restriction(const double* f, const double* u, double* f_out, int m, double dxs);

/*
 * interpolates the function f into f_out where m is that of the coarse grid.
 */
void interpolate(const double* f, double* f_out, int m, double dx);

#endif //C_SOR_3D_OPERATORS_H
