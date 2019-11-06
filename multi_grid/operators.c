#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include "multi_grid.h"

// On each length it dose (N+1)/2 to calculate the coarser gird Nlen
// If any of the lengths are less then 3 it will assert.
N_len coarsen(N_len Nlen) {
    int i = (Nlen.i + 1)/2;
    int j = (Nlen.k + 1)/2;
    int k = (Nlen.k + 1)/2;
    assert(!(i < 3 || j < 3 || k < 3));
    return (N_len) {i, j, k};
}

// calculates the lengths of a coarser gird and returns true if all lengths are not less then 3
bool can_coarsen(N_len Nlen) {
    int i = (Nlen.i + 1)/2;
    int j = (Nlen.k + 1)/2;
    int k = (Nlen.k + 1)/2;
    if (!(i < 3 || j < 3 || k < 3)) {
        return true;
    }
    return false;
}

//just multiples the lengths together to get length
int length(N_len Nlen) {
    return Nlen.i*Nlen.j*Nlen.k;
}

//needs to be check to see if correct
int loc(int i, int j, int k, N_len Nlen) {
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    assert(i*Ni*Nj + j*Nj + k < Ni*Nj*Nk);
    return i*Ni*Nj + j*Nj + k;
}


/*
 * This function reduced a fine grid to a coarser grid by iterating though the points on the fine grid that mach and at
 * each point take the avg of it and the the closest points and map that onto the coarser grid.
 * This function only iterates on the interior points of the  of the gird as the boundary conditions are fixed
 * and because this function is used on residuals they will always be zero at the boundary.
 */
void reduce(const double* f_in, double* f_out, N_len Nlen) {
    N_len Nclen = coarsen(Nlen);
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    int i, j, k;
#pragma omp parallel for num_threads(4) private(i, j, k) collapse(3)
    //iterating though all interior points, we move by 2 in order to get to points that mach on the coarser gird.
    for (i = 1; i < Ni - 1; i += 2) {
        for (j = 1; j < Nj - 1; j += 2) {
            for (k = 1; k < Nk - 1; k += 2) {
                //taking the avg of the closest points on fine and mapping to the coarse grid.
                f_out[loc(i/2, j/2, k/2, Nclen)] =
                        (f_in[loc(i, j, k, Nlen)] +
                         f_in[loc(i + 1, j, k, Nlen)] +
                         f_in[loc(i - 1, j, k, Nlen)] +
                         f_in[loc(i, j + 1, k, Nlen)] +
                         f_in[loc(i, j - 1, k, Nlen)] +
                         f_in[loc(i, j, k + 1, Nlen)] +
                         f_in[loc(i, j, k - 1, Nlen)])/7;
            }
        }
    }
}

/*
 * This function iterates thought the interior points calculates the residual from central difference method.
 * Only the interior points are done because boundary conditions are fixed and their residual is zero.
 * This function assumes that input array is filled with zeros on the boundary
 */
void residual(const double* f, const double* u, double* r, N_len Nlen, double dxs) {
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    //because the outer points are fixed by the boundary conditions, their residual is zero.
    int i, j, k;
#pragma omp parallel for num_threads(4) private(i, j, k) collapse(3)
    for (i = 1; i < Ni - 1; i++) {
        for (j = 1; j < Nj - 1; j++) {
            for (k = 1; k < Nk - 1; k++) {
                int n = loc(i, j, k, Nlen);
                r[n] = f[n] -
                       (u[n - Ni*Nj] + u[n - Nj] + u[n - 1] + -6*u[n] + u[n + 1] + u[n + Nj] + u[n + Ni*Nj])/dxs;
            }
        }
    }
}

/*
 * The restriction function is used to created the coarser problem for multi grid. it first computes the residual using
 * the residual function and then reduces that to the coarser gird.
 */
void restriction(const double* f, const double* u, double* f_out, N_len Nlen, double dxs) {
    double* r = calloc(sizeof(double), length(Nlen)); //freed at end of function
    residual(f, u, r, Nlen, dxs);//compute residual.
    reduce(r, f_out, Nlen);//reduce residuals
    free(r);
}


/*
 * This function is used by multi grid to take the solutions from the coarser problem and interpolate them to the fine
 * grid. The interpolation done trilinear interpolation. It iterates by "blocks" on the fine grid. The corners of the
 * block are points that mach on the fine gird and coarse gird. This box is made up of 6 walls, 12 line and a middle
 * points. If you imagine the box infront of you, the the back wall is the one closest to you while forward is opposite
 * of that, top is the highest and bottom is the lowest, left is the left most wall and right is opposite of that.
 * Top right is the line that connects the top and right wall and the naming is the same for the other lines.
 * We assume that the bottom, left and right wall and all lines and corners are calculated already.
 * We can make this assumption because as we iterate all plots theses parts will overlap with already done parts or
 * boundary walls. On the leftover points the trilinear interpolation is done.
 */
void interpolate(const double* f, double* f_out, N_len Nlen, double dx) {
    N_len Nclen = coarsen(Nlen);
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    //middle area and other walls.
    //we only need to do the interior areas as their should be zero
    int i, j, k;
#pragma omp parallel for num_threads(4) private(i, j, k) collapse(3)
    for (i = 0; i < Ni - 2; i += 2) {
        for (j = 0; j < Nj - 2; j += 2) {
            for (k = 0; k < Nk - 2; k += 2) {

                //last corner
                f_out[loc(i + 2, j + 2, k + 2, Nlen)] = f[loc(i/2 + 1, j/2 + 1, k/2 + 1, Nclen)];

                //line parts
                //top right
                f_out[loc(i + 2, j + 1, k + 2, Nlen)] =
                        f_out[loc(i + 2, j, k + 2, Nlen)]*(1 - dx) + f_out[loc(i + 2, j + 2, k + 2, Nlen)]*dx;

                //top forward
                f_out[loc(i + 2, j + 2, k + 1, Nlen)] =
                        f_out[loc(i + 2, j + 2, k, Nlen)]*(1 - dx) + f_out[loc(i + 2, j + 2, k + 2, Nlen)]*dx;

                //forward right
                f_out[loc(i + 1, j + 2, k + 2, Nlen)] =
                        f_out[loc(i, j + 2, k + 2, Nlen)]*(1 - dx) + f_out[loc(i + 2, j + 2, k + 2, Nlen)]*dx;

                //center of walls
                //right
                f_out[loc(i + 1, j + 1, k + 2, Nlen)] =
                        f_out[loc(i + 1, j, k + 2, Nlen)]*(1 - dx) + f_out[loc(i + 1, j + 2, k + 2, Nlen)]*dx;

                //forward
                f_out[loc(i + 1, j + 2, k + 1, Nlen)] =
                        f_out[loc(i + 1, j + 2, k, Nlen)]*(1 - dx) + f_out[loc(i + 1, j + 2, k + 2, Nlen)]*dx;

                // top
                f_out[loc(i + 2, j + 1, k + 1, Nlen)] =
                        f_out[loc(i + 2, j + 1, k, Nlen)]*(1 - dx) + f_out[loc(i + 2, j + 1, k + 2, Nlen)]*dx;

                //middle point
                f_out[loc(i + 1, j + 1, k + 1, Nlen)] =
                        f_out[loc(i + 1, j + 1, k, Nlen)]*(1 - dx) + f_out[loc(i + 1, j + 1, k + 2, Nlen)]*dx;
            }
        }
    }
}
