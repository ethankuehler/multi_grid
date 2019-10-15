//
// Created by dcrush on 10/14/19.
//

#ifndef C_SOR_3D_N_LEN_H
#define C_SOR_3D_N_LEN_H

typedef struct _N_len {
    int i;
    int j;
    int k;
} N_len;

N_len coarsen(N_len Nlen);
N_len refine(N_len Nlen);
int length(N_len Nlen);
char can_coarsen(N_len Nlen);

#endif //C_SOR_3D_N_LEN_H
