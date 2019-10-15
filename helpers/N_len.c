#include "N_len.h"
#include <assert.h>
#include <stdbool.h>

N_len coarsen(N_len Nlen) {
    int i = (Nlen.i+1)/2;
    int j = (Nlen.k+1)/2;
    int k = (Nlen.k+1)/2;
    assert(!(i < 3 || j < 3 || k < 3));
    return (N_len){i, j, k};
}

int length(N_len Nlen) {
    return Nlen.i * Nlen.j * Nlen.k;
}

char can_coarsen(N_len Nlen) {
    int i = (Nlen.i+1)/2;
    int j = (Nlen.k+1)/2;
    int k = (Nlen.k+1)/2;
    if (!(i < 3 || j < 3 || k < 3)) {
        return true;
    }
    return false;
}

