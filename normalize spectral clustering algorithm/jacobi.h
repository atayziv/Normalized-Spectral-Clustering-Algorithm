#ifndef JACOBI_H
#define JACOBI_H

#include <stdlib.h>
#include "eigen.h"

typedef struct mat_index_t
{
    size_t i;
    size_t j;
} mat_index_t;

int jacobi(const size_t n, double **mat, eigen_t *eigens);

#endif /* JACOBI_H */
