#ifndef DEBUG_H
#define DEBUG_H

#include "eigen.h"

void print_matrix(const size_t n, const size_t k, double **matrix);
void print_array(const size_t n, double *matrix);
void print_eigen(const size_t n, eigen_t *eigen);
#endif /* DEBUG_H */
