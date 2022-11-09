#ifndef EIGEN_H
#define EIGEN_H

typedef struct eigen_t
{
    double value;
    double *vector;
} eigen_t;

eigen_t *malloc_eigens(const size_t n);
void free_eigens(const size_t n, eigen_t *eigens);
int build_matrix_from_eigens(const size_t n, const size_t k, eigen_t *eigens, double **mat);
int compare_eigenvalues(const void *a, const void *b);
size_t find_eigengap_max(const size_t n, eigen_t *eigens);

#endif /* EIGEN_H */
