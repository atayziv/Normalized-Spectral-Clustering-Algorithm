#ifndef LAPLACIAN_H
#define LAPLACIAN_H

int create_laplacian_matrix(const size_t n, double **w_mat, double **d_mat, double **l_mat);
int create_normalized_laplacian_matrix(const size_t n, double **l_mat, double **n_mat);

#endif /* LAPLACIAN_H */
