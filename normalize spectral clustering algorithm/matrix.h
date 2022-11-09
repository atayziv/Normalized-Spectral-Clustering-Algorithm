#ifndef MATRIX_H
#define MATRIX_H

#define TRUE 1
#define FALSE 0

void **malloc_matrix(const size_t n, const size_t k, size_t elem_size);
void free_matrix(const size_t n, void **matrix);

void multiply_mat(const size_t n, double **left_mat, double **right_mat, double **result);
int is_diagonal(const size_t n, double **mat);

int transpose(const size_t n, double **mat, double **transposed);
void copy_matrix(const size_t n, double **src, double **dst);
void init_eye_matrix(const size_t n, double **mat);
void init_zero_matrix(const size_t n, double **mat);

void normalize_matrix(const size_t n, const size_t k, double **mat, double **normalized);
int create_diagonal_degree_matrix(const size_t n, double **d_mat, double **weigth_mat);

#endif /* MATRIX_H */
