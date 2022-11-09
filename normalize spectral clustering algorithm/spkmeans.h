#ifndef SPKMEANS_H
#define SPKMEANS_H

#include <stdlib.h>
#include "point.h"

typedef enum goal_e
{
    UNKNOWN_GOAL = -1,
    WEIGHT_MATRIX = 0,
    DIAGONAL_DEGREE_MATRIX = 1,
    NORMALIZED_GRAPH_LAPLACIAN = 2,
    JACOBI = 3,
    NORMALIZED_EIGEN_MATRIX = 4
} goal_e;

typedef enum error_e
{
    OK = 0,
    MALLOC_ERROR = 1,
    INVALID_INPUT = 2
} error_e;

int weighted_adjacency_matrix(const size_t n, double **weight_mat, double **points, const size_t dim);
int diagonal_degree_matrix(const size_t n, double **d_mat, double **weigth_mat);
int normalized_graph_laplacian(const size_t n, double **n_mat, double **w_mat, double **d_mat);
int calc_eigen_values_vectors(const size_t n, double **l_mat, double *values, double **vectors);

error_e calc_matrix(const size_t n, point_t *points, const size_t dim, goal_e goal, double **mat, size_t *k);
int kmeans(point_t *points, const size_t points_len, point_t *clusters, const size_t k, const size_t max_iter, const size_t dim, const float epsilon);

#endif /* SPKMEANS_H */
