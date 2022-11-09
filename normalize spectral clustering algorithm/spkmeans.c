#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include "spkmeans.h"
#include "input.h"
#include "point.h"
#include "eigen.h"
#include "laplacian.h"
#include "debug.h"
#include "matrix.h"
#include "jacobi.h"
#include "kmeans.h"

int weighted_adjacency_matrix(const size_t n, double **weight_mat, double **points, const size_t dim)
{
    size_t i, j;
    point_t *p_points;
    p_points = malloc_points(n, dim);
    if (p_points == NULL)
    {
        return MALLOC_ERROR;
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < dim; j++)
        {
            p_points[i].elements[j] = points[i][j];
        }
    }

    create_weight_matrix(n, weight_mat, p_points, dim);

    free_points(n, p_points);
    return OK;
}

int diagonal_degree_matrix(const size_t n, double **d_mat, double **weigth_mat)
{
    return create_diagonal_degree_matrix(n, d_mat, weigth_mat);
}

int normalized_graph_laplacian(const size_t n, double **n_mat, double **w_mat, double **d_mat)
{
    double **l_mat;
    l_mat = (double **)malloc_matrix(n, n, sizeof(double));
    if (l_mat == NULL)
    {
        return MALLOC_ERROR;
    }
    create_laplacian_matrix(n, w_mat, d_mat, l_mat);

    return create_normalized_laplacian_matrix(n, l_mat, n_mat);
}

int calc_eigen_values_vectors(const size_t n, double **l_mat, double *values, double **vectors)
{
    eigen_t *eigens;
    size_t i, j;

    eigens = malloc_eigens(n);

    if (NULL == eigens)
    {
        return MALLOC_ERROR;
    }
    jacobi(n, l_mat, eigens);
    for (i = 0; i < n; i++)
    {
        values[i] = eigens[i].value;
        for (j = 0; j < n; j++)
        {
            vectors[i][j] = eigens[j].vector[i];
        }
    }

    free_eigens(n, eigens);
    return OK;
}

int kmeans(point_t *points, const size_t points_len, point_t *clusters, const size_t k, const size_t max_iter, const size_t dim, const float epsilon)
{
    int result;
    result = fit(points, points_len, clusters, k, max_iter, dim, epsilon);
    return result;
}

error_e calc_matrix(const size_t n, point_t *points, const size_t dim, goal_e goal, double **mat, size_t *k)
{
    error_e result;
    double **w_mat, **d_mat, **l_mat, **n_mat, **u_mat, **t_mat;
    eigen_t *eigens;
    size_t i;
    w_mat = (double **)malloc_matrix(n, n, sizeof(double));
    if (NULL == w_mat)
    {
        result = MALLOC_ERROR;
        goto end;
    }
    result = create_weight_matrix(n, w_mat, points, dim);
    if (result != OK)
    {
        goto w_cleanup;
    }

    if (WEIGHT_MATRIX == goal)
    {
        copy_matrix(n, w_mat, mat);
        goto w_cleanup;
    }

    d_mat = (double **)malloc_matrix(n, n, sizeof(double));
    if (NULL == d_mat)
    {
        result = MALLOC_ERROR;
        goto w_cleanup;
    }

    result = create_diagonal_degree_matrix(n, d_mat, w_mat);
    if (result != OK)
    {
        goto d_cleanup;
    }

    if (DIAGONAL_DEGREE_MATRIX == goal)
    {
        copy_matrix(n, d_mat, mat);
        goto d_cleanup;
    }

    l_mat = (double **)malloc_matrix(n, n, sizeof(double));
    if (NULL == l_mat)
    {
        result = MALLOC_ERROR;
        goto d_cleanup;
    }
    result = create_laplacian_matrix(n, w_mat, d_mat, l_mat);
    if (result != OK)
    {
        goto l_cleanup;
    }

    n_mat = (double **)malloc_matrix(n, n, sizeof(double));
    if (NULL == n_mat)
    {
        result = MALLOC_ERROR;
        goto l_cleanup;
    }

    result = create_normalized_laplacian_matrix(n, l_mat, n_mat);
    if (result != OK)
    {
        goto n_cleanup;
    }
    if (NORMALIZED_GRAPH_LAPLACIAN == goal)
    {
        copy_matrix(n, n_mat, mat);
        goto n_cleanup;
    }

    /* if (NORMALIZED_EIGEN_MATRIX == goal) */

    eigens = malloc_eigens(n);
    if (NULL == eigens)
    {
        result = MALLOC_ERROR;
        goto n_cleanup;
    }

    jacobi(n, n_mat, eigens);
    qsort(eigens, n, sizeof(eigen_t), compare_eigenvalues);

    if (0 == *k)
    {
        *k = find_eigengap_max(n, eigens);
    }

    u_mat = (double **)malloc_matrix(n, *k, sizeof(double));
    if (NULL == u_mat)
    {
        result = MALLOC_ERROR;
        goto eigen_vectors_cleanup;
    }
    build_matrix_from_eigens(n, *k, eigens, u_mat);
    t_mat = (double **)malloc_matrix(n, *k, sizeof(double));
    if (NULL == t_mat)
    {
        result = MALLOC_ERROR;
        goto u_cleanup;
    }

    normalize_matrix(n, *k, t_mat, u_mat);
    copy_matrix(n, t_mat, mat);

    free_matrix(n, (void **)t_mat);
u_cleanup:
    free_matrix(n, (void **)u_mat);
eigen_vectors_cleanup:
    for (i = 0; i < n; i++)
    {
        if (NULL != eigens[i].vector)
        {
            free(eigens[i].vector);
        }
    }
    free(eigens);

n_cleanup:
    free_matrix(n, (void **)n_mat);
l_cleanup:
    free_matrix(n, (void **)l_mat);
d_cleanup:
    free_matrix(n, (void **)d_mat);
w_cleanup:
    free_matrix(n, (void **)w_mat);
end:
    return result;
}

int create_eigen_matrix(const size_t n, double **l_mat, eigen_t *eigens)
{
    return jacobi(n, l_mat, eigens);
}

goal_e get_goal(char *goal_str)
{
    if (strcmp(goal_str, "wam") == 0)
    {
        return WEIGHT_MATRIX;
    }
    else if (strcmp(goal_str, "ddg") == 0)
    {
        return DIAGONAL_DEGREE_MATRIX;
    }
    else if (strcmp(goal_str, "lnorm") == 0)
    {
        return NORMALIZED_GRAPH_LAPLACIAN;
    }
    else if (strcmp(goal_str, "jacobi") == 0)
    {
        return JACOBI;
    }
    else if (strcmp(goal_str, "spk") == 0)
    {
        return NORMALIZED_EIGEN_MATRIX;
    }
    else
    {
        return UNKNOWN_GOAL;
    }
}

int main(int argc, char **argv)
{
    error_e result;
    goal_e goal;
    point_t *points;
    size_t n, i, dim, k;
    double **mat;
    eigen_t *eigens;

    k = 0;
    result = OK;
    if (argc != 3)
    {
        result = INVALID_INPUT;
        goto end;
    }
    goal = get_goal(argv[1]);

    if (UNKNOWN_GOAL == goal)
    {
        result = INVALID_INPUT;
        goto end;
    }
    n = get_lines_count(argv[2]);
    if (0 == n)
    {
        result = INVALID_INPUT;
        goto end;
    }
    dim = get_dimension(argv[2]);
    if (0 == dim)
    {
        result = INVALID_INPUT;
        goto end;
    }

    mat = (double **)malloc_matrix(n, n, sizeof(double));
    if (NULL == mat)
    {
        result = MALLOC_ERROR;
        goto end;
    }

    if (JACOBI != goal)
    {
        points = malloc_points(n, dim);
        if (NULL == points)
        {
            result = MALLOC_ERROR;
            goto points_cleanup;
        }
        result = read_points(argv[2], points, n, dim);
        if (OK != result)
        {
            result = INVALID_INPUT;
            goto points_cleanup;
        }

        result = calc_matrix(n, points, dim, goal, mat, &k);
        if (OK == result)
        {
            print_matrix(n, n, mat);
        }

    points_cleanup:
        free_points(n, points);
    }
    else
    {
        result = read_matrix(argv[2], mat, n);
        if (OK != result)
        {
            result = INVALID_INPUT;
            goto mat_cleanup;
        }
        eigens = malloc_eigens(n);
        if (NULL == eigens)
        {
            result = MALLOC_ERROR;
            goto mat_cleanup;
        }

        result = create_eigen_matrix(n, mat, eigens);
        if (OK == result)
        {
            print_eigen(n, eigens);
        }

        for (i = 0; i < n; i++)
        {
            free(eigens[i].vector);
        }
        free(eigens);
    }

mat_cleanup:
    free_matrix(n, (void **)mat);
end:
    if (INVALID_INPUT == result)
    {
        printf("Invalid Input\n");
    }
    else if (MALLOC_ERROR == result)
    {
        printf("An Error Has Occurred\n");
    }
    return result;
}
