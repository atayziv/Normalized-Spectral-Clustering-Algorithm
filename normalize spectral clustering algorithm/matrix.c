/*This page groups all the operations on the matrices relevant to the project*/

#include <stdlib.h>
#include <math.h>
#include "matrix.h"

double row_norm(const size_t n, double **mat, const size_t row);

/*A function to allocate memory space for an n*k matrix*/
void **malloc_matrix(const size_t n, const size_t k, size_t elem_size)
{
    size_t i;
    void **matrix = malloc(n * sizeof(void *));
    if (NULL == matrix)
    {
        return NULL;
    }
    for (i = 0; i < n; i++)
    {
        matrix[i] = calloc(k, elem_size);
    }

    /* void *values = calloc(n * k, elem_size);
    for (i = 0; i < n; i++)
    {
        matrix[i] = values + i * n * elem_size;
    } */

    return matrix;
}

/*A function to free the memory space for a matrix*/
void free_matrix(const size_t n, void **matrix)
{
    size_t i;
    for (i = 0; i < n; i++)
    {
        free(matrix[i]);
    }
    free(matrix);
}

/*Boolean function to check if a matrix is diagonal*/
int is_diagonal(const size_t n, double **mat)
{
    size_t i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i != j && mat[i][j] != .0) /* If there is an off-diagonal element of the matrix that is different from zero, it means that the matrix is not diagonal */
            {
                return FALSE;
            }
        }
    }
    return TRUE;
}

/*A function to transpose on a matrix*/
int transpose(const size_t n, double **mat, double **transposed)
{
    size_t i, j;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            transposed[j][i] = mat[i][j];
        }
    }
    return 0;
}

/*A function to multiply matrices with maximum efficiency!*/
void multiply_mat(const size_t n, double **left_mat, double **right_mat, double **result)
{
    size_t i, j, k;
    int l_diagonal = is_diagonal(n, left_mat), r_diagonal = is_diagonal(n, right_mat);
    if (TRUE == l_diagonal && TRUE == r_diagonal)
    {
        for (i = 0; i < n; i++)
        {
            result[i][i] = left_mat[i][i] * right_mat[i][i]; /* If both matrices are diagonal, it is sufficient to multiply diagonal by diagonal */
        }
    }
    else if (TRUE == l_diagonal)
    {
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                result[i][j] = left_mat[i][i] * right_mat[i][j]; /* In case the left matrix is diagonal, it is enough to multiply its diagonal in the entire right matrix */
            }
        }
    }
    else if (TRUE == r_diagonal)
    {
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                result[i][j] = left_mat[i][j] * right_mat[j][j]; /* The same as the explanation for a diagonal left matrix only reversed (right with left reversed) */
            }
        }
    }
    else
    {
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                result[i][j] = .0;
                for (k = 0; k < n; k++)
                {
                    result[i][j] += left_mat[i][k] * right_mat[k][j]; /* Otherwise - the two matrices are not diagonal, we multiply in the naive way */
                }
            }
        }
    }
}

/*Matrix copy function*/
void copy_matrix(const size_t n, double **src, double **dst)
{
    size_t i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            dst[i][j] = src[i][j];
        }
    }
}

/*Function to create the unit matrix*/
void init_eye_matrix(const size_t n, double **mat)
{
    size_t i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            mat[i][j] = i == j ? 1.0 : 0.0;
        }
    }
}

/*Function to create the zero matrix*/
void init_zero_matrix(const size_t n, double **mat)
{
    size_t i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            mat[i][j] = 0.0;
        }
    }
}

/*function to normalize a matrix*/
void normalize_matrix(const size_t n, const size_t k, double **normalized, double **mat)
{
    size_t i, j;
    double r_sum;
    for (i = 0; i < n; i++)
    {
        r_sum = row_norm(k, mat, i);
        if (r_sum > .0)
        {
            for (j = 0; j < k; j++)
            {
                normalized[i][j] = mat[i][j] / r_sum; /* Each term is equal to the root of the sum of the terms in its same row */
            }
        }
    }
}

double row_norm(const size_t k, double **mat, const size_t row)
{
    size_t i;
    double sum = .0;
    for (i = 0; i < k; i++)
    {
        sum += pow(mat[row][i], 2.0);
    }
    return sqrt(sum);
}

/*A function to create a diagonal weight matrix*/
int create_diagonal_degree_matrix(const size_t n, double **d_mat, double **weigth_mat)
{
    size_t i, j;
    for (i = 0; i < n; i++)
    {
        d_mat[i][i] = .0;
        for (j = 0; j < n; j++)
        {
            /* if (i != j && adjency_matrix[i][j] != .0)*/
            {
                d_mat[i][i] += weigth_mat[i][j]; /* The sum of the terms in a row will enter the diagonal term in that row */
            }
        }
    }
    return 0;
}
