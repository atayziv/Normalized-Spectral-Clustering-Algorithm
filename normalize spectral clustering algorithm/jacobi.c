#include <math.h>
#include <float.h>
#include <stdio.h>
#include "jacobi.h"
#include "matrix.h"
#include "debug.h"
#include "eigen.h"

#define MAX_ITERATIONS 100
#define EPSILON 0.00001

/*Functions that return the values of θ, t, c, s according to the format in section 1.2.1-4*/
double get_tetha(double **mat, const mat_index_t index)
{
    return ((mat[index.j][index.j] - mat[index.i][index.i]) / (2.0 * mat[index.i][index.j]));
}
double get_t(const double tetha)
{
    double sign = (tetha >= .0 ? 1.0 : -1.0);
    return sign / (fabs(tetha) + sqrt(pow(tetha, 2.0) + 1));
}
double get_c(const double t)
{
    return 1.0 / sqrt(pow(t, 2.0) + 1.0);
}
double get_s(const double t, const double c)
{
    return t * c;
}

/*Function to initialize a rotation matrix cf. 1.2.1 -2*/
int init_rotation_matrix(const size_t n, double **mat, const double c, const double s, mat_index_t index)
{
    init_eye_matrix(n, mat);
    mat[index.i][index.i] = mat[index.j][index.j] = c;
    mat[index.i][index.j] = s;
    mat[index.j][index.i] = -1.0 * s;
    return 0;
}

/*Function to find the largest number off the diagonal of the matrix, for the pivot in section 1.2.1-3*/
mat_index_t find_max_off_diagonal(const size_t n, double **mat)
{
    double max = -DBL_MAX;
    size_t i, j;
    mat_index_t max_index;
    for (i = 0; i < n; i++)
    {
        for (j = i; j < n; j++)
        {
            if (i != j)
            {
                if (fabs(mat[i][j]) > max)
                {
                    max = fabs(mat[i][j]);
                    max_index.i = i;
                    max_index.j = j;
                }
            }
        }
    }
    return max_index;
}

/*Transform the matrix A to A' by Relation between A and A' (1.2.1 - 6) description*/
void transform_rotation(const size_t n, double **mat, double **result, const size_t i, const size_t j, const double s, const double c)
{
    size_t r;
    copy_matrix(n, mat, result);
    for (r = 0; r < n; r++)
    {
        if (r != i && r != j)
        {
            result[r][i] = result[i][r] = c * mat[r][i] - s * mat[r][j];
            result[r][j] = result[j][r] = c * mat[r][j] + s * mat[r][i];
        }
    }
    result[i][i] = (pow(c, 2.0) * mat[i][i]) + (pow(s, 2.0) * mat[j][j]) - (2.0 * s * c * mat[i][j]);
    result[j][j] = (pow(s, 2.0) * mat[i][i]) + (pow(c, 2.0) * mat[j][j]) + (2.0 * s * c * mat[i][j]);
    result[i][j] = result[j][i] = 0.0;
}

/*function for sum of squares of all off-diagonal elements of A and A' respectively as described in 1.2.1-5*/
double square_off_diagonal(const size_t n, double **mat)
{
    size_t i, j;
    double result = 0.0;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i != j)
            {
                result += pow(mat[i][j], 2.0);
            }
        }
    }
    return result;
}

/*Jacobian algorithm as shown in section 1.2.1*/
int jacobi(const size_t n, double **mat, eigen_t *eigens)
{
    int result;
    mat_index_t index;
    double **mat_cpy, **mat_tag, **rotation_mat, **vectors, **e_vectors, **temp;
    double tetha, t, c, s;
    double a_off_diag, a_tag_off_diag, convergence;
    double i_temp, j_temp;
    size_t i, j, l, iter;
    result = 0;
    mat_cpy = (double **)malloc_matrix(n, n, sizeof(double));
    if (NULL == mat_cpy)
    {
        result = 1;
        goto end;
    }
    mat_tag = (double **)malloc_matrix(n, n, sizeof(double));
    if (NULL == mat_tag)
    {
        result = 1;
        goto mat_cpy_cleanup;
    }
    rotation_mat = (double **)malloc_matrix(n, n, sizeof(double));
    if (NULL == rotation_mat)
    {
        result = 1;
        goto mat_tag_cleanup;
    }
    vectors = (double **)malloc_matrix(n, n, sizeof(double));
    if (NULL == vectors)
    {
        result = 1;
        goto rotation_mat_cleanup;
    }
    e_vectors = (double **)malloc_matrix(n, n, sizeof(double));
    if (NULL == e_vectors)
    {
        result = 1;
        goto vectors_cleanup;
    }
    init_eye_matrix(n, vectors); /* For the first iteration, we will set vectors to be the unit matrix */
    copy_matrix(n, mat, mat_cpy);

    iter = 0;
    a_tag_off_diag = 0.0;
    convergence = a_off_diag = square_off_diagonal(n, mat_cpy);
    while (convergence > EPSILON && iter < MAX_ITERATIONS) /* Algorithm stopping conditions as shown in section 1.2.1 -5 */
    {

        index = find_max_off_diagonal(n, mat_cpy); /* Extracting i & j (pivot indexes) */
        i = index.i;
        j = index.j;
        tetha = get_tetha(mat_cpy, index);
        t = get_t(tetha);
        c = get_c(t);
        s = get_s(t, c);

        init_rotation_matrix(n, rotation_mat, c, s, index);
        transform_rotation(n, mat_cpy, mat_tag, i, j, s, c);

        a_tag_off_diag = square_off_diagonal(n, mat_tag);
        convergence = a_off_diag - a_tag_off_diag; /* As explained in 1.2.1-5 */
        a_off_diag = a_tag_off_diag;

        temp = mat_cpy;
        mat_cpy = mat_tag;
        mat_tag = temp;

        for (l = 0; l < n; l++)
        {
            i_temp = vectors[l][i] * c - vectors[l][j] * s;
            j_temp = vectors[l][i] * s + vectors[l][j] * c;
            vectors[l][i] = i_temp;
            vectors[l][j] = j_temp;
        }

        iter++;
    }

    /*Here we will enter the eigens and the eigenvectors we got from the vectors matrix*/
    for (i = 0; i < n; i++)
    {
        eigens[i].value = mat_cpy[i][i];
        for (j = 0; j < n; j++)
        {
            eigens[i].vector[j] = vectors[j][i];
        }
    }

    free_matrix(n, (void **)e_vectors);
vectors_cleanup:
    free_matrix(n, (void **)vectors);
rotation_mat_cleanup:
    free_matrix(n, (void **)rotation_mat);
mat_cpy_cleanup:
    free_matrix(n, (void **)mat_cpy);
mat_tag_cleanup:
    free_matrix(n, (void **)mat_tag);
end:
    return result;
}
