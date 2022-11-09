
#include <stdlib.h>
#include <math.h>

#include "laplacian.h"
#include "matrix.h"

int create_laplacian_matrix(const size_t n, double **w_mat, double **d_mat, double **l_mat)
{
    size_t i, j;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            l_mat[i][j] = i != j ? -w_mat[i][j] : d_mat[i][j];
        }
    }
    return 0;
}

/*A function that generates the Laplacian matrix as explained in section 1.1.3*/
int create_normalized_laplacian_matrix(const size_t n, double **l_mat, double **n_mat)
{
    int result = 0;
    size_t i;
    double **inverse_l, **mul_mat;

    inverse_l = (double **)malloc_matrix(n, n, sizeof(double));
    if (NULL == inverse_l)
    {
        result = 1;
        goto end;
    }
    mul_mat = (double **)malloc_matrix(n, n, sizeof(double));
    if (NULL == mul_mat)
    {
        result = 1;
        goto l_mat_cleanup;
    }

    for (i = 0; i < n; i++)
    {
        inverse_l[i][i] = sqrt(1 / l_mat[i][i]);
    }
    multiply_mat(n, inverse_l, l_mat, n_mat);
    multiply_mat(n, n_mat, inverse_l, n_mat);

    free_matrix(n, (void **)mul_mat);
l_mat_cleanup:
    free_matrix(n, (void **)inverse_l);
end:
    return result;
}
