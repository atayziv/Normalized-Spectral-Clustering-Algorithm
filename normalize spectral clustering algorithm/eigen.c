#include <stdlib.h>
#include <float.h>
#include <math.h>

#include "eigen.h"

/*A function to allocate space for a matrix of eigenvectors*/
eigen_t *malloc_eigens(const size_t n)
{
    size_t i;
    eigen_t *eigens = malloc(n * sizeof(eigen_t));
    if (NULL == eigens)
    {
        return NULL;
    }
    for (i = 0; i < n; i++)
    {
        eigens[i].vector = malloc(n * sizeof(double));
    }
    return eigens;
}

/*Function for freeing memory of eigenvectors*/
void free_eigens(const size_t n, eigen_t *eigens)
{
    size_t i;
    for (i = 0; i < n; i++)
    {
        free(eigens[i].vector);
    }
    free(eigens);
}

/*A function for building a matrix of eigenvectors*/
int build_matrix_from_eigens(const size_t n, const size_t k, eigen_t *eigens, double **mat)
{
    size_t i, j;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < k; j++)
        {
            mat[i][j] = eigens[j].vector[i];
        }
    }
    return 0;
}

/*A function for comparing eigenvalues in order to sort them in descending order in section 1.3*/
int compare_eigenvalues(const void *a, const void *b)
{
    if (((eigen_t *)a)->value > ((eigen_t *)b)->value)
    {
        return -1;
    }
    else if (((eigen_t *)a)->value < ((eigen_t *)b)->value)
    {
        return 1;
    }
    return 0;
}

/*In order to determine the value of K, in this function we find the largest delta as explained in section 1.3*/
size_t find_eigengap_max(const size_t n, eigen_t *eigens)
{
    size_t res, i;
    double diff, max;
    max = -DBL_MAX;
    res = 0;
    for (i = 0; i < n / 2; i++)
    {
        diff = fabs(eigens[i].value - eigens[i + 1].value);
        if (diff > max)
        {
            max = diff;
            res = i;
        }
    }
    return res + 1;
}
