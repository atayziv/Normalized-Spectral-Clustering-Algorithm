#include <stdio.h>

#include "debug.h"
#include "eigen.h"

/*Function to print matrices*/
void print_matrix(const size_t n, const size_t k, double **matrix)
{
    size_t i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < k; j++)
        {
            printf("%.4f", matrix[i][j]);
            if (j < k - 1)
            {
                printf(",");
            }
        }
        printf("\n");
    }
}

/*Function to print arrays*/
void print_array(const size_t n, double *array)
{
    size_t i;
    for (i = 0; i < n; i++)
    {
        printf("%.4f", array[i]);
        if (i < n - 1)
        {
            printf(",");
        }
    }
    printf("\n");
}

/*Function for printing eigenvalues and eigenvectors*/
void print_eigen(const size_t n, eigen_t *eigen)
{
    size_t i, j;
    for (i = 0; i < n; i++)
    {
        printf("%.4f", eigen[i].value);
        if (i < n - 1)
        {
            printf(",");
        }
    }
    printf("\n");
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            printf("%.4f", eigen[j].vector[i]);
            if (j < n - 1)
            {
                printf(",");
            }
        }
        printf("\n");
    }
}
