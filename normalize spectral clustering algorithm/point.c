#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "point.h"

point_t *malloc_points(const size_t n, const size_t dim)
{
    point_t *points = (point_t *)malloc(n * sizeof(point_t));
    size_t i;
    for (i = 0; i < n; i++)
    {
        points[i].elements = malloc(dim * sizeof(double));
    }
    return points;
}

void free_points(const size_t n, point_t *points)
{
    size_t i;
    for (i = 0; i < n; i++)
    {
        free(points[i].elements);
    }
    free(points);
}

/*A function to calculate the Euclidean distance between 2 points*/
double calc_distance(const point_t p1, const point_t p2, const size_t dim)
{
    double sum = 0;
    size_t i = 0;
    for (i = 0; i < dim; i++)
    {
        sum += pow(p1.elements[i] - p2.elements[i], 2.0);
    }
    return sqrt(sum);
}

/*A function to create a weight matrix as required in section 1.1.1*/
int create_weight_matrix(const size_t n, double **weight_mat, point_t *points, const size_t dim)
{
    size_t i, j;
    double distance;
    for (i = 0; i < n; i++)
    {
        for (j = i; j < n; j++)
        {
            if (i != j)
            {
                distance = calc_distance(points[i], points[j], dim);
                distance = exp(-distance / 2.0);
                weight_mat[j][i] = weight_mat[i][j] = distance;
            }
            else
            {
                weight_mat[i][i] = 0.0;
            }
        }
    }
    return 0;
}
