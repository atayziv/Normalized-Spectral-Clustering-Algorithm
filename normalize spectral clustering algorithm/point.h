#ifndef POINT_H
#define POINT_H

typedef struct point_t
{
    double *elements;
} point_t;

point_t *malloc_points(const size_t n, const size_t dim);
void free_points(const size_t n, point_t *points);
double calc_distance(const point_t p1, const point_t p2, const size_t dim);
int create_weight_matrix(const size_t n, double **weight_mat, point_t *points, const size_t dim);

#endif /* POINT_H */
