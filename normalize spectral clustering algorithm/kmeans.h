#ifndef KMEANS_H
#define KMEANS_H

#include "point.h"

typedef struct cluster
{
    point_t centroid;
    size_t size;
} cluster_t;

int fit(point_t *points, const size_t points_len, point_t *clusters, const size_t k, const size_t max_iter, const size_t dim, const float epsilon);

#endif /* KMEANS_H */
