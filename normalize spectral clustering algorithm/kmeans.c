#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "kmeans.h"
#include "point.h"

#define DELIM ','
#define EPSILON 0.01
#define DEFAULT_ITER 200
#define PRINT_ITERS
#define INVALID_INPUT_ERR "Invalid Input!\n"
#define GENERAL_ERR "Invalid Input!\n"
#define FUNC_SUCCESS 0
#define MALLOC_FAILED 1
#define FUNC_FAILED 2

typedef struct clustered_point
{
    point_t point;
    size_t cluster;
} clustered_point_t;

/* Assistive function for K-means - checks the closest cluster for a specific observation (vector) */
size_t find_closest_cluster(point_t point, cluster_t *centroids, size_t k, size_t dim)
{
    size_t i;
    double distance, max_distance = DBL_MAX;
    size_t selected = k + 1;
    for (i = 0; i < k; i++)
    {
        distance = calc_distance(point, centroids[i].centroid, dim);
        if (distance < max_distance)
        {
            max_distance = distance;
            selected = i;
        }
    }
    return selected;
}


int assign_to_clusters(clustered_point_t *points, size_t points_len, cluster_t *clusters, size_t k, size_t dim)
{
    size_t closest;
    size_t i;
    for (i = 0; i < points_len; i++)
    {
        closest = find_closest_cluster(points[i].point, clusters, k, dim);
        points[i].cluster = closest;
        clusters[closest].size++;
    }
    return FUNC_SUCCESS;
}

double update_centroids(clustered_point_t *points, size_t points_len, cluster_t *clusters, size_t k, size_t dim)
{
    size_t i, j;
    size_t cluster;
    double delta, centroid_delta;
    double *axis_sum;
    point_t *new_centroids;

    delta = .0;
    new_centroids = (point_t *)calloc(k, sizeof(point_t));
    if (new_centroids == NULL)
    {
        delta = -1;
        goto end;
    }
    axis_sum = (double *)calloc(dim * k, sizeof(double));
    if (NULL == axis_sum)
    {
        delta = -1;
        goto new_centroids_cleanup;
    }
    for (i = 0; i < points_len; i++)
    {
        cluster = points[i].cluster;
        for (j = 0; j < dim; j++)
        {
            axis_sum[cluster * dim + j] += points[i].point.elements[j] / clusters[cluster].size;
        }
    }
    for (i = 0; i < k; i++)
    {
        new_centroids[i].elements = &axis_sum[i * dim];
        centroid_delta = calc_distance(clusters[i].centroid, new_centroids[i], dim);
        if (centroid_delta > delta)
        {
            delta = centroid_delta;
        }

        for (j = 0; j < dim; j++)
        {
            clusters[i].centroid.elements[j] = axis_sum[i * dim + j];
        }
        clusters[i].size = 0;
    }

    free(axis_sum);
new_centroids_cleanup:
    free(new_centroids);
end:
    return delta;
}

int fit(point_t *points, const size_t points_len, point_t *clusters, const size_t k, const size_t max_iter, const size_t dim, const float epsilon)
{
    int result = FUNC_SUCCESS;
    size_t i;
    clustered_point_t *clustered_points;
    double max_delta;

    cluster_t *k_clusters;
    k_clusters = (cluster_t *)malloc(k * sizeof(cluster_t));
    if (NULL == k_clusters)
    {
        result = MALLOC_FAILED;
        goto end;
    }
    for (i = 0; i < k; i++)
    {
        k_clusters[i].centroid = clusters[i];
        k_clusters[i].size = 0;
    }

    clustered_points = (clustered_point_t *)malloc(sizeof(clustered_point_t) * points_len);
    if (NULL == clustered_points)
    {
        result = MALLOC_FAILED;
        goto k_clusters_cleanup;
    }
    for (i = 0; i < points_len; i++)
    {
        clustered_points[i].point = points[i];
    }
    for (i = 0; i < max_iter; i++)
    {
        assign_to_clusters(clustered_points, points_len, k_clusters, k, dim);
        max_delta = update_centroids(clustered_points, points_len, k_clusters, k, dim);
        if (max_delta < .0)
        {
            return FUNC_FAILED;
        }
        else if (max_delta <= epsilon)
        {
            break;
        }
    }

    free(clustered_points);
k_clusters_cleanup:
    free(k_clusters);
end:
    return result;
}
