#define PY_SSIZE_T_CLEAN

#include <Python.h>

#include "spkmeans.h"
#include "debug.h"

static double **malloc_matrix(const size_t n, const size_t dim)
{
    size_t i;
    double **points = malloc(n * sizeof(void *));
    if (NULL == points)
    {
        return NULL;
    }
    for (i = 0; i < n; i++)
    {
        points[i] = calloc(dim, sizeof(double));
    }

    return points;
}
static void free_matrix(const size_t n, double **mat)
{
    size_t i;
    for (i = 0; i < n; i++)
    {
        free(mat[i]);
    }
    free(mat);
}

PyObject *create_result(point_t *clusters, size_t k, size_t dim)
{
    PyObject *result = PyList_New(k);
    PyObject *centroid;
    size_t i, j;
    for (i = 0; i < k; ++i)
    {
        centroid = PyTuple_New(dim);
        for (j = 0; j < dim; ++j)
        {
            PyTuple_SetItem(centroid, j, PyFloat_FromDouble(clusters[i].elements[j]));
        }
        PyList_SetItem(result, i, centroid);
    }
    return result;
}

int parse_points(PyObject *points_obj, point_t *points, const size_t points_len, const size_t dim)
{

    PyObject *cur_point;
    size_t i, j;

    for (i = 0; i < points_len; i++)
    {
        cur_point = PyList_GetItem(points_obj, i);
        for (j = 0; j < dim; j++)
        {
            points[i].elements[j] = PyFloat_AsDouble(PyTuple_GetItem(cur_point, j));
        }
    }
    return OK;
}
int parse_matrix(PyObject *in_mat, double **mat, const size_t points_len, const size_t dim)
{
    PyObject *cur_row;
    size_t i, j;

    for (i = 0; i < points_len; i++)
    {
        cur_row = PyList_GetItem(in_mat, i);
        for (j = 0; j < dim; j++)
        {
            mat[i][j] = PyFloat_AsDouble(PyList_GetItem(cur_row, j));
        }
    }
    return OK;
}

size_t get_dim(PyObject *points_obj)
{
    return PyTuple_Size(PyList_GetItem(points_obj, 0));
}

PyObject *create_py_matrix(const size_t n, const size_t k, double **mat)
{
    PyObject *cur_row, *result;
    size_t i, j;
    result = PyList_New(n);
    for (i = 0; i < n; ++i)
    {
        cur_row = PyList_New(k);
        for (j = 0; j < k; ++j)
        {
            PyList_SetItem(cur_row, j, PyFloat_FromDouble(mat[i][j]));
        }
        PyList_SetItem(result, i, cur_row);
    }
    return result;
}
static PyObject *calc(PyObject *self, PyObject *args, goal_e goal)
{
    error_e result = OK;
    PyObject *data_points = NULL, *result_obj = NULL;
    size_t dim, points_len, k;
    double **mat;
    point_t *points;
    if (!PyArg_ParseTuple(args, "Ol", &data_points, &k))
    {
        return NULL;
    }
    dim = get_dim(data_points);
    points_len = PyObject_Length(data_points);
    points = malloc_points(points_len, dim);
    if (NULL == points)
    {
        return NULL;
    }
    parse_points(data_points, points, points_len, dim);
    mat = malloc_matrix(points_len, points_len);
    if (NULL == mat)
    {
        result = MALLOC_ERROR;
        goto cleanup_points;
    }
    result = calc_matrix(points_len, points, dim, goal, mat, &k);

    if (NORMALIZED_EIGEN_MATRIX == goal)
    {
        result_obj = create_py_matrix(points_len, k, mat);
    }
    else
    {
        result_obj = create_py_matrix(points_len, points_len, mat);
    }

    free_matrix(points_len, mat);
cleanup_points:
    free_points(points_len, points);
    if (result != OK)
    {
        Py_RETURN_NONE;
    }
    return result_obj;
}

static PyObject *calc_wam(PyObject *self, PyObject *args)
{
    return calc(self, args, WEIGHT_MATRIX);
}

static PyObject *calc_ddg(PyObject *self, PyObject *args)
{
    return calc(self, args, DIAGONAL_DEGREE_MATRIX);
}

static PyObject *calc_lnorm(PyObject *self, PyObject *args)
{
    return calc(self, args, NORMALIZED_GRAPH_LAPLACIAN);
}

static PyObject *calc_spk(PyObject *self, PyObject *args)
{
    return calc(self, args, NORMALIZED_EIGEN_MATRIX);
}

static PyObject *kmeans_fit(PyObject *self, PyObject *args)
{

    PyObject *centroids, *data_points = NULL, *result = NULL;
    point_t *points, *clusters;
    size_t dim, k, points_len, max_iter;
    float epsilon;
    int fit_result = 1;
    /* Parse arguments */
    if (!PyArg_ParseTuple(args, "OOnf", &centroids, &data_points, &max_iter, &epsilon))
    {
        return NULL;
    }
    k = PyObject_Length(centroids);
    dim = get_dim(data_points);
    points_len = PyObject_Length(data_points);
    points = malloc_points(points_len, dim);
    parse_points(data_points, points, points_len, dim);
    if (NULL == points)
    {
        goto end;
    }
    clusters = malloc_points(k, dim);
    if (NULL == clusters)
    {
        fit_result = MALLOC_ERROR;
        goto points_free;
    }
    parse_points(centroids, clusters, k, dim);
    fit_result = kmeans(points, points_len, clusters, k, max_iter, dim, epsilon);
    if (0 == fit_result)
    {
        result = create_result(clusters, k, dim);
    }
    free_points(k, clusters);
points_free:
    free_points(points_len, points);
end:
    if (OK != fit_result)
    {
        Py_RETURN_NONE;
    }
    else
    {
        return result;
    }
}
static PyObject *calc_jacobi(PyObject *self, PyObject *args)
{
    PyObject *in_mat = NULL, *out_values = NULL, *out_vectors = NULL;
    size_t dim;
    double **mat, *values, **vectors;
    size_t i;
    if (!PyArg_ParseTuple(args, "O", &in_mat))
    {
        return NULL;
    }
    dim = PyObject_Length(in_mat);

    mat = malloc_matrix(dim, dim);
    parse_matrix(in_mat, mat, dim, dim);

    values = (double *)malloc(dim * sizeof(double));
    vectors = malloc_matrix(dim, dim);
    calc_eigen_values_vectors(dim, mat, values, vectors);
    out_values = PyList_New(dim);
    for (i = 0; i < dim; i++)
    {
        PyList_SetItem(out_values, i, PyFloat_FromDouble(values[i]));
    }
    out_vectors = create_py_matrix(dim, dim, vectors);

    return PyTuple_Pack(2, out_values, out_vectors);
}

static PyMethodDef spkmeansMethods[] =
    {

        {"wam",
         calc_wam,
         METH_VARARGS,
         PyDoc_STR("Calculate Weighted Adjacency Matrix.")},
        {"ddg",
         calc_ddg,
         METH_VARARGS,
         PyDoc_STR("Calculate Diagonal Degree Graph.")},
        {"lnorm",
         calc_lnorm,
         METH_VARARGS,
         PyDoc_STR("Calculate Normalized Graph Laplacian.")},
        {"jacobi",
         calc_jacobi,
         METH_VARARGS,
         PyDoc_STR("Calculate eigen values and vectors of a matrix using Jacobi algorithm.")},
        {"spk",
         calc_spk,
         METH_VARARGS,
         PyDoc_STR("Calculate the normalized eigen matrix.")},
        {"kmeans_fit",
         kmeans_fit,
         METH_VARARGS,
         PyDoc_STR("runs kmeans algorithm")},
        {NULL, NULL, 0, NULL},
};

static struct PyModuleDef moduledef =
    {
        PyModuleDef_HEAD_INIT,
        "spkm",
        NULL,
        -1,
        spkmeansMethods};

PyMODINIT_FUNC
PyInit_spkm(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m)
    {
        return NULL;
    }
    return m;
}
