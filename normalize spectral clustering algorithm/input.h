#ifndef INPUT_H
#define INPUT_H

#include <stdlib.h>
#include "point.h"

#define DELIM ','

size_t get_dimension(char *filename);
size_t get_lines_count(char *filename);
int read_points(char *filename, point_t *points, const size_t points_len, const size_t dim);
int read_matrix(char *filename, double **mat, const size_t n);
#endif /* INPUT_H */
