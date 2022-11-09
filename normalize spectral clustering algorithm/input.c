#include <stdlib.h>
#include <stdio.h>

#include "input.h"

/*The value of the dimension will be determined by the number of commas in the line (+1) in the given file*/
size_t get_dimension(char *filename)
{
    char current_c = '\0';
    size_t dim = 0;
    FILE *points_file = fopen(filename, "r");
    if (points_file == NULL)
    {
        return 0;
    }

    while ((current_c = fgetc(points_file)) != '\n')
    {
        if (current_c == DELIM)
        {
            dim++;
        }
    }

    fclose(points_file);
    return dim + 1;
}

/*n (the number of points) will be determined by the number of lines in the given text*/
size_t get_lines_count(char *filename)
{
    size_t lines = 0;
    char current_c;
    FILE *points_file = fopen(filename, "r");
    if (NULL == points_file)
    {
        return 0;
    }

    while (!feof(points_file))
    {
        current_c = fgetc(points_file);
        if (current_c == '\n')
        {
            lines++;
        }
    }
    fclose(points_file);
    return lines;
}

/*In this function we will read the text, and insert the points according to the format into the matrix*/
int read_points(char *filename, point_t *points, const size_t points_len, const size_t dim)
{
    FILE *points_file;
    size_t i, j;
    points_file = fopen(filename, "r");
    if (NULL == points_file)
    {
        return -1;
    }

    for (i = 0; i < points_len; i++)
    {
        for (j = 0; j < dim; j++)
        {
            fscanf(points_file, "%lf,", &(points[i].elements[j]));
        }
        fseek(points_file, 1, SEEK_CUR);
    }
    fclose(points_file);
    return 0;
}

/*We will enter the 'mat' matrix point by point from the text*/
int read_matrix(char *filename, double **mat, size_t n)
{
    FILE *points_file;
    size_t i, j;
    points_file = fopen(filename, "r");
    if (NULL == points_file)
    {
        return -1;
    }

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            fscanf(points_file, "%lf,", &(mat[i][j]));
        }
        fseek(points_file, 1, SEEK_CUR);
    }
    fclose(points_file);
    return 0;
}
