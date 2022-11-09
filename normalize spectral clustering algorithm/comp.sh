#!/bin/bash
# Script to compile and execute a c program

SRC_FILES="debug.c eigen.c input.c jacobi.c kmeans.c laplacian.c matrix.c point.c spkmeans.c"
# SRC_FILES="src/debug.c src/eigen.c src/input.c src/jacobi.c src/kmeans.c src/laplacian.c src/matrix.c src/point.c src/spkmeans.c"

#gcc -ansi -Wall -Wextra -Werror -pedantic-errors debug.c input.c jacobi.c matrix.c point.c spkmeans.c types.c -lm -o spkmeans
clang -ansi -Wall -Wextra -Werror -pedantic-errors $SRC_FILES -lm -o spkmeans

