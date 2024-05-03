#!/bin/bash
clear
rm *.png
rm *.dat
rm *.out

# mpicc -O3 jacobi.c
echo -e "Running..."
mpicc -O3 main.c
mpirun -np 2 ./a.out 500 1000000

echo -e "Plotting..."
gnuplot plot-error.plt
open error.png

gnuplot plot-laplace.plt
open laplace.png

echo -e "Done!\n"
