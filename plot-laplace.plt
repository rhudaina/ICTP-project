set term png size 1280,1024
unset border
unset xtics
unset ytics
set palette rgb 33,13,10
set size ratio -1

set output "laplace.png"
plot 'solution.dat' with image
