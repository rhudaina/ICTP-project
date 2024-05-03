set term pngcairo font "arial,20" size 1280,1024
set key noautotitle

set xrange [0:2000]

set output "error.png"
set title "L2-error"
set xlabel "iterations"
plot 'error.dat' w l lw 3 lc rgb "blue"
