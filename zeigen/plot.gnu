reset 
#unset key
set term x11 font "helvatica,24"
#set term postscript eps enhanced color font "helvatica,24"

set ytics nomirror
set y2tics nomirror

set log y
set format y '%.1e'
#set output ".eps"

set ylabel "time for each(s)"
set y2label "speedup"
set xlabel "matrix size N"
set title "Matrix Diagonalization CPU vs GPU"

set y2range [0:120]

set y2tics textcolor rgb "red"

plot 'comp.dat' u 1:2 w l lw 2 axis x1y1 title 'CPU-single core', 'comp.dat' u 1:3 w l lw 2 axis x1y1 title "GPU-one card", \
     'comp.dat' u 1:($2/$3) w lp lw 2 lc rgb "red" axis x1y2 title 'speedup'
