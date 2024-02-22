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
set title "Matrix Multiplication CPU vs GPU"

set y2range [0:140]

set y2tics textcolor rgb "red"

plot 'comp.dat' u 1:2 w l lw 2 lc rgb "blue" axis x1y1 title 'CPU-single core', 'comp.dat' u 1:3 w l lw 2 lc rgb "dark-violet" axis x1y1 title "GPU-one card", \
     'comp.dat' u 1:($2/$3) with lines lw 2 lc rgb "red" axis x1y2 title 'speedup'
