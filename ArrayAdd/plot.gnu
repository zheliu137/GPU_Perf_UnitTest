reset 
#unset key
set term x11 font "helvatica,24"
#set term postscript eps enhanced color font "helvatica,24"

set ytics nomirror
set y2tics nomirror

set log y
set format y '%.1e'
set log x
# set format x '10^{%T}'
# set format x "10^{%L}"
set format x "1E{%L}"
#set output ".eps"

set ylabel "time for each element(s)"
set y2label "speedup"
set xlabel "vector length N"
set title "vector summation CPU vs GPU"

set y2range [0:140]

set y2tics textcolor rgb "red"

# plot 'comp.dat' u 1:2 w l lw 2 axis x1y1 title "GPU-one card"

plot 'comp.dat' u 1:2 w l lw 2 axis x1y1 title 'CPU-single core', 'comp.dat' u 1:3 w l lw 2 axis x1y1 title "GPU-one card", \
     'comp.dat' u 1:($2/$3) w lp lw 2 lc rgb "red" axis x1y2 title 'speedup'
