reset 
#unset key
set term x11 font "helvatica,24"
#set term postscript eps enhanced color font "helvatica,24"

set ytics nomirror
set y2tics nomirror

set log y
set log y2
set format y '%.1e'
set log x

set key left top
# set format x '10^{%T}'
# set format x "10^{%L}"
set format x "1E{%L}"
#set output ".eps"

set ylabel "time for each summation(s)"
set y2label "speedup"
set xlabel "vector length N"
set title "vector summation CPU vs GPU"

#set y2range [0:100]
set yrange [1e-13:1e-7]
set y2range [0:5000]

set y2tics textcolor rgb "red"

# plot 'comp.dat' u 1:2 w l lw 2 axis x1y1 title "GPU-one card"

# plot 'comp.dat' u 1:2 w l lw 2 axis x1y1 title 'CPU-single core', 'comp.dat' u 1:3 w l lw 2 axis x1y1 title "GPU-one card", \
#      'comp.dat' u 1:($2/$3) w lp lw 2 lc rgb "red" axis x1y2 title 'speedup'
plot 'comp.dat' u 1:2 w l lw 2 axis x1y1 title 'CPU-single core', 'comp.dat' u 1:3 w l lw 2 axis x1y1 title "GPU-one card-v1", \
     'comp.dat' u 1:4 w l lw 2 axis x1y1 title "GPU-one card-v2", \
     'comp.dat' u 1:($2/$3) w lp lw 2 lc rgb "red" axis x1y2 title 'speedup-v1', \
     'comp.dat' u 1:($2/$4) w lp lw 2 linetype 1 dashtype 5 lc rgb "red" axis x1y2 title 'speedup-v2', \
     'comp.dat' u 1:5 w lp lw 2 title 'gpu-zgemv'
