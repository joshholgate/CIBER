set terminal postscript eps enhanced colour font ",20"
set output "x-profile.eps"

set size square
set key bottom right
set key spacing 1.4
set xlabel 'x / {/Symbol l}_D'
set ylabel 'e~{/Symbol f} / k_BT_e'
set xrange [0:1]

# Some gnuplot trickery is needed to interpolate between grid points
a=0
b=0

plot '12x240-microscopic.dat' u ($2==0.291667 ? $1 : 1/0):4 w p ps 1.5 lc rgb "#d95319" pt 5 t '12x240',\
     '25x500-microscopic.dat' u ($2==0.3 ? $1 : 1/0):4 w p ps 1.5 lc rgb "#edb120" pt 9 t '25x500',\
     '50x1000-microscopic.dat' u ($2==0.29 ? (a=$4, 1/0) : ($2==0.31 ? (b=$4, $1) : 1/0) ):(a+b)/2 w p ps 1 lc rgb "#7e2f8e" pt 7 t '50x1000',\
     '100x2000-microscopic.dat' u ($2==0.295 ? (a=$4, 1/0) : ($2==0.305 ? (b=$4, $1) : 1/0) ):(a+b)/2 w p ps 1 lc rgb "#77ac30" pt 13 t '100x2000',\
     -2.61777-0.0301386*cos(3.14*x) w l lw 3 lc rgb "black" t 'linear theory'

# Note: see ELIPS for code for calculation of linear theory

