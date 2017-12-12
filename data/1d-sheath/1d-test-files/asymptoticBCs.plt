set terminal postscript eps enhanced colour font ",20"
set output "asymptoticBCs.eps"

set multiplot

set size square
set key bottom right
set key spacing 1.4
set xlabel 'z / {/Symbol l}_D'
set ylabel 'e~{/Symbol f} / k_BT_e'
set xrange [15:100]

plot 'unperturbed.dat' u 1:2 w l lw 2 dt 2 lc rgb "#0072bd" t 'Exact solution',\
     -6/(x*x) w l lw 2 dt 4 lc rgb "#a2142f" t 'Approximate solution',\
     '' u 1:($2+6/(($1)*($1))) w l lw 2 dt 1 lc rgb "#77ac30" t 'Discrepancy'

set origin 0.4,0.4
set size 0.4,0.4
set xrange [50:100]
unset xlabel
unset ylabel

plot 'unperturbed.dat' u 1:2 w l lw 2 dt 2 lc rgb "#0072bd" t '',\
     -6/(x*x) w l lw 2 dt 4 lc rgb "#a2142f" t '',\
     '' u 1:($2+6/(($1)*($1))) w l lw 2 dt 1 lc rgb "#77ac30" t ''
