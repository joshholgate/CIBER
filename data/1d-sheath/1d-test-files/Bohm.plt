set terminal postscript eps enhanced colour font ",20"
set output "Bohm.eps"

set size square
set key bottom right
set key spacing 1.4
set xlabel 'z / {/Symbol l}_D'
set ylabel 'e~{/Symbol f} / k_BT_e'
#set yrange [-0.1:0]
set xrange [0:20]
set key off

plot 'unperturbed.dat' u 1:2 w l lw 3 lc rgb "black"

#plot 'RK4-loopXi-Ksq0.01.dat' every 9 u 1:6 w lp lw 1.5 lc rgb "#d95319" lt 6 t 'k^2{/Symbol l}_D^2=0.01',\
#     'RK4-loopXi-Ksq0.1.dat' every 9 u 1:6 w lp lw 1.5 lc rgb "#edb120" lt 7 t 'k^2{/Symbol l}_D^2=0.1',\
#     'RK4-loopXi-Ksq1.dat' every 9 u 1:6 w lp lw 1.5 lc rgb "#7e2f8e" pt 4 t 'k^2{/Symbol l}_D^2=1',\
#     'RK4-loopXi-Ksq10.dat' every 9 u 1:6 w lp lw 1.5 lc rgb "#77ac30" pt 5 t 'k^2{/Symbol l}_D^2=10'

#plot 'unperturbed.dat' u 1:2 w l lw 3 lc rgb "black" t '1D solution',\
#     '100x1600.dat' every 2 u ($1<0.011 ? ($2-1) : 1/0):4 w l dt 2 lw 3 lc rgb "#77ac30" t '100x1600'
#     '100x1600.dat' every 4 u ($1<0.06 ? ($2-1) : 1/0):4 w p ps 1 lc rgb "#77ac30" pt 13 t '100x1600',\
#     'RK4-soln-K9.8696.dat' u 1:($2+$7) w l lw 3 lc rgb "black" t 'linear theory'
#     'RK4-soln-K9.8696.dat' u 1:($2-$7) w l lw 2 lc rgb "#77ac30" t 'linear PT, x={/Symbol l}_D',\
#     '50x800.dat' u ($1>0.981 ? ($2-1) : 1/0):4 lc rgb "#7e2f8e" pt 7 t 'CIBER, x={/Symbol l}_D',\
#     '25x400.dat' u ($1>0.971 ? ($2-1) : 1/0):4 t 'CIBER, x=0',\
#     '12x196.dat' u ($1>0.95 ? ($2-1) : 1/0):4 t 'CIBER, x=0'

