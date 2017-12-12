set terminal postscript eps enhanced colour font ",20"
set output "energy.eps"

set size square
set key top right
set key spacing 1.4
set xlabel 't / {/Symbol l}_Dc_s^{-1}'
set ylabel 'Maximum ion energy error'
set yrange [0:0.02]
set xrange [0:500]

set arrow nohead from 0,0 to 20,0 dt 2

plot '12x240-macroscopic.dat' every 50 u 1:3 w l lw 2 dt 1 lc rgb "#0072bd" t '12x240',\
     '25x500-macroscopic.dat' every 50 u 1:3 w l lw 2 dt 2 lc rgb "#a2142f" t '25x500',\
     '50x1000-macroscopic.dat' every 50 u 1:3 w l lw 2 dt 3 lc rgb "#77ac30" t '50x1000',\
     '100x2000-macroscopic.dat' every 50 u 1:3 w l lw 2 dt 4 lc rgb "#7e2f8e" t '100x2000'

