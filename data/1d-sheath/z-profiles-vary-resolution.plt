set terminal postscript eps enhanced colour font ",20"
set output "z-profiles-vary-resolution.eps"

set size square
set key bottom right
set key spacing 1.4
set xlabel 'z / {/Symbol l}_D'
set ylabel 'e{/Symbol D}{/Symbol f} / k_BT_e'
#set yrange [-0.1:0]
set xrange [0:20]

set arrow nohead from 0,0 to 20,0 dt 2

plot 'differences/width20-12x240diff.dat' w l lw 2 dt 1 lc rgb "#0072bd" t '12x240',\
     'differences/width20-25x500diff.dat' w l lw 2 dt 2 lc rgb "#a2142f" t '25x500',\
     'differences/width20-50x1000diff.dat' w l lw 2 dt 3 lc rgb "#77ac30" t '50x1000',\
     'differences/width20-100x2000diff.dat' w l lw 2 dt 4 lc rgb "#7e2f8e" t '100x2000'
