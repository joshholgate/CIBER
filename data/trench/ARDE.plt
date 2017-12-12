set terminal postscript eps enhanced colour font ",20"
set output "ARDE.eps"

set size square
set key top left
set key spacing 1.4
set key title 'Trench depth / {/Symbol l}_D'
set xlabel 'z / {/Symbol l}_D'
set ylabel 'n / n_0' offset 1
set xrange [-32:20]
set arrow nohead from 0,0 to 0,1 lw 2 dt 2 back

set arrow nohead from -16,0 to -16,0.055 lw 2 dt 2 lc rgb "#a2142f"
set arrow nohead from -8,0 to -8,0.1 lw 2 dt 2 lc rgb "#77ac30"
set arrow nohead from -4,0 to -4,0.15 lw 2 dt 2 lc rgb "#7e2f8e"
set arrow nohead from -2,0 to -2,0.175 lw 2 dt 2 lc rgb "#edb120"
set arrow nohead from -1,0 to -1,0.2 lw 2 dt 2 lc rgb "#d95319"

plot 'depth32.dat' u ($1<0.011 && (($2)>0.101) ? ($2-32.1) : 1/0):($9) w l lw 2 lc rgb "#0072bd" t '32',\
     'depth16.dat' u ($1<0.011 && (($2)>0.101) ? ($2-16.1) : 1/0):($9) w l lw 2 lc rgb "#a2142f" t '16',\
     'depth8.dat' u ($1<0.011 && (($2)>0.101) ? ($2-8.1) : 1/0):($9) w l lw 2 lc rgb "#77ac30" t '8',\
     'depth4.dat' u ($1<0.011 && (($2)>0.101) ? ($2-4.1) : 1/0):($9) w l lw 2 lc rgb "#7e2f8e" t '4',\
     'depth2.dat' u ($1<0.011 && (($2)>0.101) ? ($2-2.1) : 1/0):($9) w l lw 2 lc rgb "#edb120" t '2',\
     'depth1.dat' u ($1<0.011 && (($2)>0.101) ? ($2-1.1) : 1/0):($9) w l lw 2 lc rgb "#d95319" t '1',\
     'depth0.dat' u ($1<0.011 && (($2-0.1)>0.101) ? ($2-0.1) : 1/0):($9) w l lw 2 lc rgb "#4dbeee" t '0'
