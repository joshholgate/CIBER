set terminal postscript eps enhanced colour font ",20"
set output "z-profiles-vary-width.eps"

set size square
set key top right
set key spacing 1.4
set xlabel 'z / {/Symbol l}_D'
set ylabel 'e{/Symbol D}{/Symbol f} / k_BT_e'
#set yrange [-0.1:0]
set xrange [0:40]

set arrow nohead from 0,0 to 40,0 dt 2

plot 'differences/width10-50x500diff.dat' w l lw 2 dt 1 lc rgb "#0072bd" t 'z_{max} = 10{/Symbol l}_D',\
     'differences/width20-50x1000diff.dat' w l lw 2 dt 2 lc rgb "#a2142f" t 'z_{max} = 20{/Symbol l}_D',\
     'differences/width30-50x1500diff.dat' w l lw 2 dt 3 lc rgb "#77ac30" t 'z_{max} = 30{/Symbol l}_D',\
     'differences/width40-50x2000diff.dat' w l lw 2 dt 4 lc rgb "#7e2f8e" t 'z_{max} = 40{/Symbol l}_D'
