set terminal postscript eps enhanced colour font ",20"
set output "trench-phi.eps"

datafile = "depth2.dat"

set contour base
set cntrparam levels disc 0
unset surface
set table 'cont.dat'
splot datafile u 1:2:3, datafile u ($1)*(-1):2:3
unset table

#set dgrid3d 500,100 gauss 0.05

set contour base
set cntrparam level incremental -10,0.5,0
unset surface
set table 'cont2.dat'
splot datafile u ($3 < 0 ? $1 : 1/0):2:4, datafile u ($3 < 0 ? ($1)*(-1) : 1/0):2:4
unset table

#set cbrange [0:2]
#set cbtics 0,0.2,2
set cblabel 'e{/Symbol f} / k_BT_e'
#unset colorbox
set xlabel 'x / {/Symbol l}_D'
set ylabel 'z / {/Symbol l}_D' offset 1

set key off
#set palette maxcolors 10
set palette rgb 33,13,10
set size ratio 1.5
set xrange [-1:1]
set yrange [-2.1:1]
set style rect fc lt -1 fs solid 0.0 noborder
set obj rect from -0.1,-2.1 to 0.1,-2 front
set arrow nohead from -0.1,-2.1 to 0.1,-2.1 lt -1 lw 2 front
set arrow nohead from -0.1,-2 to 0.1,-2 lt -1 lw 2 front

plot datafile u 1:($2-2.1):(($3 < 0)  ? $4 : 1/0) w image,\
     datafile u ($1)*(-1):($2-2.1):(($3 < 0)  ? $4 : 1/0) w image,\
     'cont2.dat' u 1:($2-2.1) w l lt -1 lw 2,\
     'cont.dat' u 1:($2-2.1) w filledcurves x1 lc rgb "white",\
     'cont.dat' u 1:($2-2.1) w filledcurves x1 lc rgb "white",\
     'cont.dat' u 1:($2-2.1) w filledcurves x1 lc rgb "white",\
     'cont.dat' u 1:($2-2.1) w filledcurves x1 lc rgb "white",\
     'cont.dat' u 1:($2-2.1) w l lt -1 lw 2#,\
#     'cont2.dat' w l lt -1 lw 2#,\
#     'cont2.dat' u (($1<0.11 && $1>0.09) ? 0 : 1/0):2:3 w labels


