set terminal postscript eps enhanced colour font ",20"
set output "field-profile.eps"

datafile = "nu1.dat"

set contour base
set cntrparam levels disc 0
unset surface
set table 'cont.dat'
splot datafile u 1:2:3, datafile u ($1)*(-1):2:3
unset table

## NOTE: this line makes the script take AGES to run, but gives slightly smoother contours for publication-stabdard plot
set dgrid3d 500,100 gauss 0.05

set contour base
set cntrparam level incremental 0, 0.2, 2
unset surface
set table 'cont2.dat'
splot datafile u ($3 < 0 ? $1 : 1/0):2:(sqrt(($5)*($5)+($6)*($6))), datafile u ($3 < 0 ? ($1)*(-1) : 1/0):2:(sqrt(($5)*($5)+($6)*($6)))
unset table

set cbrange [0:2]
set cbtics 0,0.2,2
set cblabel '|E| / (e{/Symbol l}_D/k_BT_e)'
#unset colorbox
set xlabel 'r / {/Symbol l}_D'
set ylabel 'z / {/Symbol l}_D' offset 1

set key off
#set palette maxcolors 10
set palette rgb 33,13,10
set size ratio 1
set xrange [-1:1]
set yrange [0:2]

set style rect fc lt -1 fs solid 0.0 noborder
set obj rect from -0.01,0 to 0.01,0.3 front
set arrow nohead from -0.01,0 to 0.01,0 lt -1 lw 2 front
set arrow nohead from -0.01,0.3 to 0.01,0.3 lt -1 lw 2 front

plot datafile u 1:2:(($3 < 0)  ? (sqrt(($5)*($5)+($6)*($6))) : 1/0) w image,\
     datafile u ($1)*(-1):2:(($3 < 0)  ? (sqrt(($5)*($5)+($6)*($6))) : 1/0) w image,\
     'cont2.dat' w l lt -1 lw 2,\
     'cont.dat' w filledcurves x1 lc rgb "white",\
     'cont.dat' w filledcurves x1 lc rgb "white",\
     'cont.dat' w filledcurves x1 lc rgb "white",\
     'cont.dat' w filledcurves x1 lc rgb "white",\
     'cont.dat' w l lt -1 lw 2#,\
#     'cont2.dat' w l lt -1 lw 2#,\
#     'cont2.dat' u (($1<0.11 && $1>0.09) ? 0 : 1/0):2:3 w labels


