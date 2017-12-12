set terminal postscript eps enhanced colour font ",20"
set output "enhancement-factors.eps"

set size square
set key bottom right
set key spacing 1.4
set xlabel 'Aspect ratio'
set yrange [1:4.5]
set xrange [0:5]

zeta(x) = sqrt(x*x-1)
nu(x) = zeta(x)**3 / (x*log(x + zeta(x)) - zeta(x))


plot 'enhancement-factors.dat' u 1:($2)*(-1) w l lw 3 dt 1 lc rgb "#0072bd" t 'Apex field strength' smooth bezier,\
     'enhancement-factors.dat' u 1:3 w l lw 3 dt 2 lc rgb "#a2142f" t 'Field-enhancement factor' smooth bezier,\
     nu(x) w l lw 3 dt 4 lc rgb "#77ac30" t 'Vacuum field-enhancement factor'
