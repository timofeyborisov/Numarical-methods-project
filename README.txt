to compile program - enter this in cmd: 
gcc code.c -o interpolation -lm

to run - enter this:
./interpolation

default number of dots - 101

commands for gnuplot for data visualization:

plot "f1.dat" title "first dataset" lt 7 with lp, "f1_poly_3.dat" title "Lagrange polynomial (n = 3)" lc 3 lt 7 with l
plot "f1.dat" title "first dataset" lt 7 with lp, "f1_poly_5.dat" title "Lagrange polynomial (n = 5)" lc 3 lt 7 with l
plot "f1.dat" title "first dataset" lt 7 with lp, "f1_poly_9.dat" title "Lagrange polynomial (n = 9)" lc 3 lt 7 with l
plot "f1.dat" title "first dataset" lt 7 with lp, "f1_poly_17.dat" title "Lagrange polynomial (n = 17)" lc 3 lt 7 with l

plot "f2.dat" title "second dataset" lt 7 with lp, "f2_spln_3.dat" title "cubic spline (n = 3)" lc 3 lt 7 with l
plot "f2.dat" title "second dataset" lt 7 with lp, "f2_spln_5.dat" title "cubic spline (n = 5)" lc 3 lt 7 with l
plot "f2.dat" title "second dataset" lt 7 with lp, "f2_spln_9.dat" title "cubic spline (n = 9)" lc 3 lt 7 with l
plot "f2.dat" title "second dataset" lt 7 with lp, "f2_spln_17.dat" title "cubic spline (n = 17)" lc 3 lt 7 with l