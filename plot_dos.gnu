#!/usr/bin/gnuplot
set title "Density of states for test system"
set xlabel "E (eV)"

set term pngcairo
set output "DOS_total.png"

plot "TDOS.dat" using 1:2 title "Up" with lines, \
"TDOS.dat" using 1:3 title "Down" with lines, \
"TDOS.dat" using 1:($2 + abs($3)) with lines title "Total = up + abs(down)"
