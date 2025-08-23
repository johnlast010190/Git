#!/bin/sh
 
gnuplot -p << EOF                                                               
	set grid                                                                        
	set title 'Displacement in the middle of the tube                                        
	set xlabel 'Time [s]'                                                           
	set ylabel 'Displacement [m]'
	set term pngcairo enhanced size 900,654
	set output "results/tutorials-elastic-tube-3d-displacement-all-watchpoints.png"
	plot "solid/precice-Solid-watchpoint-Tube-Midpoint.log" using 1:5 with lines title "HELYX-CalculiX circumferential", \
	     "solid/precice-Solid-watchpoint-Tube-Midpoint.log" using 1:6 with lines title "HELYX-CalculiX radial", \
	     "solid/precice-Solid-watchpoint-Tube-Midpoint.log" using 1:7 with lines title "HELYX-CalculiX axial"
EOF
