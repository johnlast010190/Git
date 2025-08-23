#/*---------------------------------------------------------------------------*\
#|       o        |
#|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
#|   o   O   o    |  Version : 4.4.0
#|    o     o     |  ENGYS Ltd. <http://engys.com/>
#|       o        |
#\*---------------------------------------------------------------------------
#License
#    This file is part of HELYXcore.
#    HELYXcore is based on OpenFOAM (R) <http://www.openfoam.org/>.
#
#    HELYXcore is free software: you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    HELYXcore is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#    for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with HELYXcore.  If not, see <http://www.gnu.org/licenses/>.

#Copyright
#    (c) 2016 OpenFOAM Foundation 
#
# Script
#     pendulumAndSpring.gnuplot
#
# Description
#     Creates an PostScript graph file of Test-pendulumAndSpring results
#
#------------------------------------------------------------------------------

reset

set xlabel "Time/[s]"
set ylabel "x"
set y2label "omega"

set ytics nomirror
set y2tics

set yrange [-1.5:1.5]
set y2range [-35:35]

set xzeroaxis

set terminal postscript eps color enhanced solid
set output "pendulumAndSpring.eps"

plot \
    "xVsTime" u 1:2 w l t "x", \
    "omegaVsTime" u 1:(57.29578*$2) w l axes x1y2 t "omega"

#------------------------------------------------------------------------------
