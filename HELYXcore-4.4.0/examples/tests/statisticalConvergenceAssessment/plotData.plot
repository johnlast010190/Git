#set terminal qt persist
set title "Convergence Assesmnet"
set xlabel "Time"
set ylabel "Value"
set key outside
set grid

# Define the input file
filename = "postProcessing/RASCC/0/confidenceInterval.dat"

# Define a function to simulate progressive plotting
stats filename using 2 nooutput
maxLines = STATS_records

set style fill transparent solid 0.3 noborder
#set yr [0.25:0.35]
#set yr [0.029:0.033]
#set yr [0.22:0.27]

# Loop to update plot dynamically
#do for [i=1:maxLines-1] {
#    plot \
#        filename every ::0::i using 2:4 w l t "Value" lc rgb "blue", \
#        filename every ::0::i using 2:7:6 w filledcu t "Confidence Interval", \
#        filename every ::0::i using 2:5 with lines title "Mean" lc rgb "green", \
#        filename every ::i::i u 3:5 w p pt 7 ps 2 t "Transient Time"
#    pause 0.05
#}

#Loop to update plot dynamically
#do for [i=1:maxLines-1] {
#    unset arrow
#
#    # Extract the transient time (column 3) from the current row (row index i-1)
#    stats filename u 3 every ::i::i nooutput
#    T_t = real(word(system(sprintf("sed -n '%dp' %s", i+1, filename)), 3))
#
#    # Set a vertical line at T_t spanning the full y-range
#    set arrow from T_t, graph 0 to T_t, graph 1 nohead lc rgb "black" lw 2
#
#    plot \
#        filename every ::0::i using 2:5 with lines title "Mean" lc rgb "green", \
#        filename every ::0::i using 2:8 w l t "Filter" lc rgb "red", \
#        filename every ::0::i using 2:7:6 w filledcu t "Confidence Interval", \
##        filename every ::0::i using ( ($2 < T_t) ? $2 : NaN ):( ($2 < T_t) ? $4 : NaN ) with lines title "Value (transient)" lc rgb "blue", \
##        filename every ::0::i using ( ($2 >= T_t) ? $2 : NaN ):( ($2 >= T_t) ? $4 : NaN ) with lines title "Value (steady)" lc rgb "magenta", \
#
#    pause 0.01
#}


T_t = real(word(system(sprintf("sed -n '%dp' %s", maxLines, filename)), 3))
set arrow from T_t, graph 0 to T_t, graph 1 nohead lc rgb "black" lw 2
plot \
    filename using 2:5 with lines title "Mean" lc rgb "green", \
    filename using 2:4 w l t "Filter" lc rgb "red", \
    filename using 2:7:6 w filledcu t "Confidence Interval", \
    filename using 2:8 w l t "Initial" lc rgb "black", \
    #filename using ( ($2 < T_t) ? $2 : NaN ):( ($2 < T_t) ? $4 : NaN ) with lines title "Value (transient)" lc rgb "blue", \
    #filename using ( ($2 >= T_t) ? $2 : NaN ):( ($2 >= T_t) ? $4 : NaN ) with lines title "Value (steady)" lc rgb "magenta", \

    pause 1000


