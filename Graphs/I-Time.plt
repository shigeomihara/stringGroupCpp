set lmargin 9
set bmargin 7
#set xlabel "{/Arial:Roman=20 Time (month/day)}" offset 0,-1 font "Arial,12"
set xlabel "{/Arial:Roman=20 Time (hour:minute)}" offset 0,-1 font "Arial,12"
set ylabel "{/Arial:Roman=20 {/Arial:Italic=20 I} (A)}" offset 1.2 font "Arial,12"
#set ylabel "{/Arial:Roman=20 {/Arial:Italic=20 V} (V)}" offset 1.5 font "Arial,12"
#set ylabel "{/Arial:Roman=20 {/Arial:Italic=20 G} (kW/m^2)}" offset 1.5 font "Arial,12"
#set y2label "{/Arial:Roman=20 {/Arial:Italic=20 T} ({/Symbol \260}C)}" offset 1.5 font "Arial,12"

zyougen = 1.58400
#zyougen = 6.0
kagen = 0.00000

set grid
set xdata time

set timefmt "%Y/%m/%d/%H:%M"
#set xrange ["2019/5/1":"2019/9/30"]
#set xrange ["2019/5/7/9:0":"2019/5/7/16:0"]
#set yrange [kagen:zyougen]

#set timefmt "%Y/%m/%d/%H:%M"
#set format x "%m/%d"
set format x "%H:%M"
set xtics rotate by -90 font "Arial,11"
set ytics font "Arial,11"
#set y2tics font "Arial,11"

set terminal png enhanced font "Arial,10"
set output "I-TimeSim.png"

set datafile separator comma
#set datafile separator ", "

#plot "GTIVTimeSim.dat" using 5:1 with linespoints pt 7 lc "red" ps 0.5 notitle
plot "../Dat/GTIVTimeSim.dat" using 5:($3) with linespoints pt 7 lc "blue" ps 0.5 axes x1y1 title "{/Arial:Italic=20 Vmeas}" at 0.8, 0.9,\
 "../Dat/GTIVTimeSim.dat" using 5:($6) with linespoints pt 7 lc "red" ps 0.5 axes x1y2 title "{/Arial:Italic=20 Vsim}" at 0.8, 0.83
#plot "PstcTime0.dat" using 1:2 with linespoints pt 7 lc "blue" ps 0.5 notitle
