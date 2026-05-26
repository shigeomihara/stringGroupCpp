set lmargin 9
set bmargin 7
set xlabel "{/Arial:Roman=20 Time (month/day)}" offset 0,-1 font "Arial,12"
set ylabel "{/Arial:Roman=20 {/Arial:Italic=20 G} (kW/m^2)}" offset 1.5 font "Arial,12"

zyougen = 1.58400
#zyougen = 6.0
kagen = 0.00000

set grid
set xdata time

set timefmt "%Y/%m/%d"
#set xrange ["2019/5/1":"2019/9/30"]
set xrange ["2019/5/1":"2019/5/10"]
set yrange [kagen:zyougen]

set timefmt "%Y/%m/%d/%H:%M"
set format x "%m/%d"
set xtics rotate by -90 font "Arial,11"
set ytics font "Arial,11"

set terminal png enhanced font "Arial,10"
set output "GTime.png"

#set datafile separator "whitespace|comma"

plot "../GTIVTime.dat" using 10:2 with linespoints pt 7 lc "blue" ps 0.5 notitle
#plot "PstcTime0.dat" using 1:2 with linespoints pt 7 lc "blue" ps 0.5 notitle
