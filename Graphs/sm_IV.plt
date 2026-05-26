#set bmargin at screen 0.3

#set xlabel "{/Arial:Italic=20 {/Symbol l} (nm)}" offset -0.5
#set xlabel "{/Arial:Italic=20 x ({/Symbol m}m)}" offset -0.5
set xlabel "{/Arial:Italic=20 V} {/Arial:Roman=20 (V)}" offset -0.5
#set xlabel "{/Arial:Italic=20 n (cm^{-3})}" offset 0, graph -0.1

#set ylabel "{/Arial:Italic=20 {/Symbol j}} {/Arial:Roman=20 (V)}" offset 0
#set ylabel "{/Arial:Italic=20 E} {/Arial:Roman=20 (V/m)}" offset 0
set ylabel "{/Arial:Italic=20 I} {/Arial:Roman=20 (A)}" offset 0
#set ylabel "{/Arial:Italic=20 n} {/Arial:Roman=20 (/cm^3)}" offset 0
#set ylabel "{/Arial:Italic=20 p} {/Arial:Roman=20 (/cm^3)}" offset 0
#set ylabel "{/Arial:Italic=20 J_n} {/Arial:Roman=20 (A/cm^2)}" offset 0
#set ylabel "{/Arial:Italic=20 J_p} {/Arial:Roman=20 (A/cm^2)}" offset 0
#set ylabel "{/Arial:Italic=20 I_0 {/Arial:Roman=20 (s^{-1}m^{-3})}" offset 0.5
#set ylabel "{/Arial:Italic=20 I_0 {/Arial:Roman=20 (s^{-1}cm^{-2}nm^{-1})}" offset 0.5
#set ylabel "{/Arial:Italic=20 G {/Arial:Roman=20 (s^{-1}cm^{-3})}" offset 0.5
#set ylabel "{/Arial:Italic=20 U {/Arial:Roman=20 (s^{-1}cm^{-3})}" offset 0.5

set colorsequence classic
set grid
#set logscale xy

set xzeroaxis
set yzeroaxis

set yrange [-1:12]

#set xtics rotate by 90 offset 0, graph -0.25

#set key left top
#set key right top

set terminal png enhanced font "Arial,16"

#set output "result_phi.png"
set output "sm_IV.png"
#set output "result_n.png"
#set output "result_p.png"
#set output "result_Jn.png"
#set output "result_Jp.png"
#set output "generation_alpha.png"
#set output "generation_I0.png"
#set output "generation_G.png"
#set output "recombination_U-n.png"

#plot "result.dat" using ($2*1.0E6):($8) with linespoints title "{/Arial:Italic=20 {/Symbol j}}"
plot "sm_IV.dat" using ($2):($4) with linespoints notitle
#plot "result.dat" using ($2*1.0E6):($10*1.0E-6) with linespoints title "n"
#plot "result.dat" using ($2*1.0E6):($12*1.0E-6) with linespoints title "p"
#plot "result.dat" using ($2*1.0E6):($4*1.0E-4) with linespoints title "J_n"
#plot "result.dat" using ($2*1.0E6):($6*1.0E-4) with linespoints title "J_p"
#plot "generation_alpha.dat" using ($2):($4) with linespoints notitle
#plot "generation_I0.dat" using ($2):($4*1.0E-13) with linespoints notitle
#plot "generation_G.dat" using ($2*1.0E6):($4*1.0E-6) with linespoints notitle
# plot "recombination_U-n_pE17.dat" using ($2):($4) with linespoints title "p=1.0E17",\
# "recombination_U-n_pE20.dat" using ($2):($4) with linespoints title "p=1.0E20",\
# "recombination_U-n_pE23.dat" using ($2):($4) with linespoints title "p=1.0E23"
