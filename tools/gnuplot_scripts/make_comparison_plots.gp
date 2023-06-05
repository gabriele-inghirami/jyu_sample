set term pos eps enh col font "Helvetica, 20"

set xrange [0:3.8]
set key top right


do for [i=1:318] {
id_th=sprintf("%03d", i )
id_mc=sprintf("%03d", i+1)

file_th="reference/".id_th.".dat"
file_mc_wd="fias_wdec_partid_".id_mc.".dat"

t1="Th.dec"
t4="MC dec"


set out "plots/dN_".id_th.".eps"
set ylabel "<dN/(dy 2 pT dpT)>" 
#set ylabel "<dN/(2 pT dpT)>" 
plot file_th u 1:7 w l lw 2 lc "black" dt 1 t t1, file_mc_wd u 1:2 w p pt 7 lc "blue" t t4

#set logscale y
#set out "plots/dN_".id_th."_logscale.eps"
#plot file_th u 1:7 w l lw 2 lc "black" dt 1 t t1, file_mc_wd u 1:2 w p pt 7 lc "blue" t t4

set out "plots/v2_".id_th.".eps"
set ylabel "v_2"
plot file_th u 1:8 w l lw 2 lc "black" dt 1 t t1, file_mc_wd u 1:3 w p pt 5 lc "magenta" t t4

set out "plots/v3_".id_th.".eps"
set ylabel "v_3"
plot file_th u 1:9 w l lw 2 lc "black" dt 1 t t1, file_mc_wd u 1:4 w p pt 5 lc "magenta" t t4


unset logscale
}
