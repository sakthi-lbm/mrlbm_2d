reset
set terminal postscript eps color enhanced "Times-Roman, 22"
set output 'p_mean_comp_smooth.eps' 
set size 1,1.2
#set key box opaque
set ylabel "C_p" offset 3,0,0
set xlabel "{/Symbol q} (deg)"
plot "pmean.dat" u 1:5 with l lw 3 lc -1 title 'Present (velocity interp)',\
"lit_hwang_2007.dat" u 1:2 w p pt 7 title 'Hwang and Yang, 2007',\
"lit_vu_2015.dat" u 1:2 pt 5 w p title 'Vu et. al, 2015'

reset
set terminal postscript eps color enhanced "Times-Roman, 22"
set output 'p_mean_comp.eps' 
set size 1,1.2
#set key box opaque
set ylabel "C_p" offset 3,0,0
set xlabel "{/Symbol q} (deg)"
plot "pmean.dat" u 1:4 with l lw 3 lc -1 title 'Present (velocity interp)',\
"lit_hwang_2007.dat" u 1:2 w p pt 7 title 'Hwang and Yang, 2007',\
"lit_vu_2015.dat" u 1:2 pt 5 w p title 'Vu et. al, 2015'
