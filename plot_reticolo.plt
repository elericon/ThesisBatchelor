set datafile sep '\t'
set grid


set style line 2 lc rgb 'black' pt 7 ps 1
set term qt 0 title "reticolo" font "Helvetica" position 50,0
set ytics add ("1 " 1)
set xrange [0:14]
set yrange [0:14]
plot "out/reticolo.txt" i 0 u 1:2 w p title "" ls 2

pause-1