set datafile sep '\t'
set grid
set key font ",13"
set ytics add ("1 " 1)

N = 3

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5

set term qt 0 title "magnetizzazione media - Temperatura" font "Helvetica" position 50,0
set xlabel "T" font ",14"
set ylabel "<M>" font ",14"
plot for [k=0:N] "out/risultati.txt" index k u 1:2:3 with yerrorlines lt k+1 title columnheader(1)

set term qt 1 title "chi - Temperatura" font "Helvetica" position 350,0
set xlabel "T" font ",14"
set ylabel "{/Symbol c}" font ",14"
#set xrange [0.5:3]
set xrange [1:4]
plot for [k=0:N] "out/risultati.txt" index k u 1:4:5 with yerrorlines lt k+1 title columnheader(1)

set term qt 2 title "calore specifico (per spin) - Temperatura" font "Helvetica" position 750,350
set xlabel "T" font ",14"
set ylabel "c_v per spin" font ",14"
plot for [k=0:N] "out/risultati.txt" index k u 1:9:10 with yerrorlines lt k+1 title columnheader(1)

set term qt 3 title "Binder ratio - Temperatura" font "Helvetica" position 750,0
set xlabel "T" font ",14"
set ylabel "U" font ",14"
#set xrange [0.5:3]
set xrange [2.15:2.4]
plot for [k=0:N] "out/risultati.txt" index k u 1:6:7 with yerrorlines lt k+1 title columnheader(1)
unset xrange

set term qt 4 title "dimensione cluster normalizzata - Temperatura" font "Helvetica" position 350,350
set xlabel "T" font ",14"
set ylabel "Mean dimension per spin" font ",14"
plot for [k=0:N] "out/risultati.txt" index k u 1:11:12 with yerrorlines lt k+1 title columnheader(1)

set term qt 5 title "chi - Temperatura" font "Helvetica" position 50,350
set xlabel "t*L" font ",14"
set ylabel "{/Symbol c} scaled" font ",14"
set xrange [-6:12]
plot for [k=0:N] "out/results_scaled.txt" index k u 1:2:3 with yerrorbars lt k+1 title columnheader(1)
unset xrange

pause-1