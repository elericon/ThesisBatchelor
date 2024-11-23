set datafile sep '\t'
set grid
set logscale x
set ytics add ("1 " 1)

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 1 pi -1 ps 1.5
set term qt 0 title "blocking Magnetizzazione" font "Helvetica" position 50,0
set title "Magnetization error Blocking" font ",16"
set xlabel "B" font ",14"
set ylabel "{/Symbol s}" font ",14"
plot "out/blocking.txt" i 0 u 1:2 w lp ls 1 title ""

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 1 pi -1 ps 1.5
set term qt 3 title "jackknife Magnetizzazione" font "Helvetica" position 750,0
set title "Magnetization error Jackknife" font ",16"
set xlabel "B" font ",14"
set ylabel "{/Symbol s}" font ",14"
plot "out/jackknife.txt" i 0 u 1:2 w lp ls 1 title ""

pause-1

