set datafile sep '\t'
set grid

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5
set term qt 0 title "<M>" font "Helvetica" position 50,0
set ytics add ("1 " 1)
set ylabel "<M>" font ",14"
set xlabel "Numero iterazioni" font ",14"
plot "out/osservabili.txt" i 0 u 1:2 w l lw 2 title ""

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5
set term qt 1 title "chi" font "Helvetica" position 450,0
set ytics add ("1 " 1)
set ylabel "{/Symbol c}" font ",14"
set xlabel "Numero iterazioni" font ",14"
plot "out/osservabili.txt" i 0 u 1:4 w l lw 2 title ""

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5
set term qt 2 title "binder ratio" font "Helvetica" position 750,0
set ytics add ("1 " 1)
set ylabel "U" font ",14"
set xlabel "Numero iterazioni" font ",14"
plot "out/osservabili.txt" i 0 u 1:6 w l lw 2 title ""

set term qt 3 title "Dimensione cluster - Temperatura" font "Helvetica" position 50,350
set xlabel "Numero iterazioni" font ",14"
set ylabel "<Dimensione>" font ",14"
set ytics add ("1 " 1)
plot "out/dati.txt" u 1:2 w l lw 2 title ""

set term qt 4 title "calore specifico - Temperatura" font "Helvetica" position 450,350
set xlabel "Numero iterazioni" font ",14"
set ylabel "cv" font ",14"
set ytics add ("1 " 1)
plot "out/dati.txt" u 1:4 w l lw 2 title ""


pause-1