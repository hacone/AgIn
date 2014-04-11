set grid
set terminal png size 6000, 600
set output "segment-arg1.png"
set multiplot
set pointsize 0.3

set lmargin at screen 0.01
set nokey

set xrange [arg1 : arg1 + arg2]
set xtics 100000 format ""
set mxtics 20

set size 1.0, 0.7
set origin 0.0, 0.3
set yrange [-1.2:1.2]

plot "./BisulfiteSegment7.dat" u 1:(0.08*($3)+0.04) title "Bisulfite7" pt 1 lc rgb "magenta",\
     "./IPDSegment.dat" u 1:(-0.08*($3)-0.04) title "Prediction-5" pt 1 lc rgb "blue",\
     "./BisulfiteSegment7.dat" u 1:($2+0.1) title "site mC level" pt 1 lc -1,\
     "./IPDSegment.dat" u 1:(-1.0*($2)) title "site-wise prediction from IPDs" pt 1 lc -1

set tmargin 0
set logscale y
set yrange [0.001:1000]
set size 1.0, 0.3
set origin 0.0, 0.0
set xtics format "%.1s %cbp"

plot "./IPDConfval.dat" u 1:2 title "liability of segmentation" w p lw 3 lc rgb "blue", 1
