# Plot seismogram data
#=====================

set xlabel 'Time in s'
set ylabel 'Amplitutude'

set xrange[0:2]
set yrange[-0.05:0.05]

set output 'seismogram.png'
set terminal pngcairo enhanced font 'Verdana,14' size 1920,1080

set key off # unset key

set xtics 0.2
set ytics 0.01

set style line 1 linecolor rgbcolor '#009474' # grün türkis

set title 'Seismogram at [70,70,70]'
plot 'seismogram.mtx' every::2 using ($0*0.002):1 with lines ls 1# plots 1 column of seismogram.mtx skipping the first two lines
