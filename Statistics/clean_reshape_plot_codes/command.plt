set xlabel "time(t) in units of tau"
set ylabel "<r^2>(t)"
set term png
set output "100_<r_2>(t).png"
plot "r_sq_samp_av.txt" using 1:2 with lines
set term x11
