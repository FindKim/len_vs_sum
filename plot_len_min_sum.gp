#
# Histogram of the %min sum of rare codon clusters divided by bin length
# Plots to histograms: results for significant sequences and non-sig
#
#


clear
reset

set terminal postscript
set output "results/sigOrfs_len_min_sum_p1e-05.ps"

set xtics nomirror out
set ytics nomirror out
set xrange [0:300]
set yrange [-20000:0]
#set auto x
#set auto y

set xlabel "Rare Codon Cluster Length"
set ylabel "Rare Codon Cluster %Min Sum"

set style fill solid border -1
#set style data histogram
set datafile separator ","

set key off
get_point(x,y) = sprintf("(%.0f,%.2f)", x, y)
get_point(x,y) = sprintf("(%.0f,%.2f)", x, y)

set title "rare codon cluster: length against %min sum\nusing significant sequences\nwith p-Value 1e-05"
plot 'results/sig_len_min_sum_p1e-05.txt' u 1:2

#set auto x
#set xrange [0:4000]
#set yrange [0:100]

set title "rare codon cluster: length against %min sum\nusing non-significant sequences\nwith p-Value 1e-05"
plot 'results/nsig_len_min_sum_p1e-05.txt' u 1:2

#pause -1
