#!/usr/bin/python

import os

os.system("sed -i 's/_p1e-05/_p\.0001/g' plot_len_min_sum.gp")
os.system("sed -i 's/p-Value .*\"/p-Value \.0001\"/g' plot_len_min_sum.gp")
os.system("gnuplot plot_len_min_sum.gp")
os.system("ps2pdf results/sigOrfs_len_min_sum_p.0001.ps sigOrfs_len_min_sum_p.0001.pdf")

os.system("sed -i 's/_p\.0001/_p\.001/g' plot_len_min_sum.gp")
os.system("sed -i 's/p-Value .*\"/p-Value \.001\"/g' plot_len_min_sum.gp")
os.system("gnuplot plot_len_min_sum.gp")
os.system("ps2pdf results/sigOrfs_len_min_sum_p.001.ps results/sigOrfs_len_min_sum_p.001.pdf")

os.system("sed -i 's/_p\.001/_p\.01/g' plot_len_min_sum.gp")
os.system("sed -i 's/p-Value .*\"/p-Value \.01\"/g' plot_len_min_sum.gp")
os.system("gnuplot plot_len_min_sum.gp")
os.system("ps2pdf results/sigOrfs_len_min_sum_p.01.ps results/sigOrfs_len_min_sum_p.01.pdf")

os.system("sed -i 's/_p\.01/_p\.05/g' plot_len_min_sum.gp")
os.system("sed -i 's/p-Value .*\"/p-Value \.05\"/g' plot_len_min_sum.gp")
os.system("gnuplot plot_len_min_sum.gp")
os.system("ps2pdf results/sigOrfs_len_min_sum_p.05.ps results/sigOrfs_len_min_sum_p.05.pdf")

os.system("sed -i 's/_p\.05/_p1e-05/g' plot_len_min_sum.gp")
os.system("sed -i 's/p-Value .*\"/p-Value 1e-05\"/g' plot_len_min_sum.gp")
os.system("gnuplot plot_len_min_sum.gp")
os.system("ps2pdf results/sigOrfs_len_min_sum_p1e-05.ps results/sigOrfs_len_min_sum_p1e-05.pdf")

os.system("rm results/*.ps")
