#!/usr/bin/python

import os
import re

ns_pattern = re.compile('nsig_len_min_sum_p.+\.txt')
s_pattern = re.compile('sig_len_min_sum_p.+\.txt')

LENGTH,SUM = range(2)	# 0, 1

len_out = open('/afs/crc.nd.edu/user/k/kngo/RCClust/results/avg_rcclust_length.txt', 'w') # overwrites if exits

sum_out = open('/afs/crc.nd.edu/user/k/kngo/RCClust/results/avg_rcclust_min_sum.txt', 'w') # overwrites if exits

len_out.write('Calculated average length for rare codon clusters differing by p-Value cut offs.')
len_out.write('\n')
len_out.write('\n')

sum_out.write('Calculated average length for rare codon clusters differing by p-Value cut offs.')
sum_out.write('\n')
sum_out.write('\n')

for path, subdirs, files in os.walk("/afs/crc.nd.edu/user/k/kngo/RCClust/results/"):

	for filename in files:
		print filename
		len_out.write(filename)
		len_out.write('\n')
		sum_out.write(filename)
		sum_out.write('\n')
		
		# Avg variables
		s_rcclust_len = 0
		s_rcclust_sum = 0
		s_rcclust_count = 0
		ns_rcclust_len = 0
		ns_rcclust_sum = 0
		ns_rcclust_count = 0
		
		# Nonsigificant files
		if ns_pattern.match(filename):

			file_path = os.path.join(path, filename)
			with open(file_path) as f:
				header = f.readline()
				for line in f:
					line = line.strip()		# Removes newline char
					line = line.split(",")	# Splits floato array of 3
					ns_rcclust_sum += float(line[SUM])
					ns_rcclust_len += float(line[LENGTH])
					ns_rcclust_count += 1

			ns_rcclust_len_avg = ns_rcclust_len / ns_rcclust_count
			ns_rcclust_sum_avg = ns_rcclust_sum / ns_rcclust_count
			print (ns_rcclust_count)
			len_out.write('nsig avg length: ' + str(round(ns_rcclust_len_avg,2)))
			len_out.write('\n')
			sum_out.write('nsig avg sum: ' + str(round(ns_rcclust_sum_avg,2)))
			sum_out.write('\n')
			
			
		# Significant files
		elif s_pattern.match(filename):

			file_path = os.path.join(path, filename)
			with open(file_path) as f:
				header = f.readline()
				for line in f:
					line = line.strip()		# Removes newline char
					line = line.split(",")	# Splits floato array of 3

					s_rcclust_sum += float(line[SUM])
					s_rcclust_len += float(line[LENGTH])
					s_rcclust_count += 1

			s_rcclust_len_avg = s_rcclust_len / s_rcclust_count
			s_rcclust_sum_avg = s_rcclust_sum / s_rcclust_count	
			print (s_rcclust_count)

			len_out.write('sig avg length: ' + str(round(s_rcclust_len_avg,2)))
			len_out.write('\n')
			sum_out.write('sig avg sum: ' + str(round(s_rcclust_sum_avg,2)))
			sum_out.write('\n')
