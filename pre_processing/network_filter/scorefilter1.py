#!/usr/bin/python

# filters the given network by a given score cutoff

first  = raw_input('Enter Input file: ')
second = raw_input('Enter output file: ')
value = raw_input('Enter selection value: ')
value1= float(value)

out = open(second, 'wb+')

f= open(first, 'r')

for line in f:
	#print line
	splitLine = line.split()
	if splitLine[2]== 'combined_score':
		out.write (line)
	#print splitLine[2]
	#print '\n'
	else:
		#print splitLine[2]
		con = float(splitLine[2])
		if con >= value1 :
			out.write (line)