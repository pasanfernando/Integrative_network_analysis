m = raw_input('Enter input file: ')
ou = raw_input('Enter output file: ')

dict={}
f= open(m, 'r') # open the input file and store it in a hash
for line in f:
	split = line.split()
	x1= split[0]
	x2= split[1]
	x3= split[2]
	dict[x1,x2]= (x3)
	#print line

f.close()
f= open(m, 'r')
out= open(ou, 'wb+')

for line in f:
	split = line.split()
	x1= split[0]
	x2= split[1]
	x3= split[2]
	if (x1,x2) in dict and (x2,x1) in dict:
		print line
		out.write(line)
		print 'yes'
		if x1==x2:
			del dict[(x1, x2)]
		else:
			del dict[(x1,x2)]
			del dict[(x2,x1)]
	
	elif (x1,x2) in dict:
		#print line
		out.write(line)
		del dict[(x1,x2)]
		
	
	
print dict

			