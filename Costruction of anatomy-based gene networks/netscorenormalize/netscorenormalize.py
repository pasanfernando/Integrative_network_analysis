

'''
Normalize the score distribution of a network using min-max normalization to range it up between 0 and 1
use this to normalize networks like resnik
'''

#Opening the network to store the min and max of the distribution

netfile = open('zebrafulldupremoved_network.txt','r')

#list to store the score distribution
scordist =[]

for line in netfile:
    if line != '\n' and 'combined_score' not in line:
        a = line.strip().split(' ')
        scordist.append(float(a[2]))

netfile.close()

# getting the min and max
minval = min(scordist)
maxval = max(scordist)
print 'min value:',minval
print 'max value:',maxval

# Opening the original network

in1 = open('zebrafulldupremoved_network.txt', 'r')
# opening the preprocessed network output file
pre_net = open('normzebrafulldupremoved_network.txt', 'wb+')

# writing the header line
pre_net.write('protein1 protein2 combined_score\n')

for line in in1:
    if line != '\n':
        # print line
        if 'combined_score' not in line:
            templist = []
            line = line.strip('\n')
            a = line.split()
            # print a[0],a[1],a[2]
            gene1 = a[0].lower()
            # print gene1
            gene2 = a[1].lower()
            normscore = (float(a[2])-minval)/(maxval-minval)
            pre_net.write('%s %s %s\n' % (gene1, gene2, normscore))


in1.close()
pre_net.close()

