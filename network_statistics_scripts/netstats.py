

'''
Prints the number of genes and the number of interactions in any input nertwork
'''



# Opening the original network

in1 = open('newlinstring0.33preprocessed.txt', 'r')

genelist =[]

# a list for interactions
intlist =[]

# list for duplicated interactions
duplicated =[]
for line in in1:
    if line != '\n':
        # print line
        if 'combined_score' not in line:

            line = line.strip('\n')
            a = line.split()
            # print a[0],a[1],a[2]
            gene1 = a[0].lower()
            # print gene1
            gene2 = a[1].lower()
            genelist.append(gene1)
            genelist.append(gene2)
            # if gene1+gene2 in intlist:
            #     duplicated.append(gene1+gene2)
            # elif gene2+gene1 in intlist:
            #     duplicated.append(gene2+gene1)
            # else:
            intlist.append(gene1+gene2)


print 'number of genes in the network:',len(set(genelist))

print 'number of interactions in the network',len(intlist)
# print 'number of duplicated interactions in the network',len(duplicated)
# print duplicated
