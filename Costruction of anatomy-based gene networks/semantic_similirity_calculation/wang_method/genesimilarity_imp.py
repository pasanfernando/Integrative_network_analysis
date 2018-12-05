import networkx as nx
from collections import defaultdict

import wang as ws

'''This method use wang gene similarity method to calculate the similarity between two genes'''

###############################################################################################################
def sumcalculate(l1,l2):
    maxsum =[]
    for i in l1:
        list1 = []
        for j in l2:
            #print i+j
            if i+j in dytable:
                list1.append(dytable[i+j])
                #print '1'
            elif j+i in dytable:
                list1.append(dytable[j+i])
                #print '2'
            else:
                simvalue = ws.wangsim(i,j)
                dytable[i+j] = simvalue
                list1.append(simvalue)

                #print simvalue

        #print list1
        maxsum.append(max(list1))
        #print maxsum
    #print 'the sumis',sum(maxsum)
    #print len(maxsum)
    return sum(maxsum),len(maxsum)
#####################################################################################################



# ifile = open('geneprofiles.txt', 'r')
#
#
# gp = collections.defaultdict(list)
#
# # reading the input file and storing in multi dictionary as gene profiles
# for line in ifile:
#     if line != '\n':
#         x = line.split()
#         if x[0]!= 'gene_symbol':
#             for i in range(1,len(x)):
#                 gp[x[0]].append(x[i])

# opening the annotation file
genelist = open('zebra_randomfuncfiltered_profiles.txt', 'r')

# defining a multiple dictionary to store functions: genes
genefunc = defaultdict(list)
# defining a multiple dictionary to store gene: functions
genedic = defaultdict(list)

# Then convert the above input file into function: gene format

# this indicator is required for leave one out method

# reading the gene file
for line in genelist:
    if line != '\n':
        line = line.strip()
        a = line.split('\t') # splitting the line by tab to separate gene name vs functions
        #print a
        if a[1] != 'gene_name': # excluding the header
            a[1] = a[1].lower()  # convert to lower case to avoid mismatches
            # print a[0]
            # here, we are evaluating, so only the proteins with functions will be considered
            if len(a) > 1:  # ths is required because there are proteins without one single function
                b = a[3].split(',')  # seperating the functions by comma
                # store each function in the multi dic by function as the key and genes as values
                for i in b:
                    genedic[a[1]].append(i)  # this dictionary stores genes as keys
                    genefunc[i].append(a[1])

print genedic


#print gp
#print gp['g7']

def wanggenesim(g1,g2):
    gene1 = genedic[g1]
    #print gene1

    gene2 = genedic[g2]
    #print gene2

    sumgene1,m = sumcalculate(gene1,gene2)
    # print sumgene1
    # print m

    sumgene2,n = sumcalculate(gene2,gene1)
    # print sumgene2
    # print n

    genesimilarity = (sumgene1+sumgene2)/(m+n)
    return genesimilarity

dytable = {}

#
# a = raw_input('enter first gene:')
# b = raw_input('enter second gene:')
#
# sim = wanggenesim(a,b)
# print 'the similarity betweeen the two genes is:',sim
#print dytable

# a dictionary to store pairwise score
pairwise ={}

out = open('wang_Rfilterednetwork.txt','wb+')
out.write('protein1 protein2 combined_score\n')

counter =0
# calculating the pairwise gene similarity
for i in genedic:
    for j in genedic:
        # do not calculate the similarity betweeen the same gene
        if i != j:
            if i+j not in pairwise and j+i not in pairwise:
                counter = counter +1
                print counter
                sim = wanggenesim(i, j)
                pairwise[i+j]=sim
                # only writing scores that are higher than 0
                if sim> 0:
                    out.write('%s %s %s\n'%(i,j,sim))





