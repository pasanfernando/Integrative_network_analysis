
'''This code counts the neighborhood genes with a given phenotype. Used to compare the modularity of a given module
Compare it with the network background
This will generate two output files, one for the candidates and one for the non candidates detailing their neghborhood counts

Assumption: all genes in the gene list should be in the network. The genelist input should be reconciled.
'''

import networkx as nx
#import matplotlib.pyplot as plt

############################################################################# ###########
### reading the interaction file########################################################
#input1 = raw_input('enter the interaction file:')
#input2 = raw_input('enter the genelist file:')

in1 = open('mslick0.3integrated_network.txt', 'r')
G = nx.Graph()
for line in in1:
    if line != '\n':
        #print line
        if 'combined_score' not in line:
            line = line.strip('\n')
            a = line.split()
            #print a[0],a[1],a[2]
            G.add_edge(a[0].lower(),a[1].lower(), weight = a[2])


#print G.edges(data=True)

genelist = []
in2 = open('modulegenes.txt', 'r')

for line in in2:
    if line != '\n':
        line = line.strip()
        genelist.append(line.lower())
#print genelist

total={}
ccount={}

# Checking how many candidate genes that are not in the network
notinnet = set(genelist)- set(G.nodes())

print 'how many candidate genes that are not in the network:',len(notinnet)

#####################################################################################
#### this code block is for candidate genes. Comment out if used for noncandidates
# for i in genelist:
#     if G.has_node(i):
#         totalc =0
#         ccountc=0
#         nb = nx.all_neighbors(G, i)
#         #print i
#         for j in nb:
#             #print j
#             totalc=totalc+1
#             if j in genelist:
#                 ccountc=ccountc+1
#
#         total[i]=totalc
#         ccount[i]=ccountc
#
# for i in genelist:
#     if G.has_node(i):
#         out.write('%s %s %s %f\n' % (i, ccount[i], total[i], float(ccount[i]) / float(total[i])))
############################################################################
#print total
#print ccount
##### Counts the candidate genes in the neighborhood for all the genes in the network including candidates and noncandidates

## please comment this out if used for candidate genes
for i in G.nodes():

    totalc =0
    ccountc=0
    nb = nx.all_neighbors(G, i)
    #print i
    for j in nb:
        #print j
        totalc=totalc+1
        if j in genelist:
            ccountc=ccountc+1

    total[i]=totalc
    ccount[i]=ccountc

# writing the output files
# For noncandidates
noncan = open('noncandidatesneighbourcount.txt', 'wb+')
noncan.write('gene_name candidate_count total_count candidate_percentage\n')

#for candidates

can = open('candidatesneighbourcount.txt', 'wb+')
can.write('gene_name candidate_count total_count candidate_percentage\n')

# For noncandidates
for i in G.nodes():
    if i in genelist:
        can.write('%s\t%s\t%s\t%f\n' % (i, ccount[i], total[i], float(ccount[i]) / float(total[i])))

    else:
        noncan.write('%s\t%s\t%s\t%f\n' % (i, ccount[i], total[i], float(ccount[i]) / float(total[i])))

noncan.close()
can.close()
##############################################################################



