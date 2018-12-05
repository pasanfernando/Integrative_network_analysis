import networkx as nx
import matplotlib.pyplot as plt

"""  code to gene function prediction based on the higishaki methods
    This code is not the full method, but a method modified by Pasan Fernando to focus on one function only"""

############################################################################# ###########

#input1 = raw_input('enter the interaction file:')
#input2 = raw_input('enter the genelist file:')

### reading the interaction file########################################################
in1 = open('network.txt', 'r')
G = nx.Graph()
for line in in1:
    if line != '\n':
        #print line
        if 'combined_score' not in line:
            line = line.strip('\n')
            a = line.split()
            #print a[0],a[1],a[2]
            G.add_edge(a[0].lower(),a[1].lower(), weight = a[2])

# visualizing the network; please turn this off when loading large networks
# nx.draw_networkx(G, with_labels=True)
# plt.draw()
# plt.show()

# reading the gene list file
genelist = []
in2 = open('genelist.txt', 'r')

for line in in2:
    if line != '\n':
        line = line.strip()
        a = line.split(',')
        for i in a:
            i = i.strip()
            genelist.append(i.lower())
print genelist

######### calculate the function frequency or phi #########
#total proteins in the network
totps = nx.number_of_nodes(G)

#total proteins for the given function
#funps = len(genelist) # you are assuming all genes matches the ones in the network need to check one by one


genes_in_network = nx.nodes(G)
print genes_in_network
#making sure all the genes in the input data file are in the network

notfound = set(genelist)-set(genes_in_network) # genes that are not found in the network
print 'genes that are not found in the network:',notfound
found = set(genelist) & set(genes_in_network) # genes that are found in the network
print 'genes that are found in the network:',found

funps = len(found)


fracps = float(funps)/totps

print fracps
#dictionaries to store total and function counts, chi square value
total={}
ccount={}
chisq = {}

#counting the genes with given function in immediate neighborhood of unannotated genes
for i in G.nodes():
    if i not in genelist:
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

        #calculate expected frequency values
        ef = fracps*totalc
        chi = ((ccountc-ef)**2)/ef
        chisq[i] = chi

print chisq

chilist = chisq.values()

print chilist

plt.hist(chilist, color='b', histtype='bar', bins=50, normed=True)
plt.show()