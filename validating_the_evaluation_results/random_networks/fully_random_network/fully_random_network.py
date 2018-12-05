

'''
Generates a completely random network from given number of nodes and edges
Replaces the nodes by the gene names in the profile

This code reads the original network, extracts the names and generates a random network by the same names


The number of nodes must exactly match the number of genes in the original network; please double check
'''


import networkx as nx
import random


# reading the original network and populating the geneset
netfile = open('wangstring_network.txt', 'r')

geneset=[]
G = nx.Graph()
for line in netfile:
    if line != '\n':
        #print line
        if 'combined_score' not in line:
            line = line.strip('\n')
            a = line.split()
            geneset.append(a[0].lower())
            geneset.append(a[1].lower())
            G.add_edge(a[0].lower(), a[1].lower(), weight=float(a[2]))

# extracting the number of nodes and number of edges from the original network
numnodes= nx.number_of_nodes(G)
numedges= nx.number_of_edges(G)

print 'number of nodes:',numnodes
print 'number of edges:',numedges
# removing duplicate genes if there are any
geneset=list(set(geneset))

print 'number of genes in the profile:',len(geneset)

G = nx.gnm_random_graph(numnodes,numedges)


#geneset = ['a','b','c','d','f']

genemap ={}

for i in nx.nodes(G):
    genemap[i]=geneset[i]

nodes = nx.nodes(G)
print nodes


# if you are using the built in relabeling you are duplicating the task
# re-labeling the nodes
#
# G = nx.relabel_nodes(G,genemap)
#
# nodes = nx.nodes(G)
# print nodes
# #

# generating all the edges in the re-labeled network
edges = nx.edges(G)
#print edges


# writing the randomized network
ran = open('randomwangstring_network.txt','wb+')

ran.write('protein1 protein2 combined_score\n')

counter=0
for i in edges:
    counter=counter+1
    print counter
    # converting the set to list so it can be indexed
    listi= list(i)
    #print listi
    ran.write('%s %s %s\n'%(genemap[listi[0]],genemap[listi[1]],random.random()))



