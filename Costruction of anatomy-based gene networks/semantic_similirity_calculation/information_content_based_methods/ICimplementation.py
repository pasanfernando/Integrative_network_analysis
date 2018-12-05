import networkx as nx
import operator
import math
import re
from collections import defaultdict


####### reading Uberon #########################


ub = nx.DiGraph()
uberon = {}
name = {}

#reading the uberon file
p = open('uberon.obo', 'r')

for line in p:
    if line.startswith('[Typedef]'):
        break

    if line.startswith('id:'):
        x = line.split()
        y = x[1].replace(':', '_').lower()

        if ub.has_node(y) == False:
            ub.add_node(y)

    if line.startswith('name:'):
        k = re.sub('name: ','',line,1)
        z = k.strip()
        d = z.replace(" ","_")
        uberon[y] = d
        name[z] = y

    if line.startswith('is_a: '):
        z = line.split()
        k = z[1].replace(':', '_').lower()
        # print k
        ub.add_edge(k, y)

    # elif line.startswith('relationship: ') or line.startswith('intersection_of:'):
    #     if 'UBERON:' in line:
    #         line = line.strip()
    #         m = 'UBERON:'+re.search('UBERON:(.+?) ', line).group(1)
    #         ub.add_edge(m, y, color = 'red')


#G = nx.drawing.nx_agraph.read_dot('name1.dot')
# assingning uberon graph to G
G= ub
# for n,d in G.in_degree().items():
#     if d==0:
#         print 'root node is',n
################################################


# nx.write_dot(G,'name.dot')
# count ={}
# fcount ={}
# dictionary to store information content
ic ={}
#G = nx.read_dot('name1.dot')
# print nx.nodes(G)
# k = nx.descendants(G, 'B')
# a = nx.ancestors(G, 'H')
# print k
# print a

#Reading the input gene function annotation file

genelist = open('zebrafish_anatomy_profiles.txt', 'r')

# defining a multiple dictionary to store functions: genes
genefunc = defaultdict(list)
# defining a multiple dictionary to store gene: functions
genedic = defaultdict(list)

# reading the gene file in implementation
for line in genelist:
    if line != '\n':
        line = line.strip()
        a = line.split('\t') # splitting the line by tab to separate gene name vs functions
        if a[1] != 'gene_name': # excluding the header
            a[1] = a[1].lower() # convert to lower case to avoid mismatches
            #print a[0]
            # here, we are evaluating, so only the proteins with functions will be considered
            if len(a)>1: # ths is required because there are proteins without one single function
                b = a[3].split(',') # seperating the functions by comma
                # store each function in the multi dic by function as the key and genes as values
                for i in b:
                    genedic[a[1]].append(i) # this dictionary stores genes as keys
                    genefunc[i].append(a[1])

print 'total number of genes:', len(genedic)
#defaultdics for genes propagation
count = defaultdict(set)
fcount =defaultdict(set)
for i in genefunc:

    count[i]=set(genefunc[i])
    fcount[i] =set(genefunc[i])

#print count

# fcount =count very unique phenomenon check this out

#print fcount

nlist = nx.nodes(G)

#print nlist

for i in nlist:
    if not i in fcount:
        fcount[i]=set()

#print fcount

#print count
# propagating upwards the number of genes
for i in count:
    alist = list(nx.ancestors(G, i))
    #print alist
    #print i,count[i]
    if  alist:
        for k in alist:
            fcount[k] = fcount[k]|count[i]
            # print fcount
            # print count

#
# for i in count:
#     G.node[i]['label']= i+'('+ str(count[i])+')'
# A = nx.to_agraph(G)
# A.layout('dot')
# A.draw('before.png')

# for i in fcount:
#     if i in count:
#         x= G.node[i]['label']
#         G.node[i]['label']= x+',('+str(fcount[i])+')'
#     else:
#         G.node[i]['label']=i+'(0),('+str(fcount[i])+')'
#
# A = nx.to_agraph(G)
# A.layout('dot')
# A.draw('after.png')
# x= G.node['A']['label']
# print x

##calculating the information content
# calculating the max genecount value
genecountlist=[]
for i in fcount:
    a =len(fcount[i])
    genecountlist.append(a)


mx = max(genecountlist)
print 'maxumuum genecount:',mx

tot= float(mx)
for i in fcount:
    numgenes = len(fcount[i])
    x = (float(numgenes)+1)/(tot+1)
    # print 'hellp',x
    y = math.log10(1/x)
    # print y
    ic[i]=y

#print ic

# for i in ic:
#     G.node[i]['label']= i+'('+ str(fcount[i])+'),('+str(ic[i])+')'
# A = nx.to_agraph(G)
# A.layout('dot')
# A.draw('ic.png')

#writing the information content and the id in a text file
out = open('icdata.txt', 'wb+')
for i in ic:
    out.write('%s %s %s\n'%(i,ic[i],len(fcount[i])))