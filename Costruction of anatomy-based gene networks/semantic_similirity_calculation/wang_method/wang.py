
import networkx as nx
import collections
from collections import defaultdict
import re

# I commented all the print statements and drawing statements for subgraphs; remove the comments if you check the functionality of the codes

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
################################################
#print name['parietal bone']
##########################################################################################################

def wang(x):
    # input a specific node for the analysis

    #get all the ancestors for the given node
    name = x + 'sub.dot'
    a = nx.ancestors(G, x)
    a = list(a)
    a.append(x)
    sb = nx.subgraph(G, a)
    #nx.drawing.nx_agraph.write_dot(sb,name)

    #start with the lowest node ( input node)
    # create a new dictionary for sv(a) and store the value for input node (1)
    sv = {}
    l2 = []
    sv[x]= 1
    l2.append(x)
    d1=  wang1(l2,sb,sv)
    # print d1
    k = sum(d1.itervalues())
    # print k
    return d1,k

    #wang1(x,sb,sv)


#get the immediate ancestors calculate their sv(a) by considering the sv(a) values of their children
def wang1(y,z,d):
    pt = []
    for j in y:
        kt = z.predecessors(j)
        pt = pt + kt
    a3 = []
    #print pt
    for m in pt:
        a2 = nx.ancestors(z, m)
        a1 = list(a2)
        a3= a3+a1
    #print a3
    a3 = list(set(a3))

    if a3:
        pt = list(set(pt)-set(a3))

        # for m1 in a3:
        #     if m1 in pt:
        #         pt.remove(m1)


    if not pt:
        # print 'list is empty'

        # print d
        return d

    else:
        #print pt
        for i in pt:
            ch =z.successors(i)
            # print ch
            li = []
            for j in ch:
                v = d[j]*0.8
                li.append(v)
            # print li
            m = max(li)
            d[i]= m
        # for i in pt:
    return wang1(pt,z,d)


    # print d
    # print sum(d.itervalues())
# wang1('I')
# s,m=wang('D')
# print s
# print m

######### calculating the semantic similarity between two terms a, b

# a = 'G'
# b = 'F'

#finding the intersection of ancestors for the two terms
def wangsim(a,b):
    a1 = nx.ancestors(G, a)
    a2 = nx.ancestors(G, b)
    #print a1
    #print a2
    a1 = list(a1)
    a1.append(a)
    a1 = set(a1)
    a2 = list(a2)
    a2.append(b)
    a2 = set(a2)
    a3 = a1 & a2
    #print a3
    a3l =list(a3)

    # calculating the S(A) values for the two terms
    ka,sa=wang(a)
    kb,sb=wang(b)

    #print sa,sb
    #print ka,kb
    sigma =0
    for i in a3l:
        va = ka[i]
        vb = kb[i]
        sigma = sigma + va + vb

    #print sigma

    sim = sigma/(sa+sb)
    #print sim
    return sim