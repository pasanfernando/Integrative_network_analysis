__author__ = 'compaq'

import networkx as nx
import operator
import re

'''
Before running the code, make sure the lowest IC value is not less than -1; if it is less than -1, it will be not recognized as a common ancestor
'''

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


# You only can enter capital letters: change this please
# x = raw_input('enter first node:')
# y= raw_input('enter second node:')

# below is the original read dot method, which doesnt work in mac
#G = nx.read_dot('name.dot')



# this method is specific for IC based algorithms; it picks the most reason common ancestor based on highest IC if there are 2 candidates
# If you just want to pick the common ancestor based on shortest path, modify this method
def commonancestor(x,y):
    #take two nodes get ancestors of both of them
    a = nx.ancestors(G, x)
    #print a
    # adding itself to the list of ancestors
    a.add(x)
    b = nx.ancestors(G, y)
    #print b
    b.add(y)
    # check wether one is an ancestor of another
    if x in b :
        return x
    elif y in a:
        return y
    #convert the ancestors to sets
    else:
        dis = {}
        a = set(a)
        b = set(b)

    #get the intersection into another list
        c = a & b
        c =list(c)

    #pick one node from the initial list, if the intersection is more than one get the shortest paths to all the items and store it in a hash
        if len(c) == 1:
            return c[0]
        else:
            #print 'yayayayayay'
            #print c
            for i in c:
                l = nx.shortest_path_length(G,i,x)
                dis[i]= l
            #get the key with the shortest distance (_path) which is the common ancestor
            #return min(dis.iteritems(), key=operator.itemgetter(1))[0]
            # there can be more than one ancestor with the shortest path; pick the one with highest IC

            # make sure all your ic values are greater than -1, if they are not they won't be selected
            maxic = -1
            term_maxic = None
            #print dis
            for key in dis:
                if dis[key] == min(dis.values()):
                    if ic[key]> maxic:
                        #print key
                        maxic= ic[key]
                        term_maxic =key

            #print term_maxic
            return term_maxic






infile = open('icdata.txt', 'r')
ic ={}
for line in infile:
    if line != '\n':
        a,b,c= line.split()
        b= float(b)
        ic[a]=b

def icsimilarity(x,y):
    z = commonancestor(x,y)
    return ic[z]

def icmethods(x,y,type):

    print 'common ancestor:',commonancestor(x,y)
    ics = icsimilarity(x,y)
    print ics

    if type=='res':
        return ics

    else:

        lin = (2*ics)/(ic[x]+ic[y])
        print 'similarity based on lin method',lin

        if type=='lin':
            return lin

        elif type=='slick':
            scl = lin* (1+ics)
            print 'similarity based onschlicker method',scl
            return scl

        else:
            print 'not a valid method; please pass the correct parameter!!'
