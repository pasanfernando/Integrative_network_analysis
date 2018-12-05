import networkx as nx
from collections import defaultdict
from operator import itemgetter
from collections import OrderedDict
import numpy as np
import re


#import matplotlib.pyplot as plt

""" Randomizes the profiles based on genes to compare the performance with the actual profile
    The number of functions annotated to a gene won't change, only number of genes per function might change from the original
    Single function method requires function based randomization
    
    Update
    when randomly selecting functions from the function pool, use the entire Uberon ontology, instead of getting only the functions in the profile
    """

############################################################################# ###########


#Reading the input gene function annotation file

genelist = open('zebrafish_anatomy_profiles.txt', 'r')

# defining a multiple dictionary to store functions: genes
genefunc = defaultdict(list)
# defining a multiple dictionary to store gene: functions
genedic = defaultdict(list)

# a nested dictionary to store zfin ids and ensemble ids
geneids ={}
# Then convert the above input file into function: gene format

#######################################################################
# reading the Uberon ontology
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

funcpool=[]
# Generating the function pool
for i in ub.nodes():
    if i.startswith('uberon'):
        funcpool.append(i)

print 'number of uberon terms in the function pool:', len(funcpool)

# this indicator is required for leave one out method



# reading the gene file in implementation
for line in genelist:
    if line != '\n':
        line = line.strip()
        a = line.split('\t') # splitting the line by tab to separate gene name vs functions
        if a[1] != 'gene_name': # excluding the header
            a[1] = a[1].lower()  # convert to lower case to avoid mismatches
            # print a[0]
            # here, we are evaluating, so only the proteins with functions will be considered
            if len(a) > 1:  # ths is required because there are proteins without one single function
                b = a[3].split(',')  # seperating the functions by comma
                geneids[a[1]] = {}
                geneids[a[1]]['ens_id']= a[2].lower()
                geneids[a[1]]['zfin_id'] = a[0].lower()

                # store each function in the multi dic by function as the key and genes as values

                for i in b:
                    genedic[a[1]].append(i)  # this dictionary stores genes as keys
                    genefunc[i].append(a[1])

# getting all the functions in the original profile
#funcpool = genefunc.keys()
print 'number of genes in the original file',len(genedic)
print 'number of functions in the original file',len(genefunc)
#print funcpool[5404]


# new dictionaries to store randomized annotations
# defining a multiple dictionary to store functions: genes
newgenefunc = defaultdict(list)
# defining a multiple dictionary to store gene: functions
newgenedic = defaultdict(list)
# going trough different functions
for i in genedic:
    genel= len(genedic[i])
    #print genedic[i]
    # generating random numbers
    randarr=np.random.choice((len(funcpool)-1),genel,replace=False)
    # appending the genes randomly
    for j in randarr:
        # populating the genefunc dictionary
        newgenedic[i].append(funcpool[j])
        #populating the gene dictionary
        newgenefunc[funcpool[j]].append(i)

    #print newgenedic[i]


# writing the randomized by function profile
genefile = open('zebra_gene_randomizedprofiles.txt', 'wb+')

genefile.write('zfin_id\tgene_name\tensemble_id\tuberon_annotations\n')


for i in newgenedic:
    # getting the functions
    v = newgenedic[i]

    genefile.write('%s\t%s\t%s\t'%(geneids[i]['zfin_id'],i,geneids[i]['ens_id']))

    # iterating through the annotations
    for j in v:
        if j.startswith('uberon'):
            if j == v[-1]:
                genefile.write('%s\n' % (j))
            else:
                genefile.write('%s,' % (j))

print 'Number of genes in the new file:', len(newgenedic)
print 'Number of functions in the new file:', len(newgenefunc)