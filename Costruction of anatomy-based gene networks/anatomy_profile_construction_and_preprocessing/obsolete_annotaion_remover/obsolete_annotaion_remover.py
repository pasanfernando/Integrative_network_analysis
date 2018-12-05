__author__ = 'compaq'

import networkx as nx
import operator
import re
from collections import defaultdict

'''
Remove the annotations that are attached to obsolete Uberon terms
'''

####### reading Uberon #########################


ub = nx.DiGraph()
uberon = {}
name = {}
obsolete =[]

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

    if line.startswith('is_obsolete: '):
        z = line.split()
        # extracting the obsolete terms
        if z[1]== 'true':
            obsolete.append(y)
        # print k


    # elif line.startswith('relationship: ') or line.startswith('intersection_of:'):
    #     if 'UBERON:' in line:
    #         line = line.strip()
    #         m = 'UBERON:'+re.search('UBERON:(.+?) ', line).group(1)
    #         ub.add_edge(m, y, color = 'red')

# removing the obsolete annotations from the profiles
# opening the annotation file
genelist = open('mouse_anatomy_profilesduplicate_removed.txt', 'r')

# defining a multiple dictionary to store functions: genes
genefunc = defaultdict(list)
# defining a multiple dictionary to store gene: functions
genedic = defaultdict(list)

# a nested dictionary to store
geneids ={}

# Then convert the above input file into function: gene format

# this indicator is required for leave one out method

#obsolete genes annotations
obsannotations=[]

# reading the gene file
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
                geneids[a[1]]['ens_id'] = a[2].lower()
                geneids[a[1]]['zfin_id'] = a[0].lower()
                # store each function in the multi dic by function as the key and genes as values
                for i in b:
                    if i not in obsolete:
                        genedic[a[1]].append(i)  # this dictionary stores genes as keys

                    else:
                        obsannotations.append(i)

#writing the new profiles
genefile = open('mouse_finalfull_profiles.txt', 'wb+')
genefile.write('zfin_id\tgene_name\tensemble_id\tuberon_annotations\n')

for i in genedic:

    v = genedic[i]


    genefile.write('%s\t%s\t%s\t'%(geneids[i]['zfin_id'],i,geneids[i]['ens_id']))



    # iterating through the annotations
    for j in v:
        if j.startswith('uberon'):
            if j == v[-1]:
                genefile.write('%s\n' % (j))
            else:
                genefile.write('%s,' % (j))

print 'number of obsolete annotations:',len(obsannotations)
print obsannotations