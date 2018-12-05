
'''

Extracting the genes associated with the parts of a given anatomical entity
 e.g,: Parts of the pectoral fin


This should include 'is a' and 'part of' relationships

Outputs the details of the newly extracted genes

Also outputs a gene list that includes both new genes and original genes associated with the original term
    Use this for the module detection procedure

There is also a user prompt for validating the predicted genes in the module
If this option is selected the predicted genes will be checked against the new genes that are annotated to subterms
and will be printed in the gene details file

Please change the term of interest in the code for any uberon function you like

'''

import networkx as nx
import itertools as it
import re
from collections import defaultdict


# Load the Uberon data into a graph data structure including the part of relationships

G = nx.DiGraph()
uberon = {}
name = {}


#reading the uberon file
p = open('uberon.obo', 'r')

for line in p:
    if line.startswith('[Typedef]'):
        break

    if line.startswith('id:'):
        x = line.split()
        y = x[1].lower().replace(":","_")

        if G.has_node(y) == False:
            G.add_node(y)

    if line.startswith('name:'):
        k = re.sub('name: ','',line,1)
        z = k.strip()
        d = z.replace(" ","_").lower()
        uberon[y] = d
        name[z] = y

    if line.startswith('is_a: '):
        z = line.split()
        k = z[1].lower().replace(":","_")

        # print k
        G.add_edge(k, y)

    # elif line.startswith('relationship: ') or line.startswith('intersection_of:'):
    # if line.startswith('relationship: ') or line.startswith('intersection_of:'):
    # if line.startswith('relationship: part_of'):
    #     if 'UBERON:' in line:
    #         line = line.strip()
    #         m = 'UBERON:'+re.search('UBERON:(.+?) ', line).group(1)
    #         G.add_edge(m.lower().replace(":","_"), y, color='red')

    if line.startswith('relationship: part_of'):
        if 'UBERON:' in line:
            line = line.strip()
            m = 'UBERON:'+re.search('UBERON:(.+?) ', line).group(1)
            G.add_edge(m.lower().replace(":","_"), y, color='red')
    #
    if line.startswith('relationship: has_potential_to_develop_into'):
        if 'UBERON:' in line:
            line = line.strip()
            m = 'UBERON:'+re.search('UBERON:(.+?) ', line).group(1)
            G.add_edge(m.lower().replace(":","_"), y, color='red')

inp = raw_input('Do you want to perform a validation of predicted genes (y or n):')

# retrieving the Uberon term of interest (use the colon between uberon and the number

function = 'uberon_0000152'
# please enter the term search name
termname='pelvic'


# getting the decsendents of the the given function

child = nx.descendants(G, function)

# for i in child:
#     ch1=nx.descendants(G, i)
#     child=ch1 | child

print child
# for i in child:
#     print uberon[i]


#Opening anatomical profile file
genelist = open('newzebrafish_anatomy_profiles.txt', 'r')

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

# a complex dictionary to store genes annotated to subparts
childgenes = defaultdict(list)

# a complex dictionary to store genes annotated to search term genes
termgenes = defaultdict(list)

# termlist = genedic['tbx4']
# for i in termlist:
#     print uberon[i]

# finding terms with the given searchterm
for i in genefunc:
    if termname in uberon[i]:
        if i not in child:
            print uberon[i]
            l2=genefunc[i]
            for j in l2:
                termgenes[j].append(uberon[i])





# retrieving the genes annotated to each subpart of the original term

for i in child:
    tname = uberon[i]
    l1 = genefunc[i]
    for j in l1:
        childgenes[j].append(tname)

#print childgenes

# checking how many of these genes are asscoiated with the original term

# original genes
origenes = genefunc[function]
print set(origenes)
extragenes = set(childgenes.keys()) - set(origenes)
olgenes = set(childgenes.keys()) & set(origenes)



totgenes = set(childgenes.keys()) | set(origenes)

# for search term genes
termextragenes = set(termgenes.keys())-totgenes
termcommongenes =set(termgenes.keys())&totgenes
#print extragenes
finaltotgenes = totgenes | set(termgenes.keys())


# writing the details about the genes associated with the parts
gstats = open('gene_details.txt', 'wb+')

gstats.write('New genes that are not annotated to original term: %s\n'%(len(extragenes)))
gstats.write('gene_name\taccociated_terms\n')
for i in extragenes:
    gstats.write('%s\t'%(i))
    terms = childgenes[i]
    for j in terms:
        if j == terms[-1]:
            gstats.write('%s\n' % (j))
        else:
            gstats.write('%s,' % (j))

gstats.write('\n')
gstats.write('\n')
gstats.write('Genes that are both annotated to the original term and the subterms: %s\n'%(len(olgenes)))
gstats.write('gene_name\taccociated_terms\n')
for i in olgenes:
    gstats.write('%s\t'%(i))
    terms = childgenes[i]
    for j in terms:
        if j == terms[-1]:
            gstats.write('%s\n' % (j))
        else:
            gstats.write('%s,' % (j))

gstats.write('\n')
gstats.write('\n')

gstats.write('Number of total genes including extra genes and the overlapped genes (only subparts): %s\n'%(len(totgenes)))

gstats.write('\n')
gstats.write('\n')
gstats.write('New genes that are only associated to terms with search phrase: %s\n'%(len(termextragenes)))
gstats.write('gene_name\taccociated_terms\n')
for i in termextragenes:
    gstats.write('%s\t'%(i))
    terms = termgenes[i]
    for j in terms:
        if j == terms[-1]:
            gstats.write('%s\n' % (j))
        else:
            gstats.write('%s,' % (j))

gstats.write('\n')
gstats.write('\n')
gstats.write('Genes that are both annotated to the original term and the subterms search terms: %s\n'%(len(termcommongenes)))
gstats.write('gene_name\taccociated_terms\n')
for i in termcommongenes:
    gstats.write('%s\t'%(i))
    terms = termgenes[i]
    for j in terms:
        if j == terms[-1]:
            gstats.write('%s\n' % (j))
        else:
            gstats.write('%s,' % (j))

gstats.write('Final number of genes: %s\n' % (len(finaltotgenes)))
# writing the total genes as a list
glist = open('genelist_including_parts.txt', 'wb+')
for i in finaltotgenes:
    glist.write('%s\n' % (i))

glist.close()

# if the prediction validation option is selected, perform the prediction validation
if inp== 'y':
    # opening the predicted gene list file
    genelist = open('only_predicted_genes.txt', 'r')

    # a list to store predicted genes
    predgenes =[]

    for line in genelist:
        if line != '\n':
            line = line.strip()
            predgenes.append(line)

    # validating the predicted genes by getting the intersection
    confgenes = set(predgenes)& extragenes

    # writing the confirmed genes in the genes details file
    gstats.write('Number of predicted genes that are validated or overlapped with new genes: %s\n' % (len(confgenes)))
    gstats.write('gene_name\taccociated_terms\n')
    for i in confgenes:
        gstats.write('%s\t' % (i))
        terms = childgenes[i]
        for j in terms:
            if j == terms[-1]:
                gstats.write('%s\n' % (j))
            else:
                gstats.write('%s,' % (j))

gstats.close()
