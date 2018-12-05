'''
Generate gene ontology profiles for zebrafish
'''
from collections import defaultdict
import re


# defining a multiple dictionary to store gene: functions
genedic = defaultdict(list)

# defining a multiple dictionary to store functions: genes
genefunc = defaultdict(list)

# a dictionary to store zfin ids
geneids ={}

# read the gene ontology profiles


genelist = open('gene_association.zfin', 'r')

# reading the gene file in implementation
# do not contain a header line so don't worry about excluding the header
for line in genelist:
    if line != '\n':
        line = line.strip()
        a = line.split('\t') # splitting the line by tab

        # only considering the biological process
        if a[8]== 'P':
            go= a[4].lower().replace(':','_').strip('\n')
            genedic[a[2].lower()].append(go)
            geneids[a[2].lower()]=a[1].lower()


print genedic
print 'Number of genes in the profile:',len(genedic)

# writing the gene ontology profile


genefile = open('zebra_goprofiles.txt', 'wb+')

genefile.write('zfin_id\tgene_name\tensemble_id\tgo_annotations\n')


for i in genedic:
    # getting the functions
    v = genedic[i]
    # removing duplicates
    v = list(set(v))
    # filling ensemble column with dummy
    genefile.write('%s\t%s\tdummy\t'%(geneids[i],i))

    # iterating through the annotations
    for j in v:
        genefunc[j].append(i)
        if j.startswith('go'):
            if j == v[-1]:
                genefile.write('%s\n' % (j))
            else:
                genefile.write('%s,' % (j))


print 'Number of functions in the new file:', len(genefunc)

