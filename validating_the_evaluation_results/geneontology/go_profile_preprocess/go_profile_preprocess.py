'''
When pre-processing, use the reconcilled anatomy profiles file. Because it already contains the genes that are common to both anatomy and the string network

Inputs: reconciled anatomy profiles file, raw GO profiles file
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


genelist = open('zebra_goprofiles.txt', 'r')

# reading the GO profiles

for line in genelist:
    if line != '\n':
        line = line.strip()
        a = line.split('\t') # splitting the line by tab to separate gene name vs functions
        if a[1] != 'gene_name': # excluding the header
            a[1] = a[1].lower() # convert to lower case to avoid mismatches
            geneids[a[1].lower()] = a[0].lower()
            #print a[0]
            # here, we are evaluating, so only the proteins with functions will be considered
            if len(a)>1: # ths is required because there are proteins without one single function
                b = a[3].split(',') # seperating the functions by comma
                # store each function in the multi dic by function as the key and genes as values
                for i in b:
                    genedic[a[1]].append(i) # this dictionary stores genes as keys

genelist.close()

# Opening the reconciled anatomy ontology profile
anprofile = open('zebrafish_anatomy_profilesmodified.txt', 'r')

anatomygeneset =[]

for line in anprofile:
    if line != '\n':
        line = line.strip()
        a = line.split('\t') # splitting the line by tab to separate gene name vs functions
        if a[1] != 'gene_name': # excluding the header
            a[1] = a[1].lower() # convert to lower case to avoid mismatches
            anatomygeneset.append(a[1])
anprofile.close()
print 'number of genes in the reconciled anatomy profile:',len(anatomygeneset)
# writing the gene ontology profile


genefile = open('zebra_goprofiles_modified.txt', 'wb+')

genefile.write('zfin_id\tgene_name\tensemble_id\tgo_annotations\n')



# a list for matched genes
matchgenes =[]

for i in genedic:
    # getting the functions
    if i in anatomygeneset:
        matchgenes.append(i)

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



print 'Number of genes in the pre-processed file:',len(matchgenes)
print 'Number of functions in the new file:', len(genefunc)
print 'Number of mismatched genes:', len(set(anatomygeneset)-set(matchgenes))
print set(anatomygeneset)-set(matchgenes)

