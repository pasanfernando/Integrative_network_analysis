'''
This code filters the anatomy profile file with a given cutoff for number of functions per gene and number of genes per function
Two output files will be created for the to filtering method
'''
from collections import defaultdict
import re

#obtaining a cutoff value from the user
filterval = int(raw_input('Please enter a cutoff value:'))

#Reading the input gene function relations file

genelist = open('mouse_anatomy_profilesmodified.txt', 'r')

# defining a multiple dictionary to store functions: genes
genefunc = defaultdict(list)
# defining a multiple dictionary to store gene: functions
genedic = defaultdict(list)

# a nested dictionary to store
geneids ={}

# Then convert the above input file into function: gene format

# this indicator is required for leave one out method

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
                geneids[a[1]]['ens_id']= a[2].lower()
                geneids[a[1]]['zfin_id'] = a[0].lower()

                # store each function in the multi dic by function as the key and genes as values

                for i in b:
                    genedic[a[1]].append(i)  # this dictionary stores genes as keys
                    genefunc[i].append(a[1])

print 'Total number of genes in the original profile file:',len(genedic.keys())
print 'Total number of functions in the original profile file:',len(genefunc.keys())

genelist.close()

# output files for gene and function filtered profiles
genefile = open('mouse_anatomy_profilesgenefilter.txt', 'wb+')
funcfile = open('mouse_anatomy_profilesfuncfilter.txt', 'wb+')

genefile.write('zfin_id\tgene_name\tensemble_id\tuberon_annotations\n')
funcfile.write('zfin_id\tgene_name\tensemble_id\tuberon_annotations\n')

# a list to count the number of genes filtered by gene cutoff
genefilters=[]

# a list to count the number of functions filtered by the cutoff
funcfiltered=[]

for i in genedic:

    v = genedic[i]

    if len(v)>= filterval:
        genefile.write('%s\t%s\t%s\t'%(geneids[i]['zfin_id'],i,geneids[i]['ens_id']))
        genefilters.append(i)


        # iterating through the annotations
        for j in v:
            if j.startswith('uberon'):
                if j == v[-1]:
                    genefile.write('%s\n' % (j))
                else:
                    genefile.write('%s,' % (j))




# writing the function filtered file

for i in genedic:
    v = genedic[i]

    temp=[]
    # checking each function whether it has number of annotations higher than the cutoff
    # if yes store them in a temp list
    for j in v:
        if len(genefunc[j])>=filterval:
            funcfiltered.append(j) # appending for counting the number of functions higher than the cutoff
            temp.append(j)

    # write the gene information only if it has functions with higher than the cutoff

    if len(temp)>0:
        funcfile.write('%s\t%s\t%s\t' % (geneids[i]['zfin_id'], i, geneids[i]['ens_id']))

        # iterating through the annotations
        for j in temp:
            if j.startswith('uberon'):
                if j == temp[-1]:
                    funcfile.write('%s\n' % (j))
                else:
                    funcfile.write('%s,' % (j))




print ' Number of genes remained after filtering:',len(set(genefilters))
print ' Number of functions remained after filtering:',len(set(funcfiltered))








