'''
This code randomly removes a user given number of functions from the original profile
Then also creates a new profile from the removed functions
Change the first two variables accordingly
'''
from collections import defaultdict
import re
import numpy as np

# selecting only the functions that have more genes than the threshold
functhreshold =10

# the number of functions that needs to be filtered
numfunc =30

#Reading the input gene function relations file

genelist = open('zebrafish_anatomy_profiles.txt', 'r')

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



# a list for the function pool
func_pool =[]

for i in genefunc:
    if len(genefunc[i])>=functhreshold:
        func_pool.append(i)

# a new multiple dictionary to store gene: functions for the randomly removed functions
newgenedic = defaultdict(list)


# output files for function filtered profiles

funcfile = open('zebra_randomfuncfiltered_profiles.txt', 'wb+')


funcfile.write('zfin_id\tgene_name\tensemble_id\tuberon_annotations\n')

# a list to count the number of functions filtered by the cutoff
funcfiltered=[]


# randomly selecting a defined number of functions
randarr=np.random.choice((len(func_pool)-1),numfunc,replace=False)

# a list to store randomly selected functions
randfuncs=[]

for i in randarr:
    randfuncs.append(func_pool[i])

# writing the function filtered file

for i in genedic:
    v = genedic[i]

    temp=[]
    # checking each function whether it has number of annotations higher than the cutoff
    # if yes store them in a temp list
    for j in v:
        if j in randfuncs:
            newgenedic[i].append(j)
        else:
            temp.append(j)

    # write the gene information only if it has functions with higher than the cutoff
    # removing duplicates
    temp = list(set(temp))
    if len(temp)>0:
        funcfile.write('%s\t%s\t%s\t' % (geneids[i]['zfin_id'], i, geneids[i]['ens_id']))

        # iterating through the annotations
        for j in temp:
            if j.startswith('uberon'):
                if j == temp[-1]:
                    funcfile.write('%s\n' % (j))
                else:
                    funcfile.write('%s,' % (j))

print 'the randomly selected functions are given below'
print randfuncs
print len(randfuncs)
# writing a new profile only for the functions that was removed randomly; can be used for multifunction evaluation

genefile = open('zebra_onlyrandomlypicked_profile.txt', 'wb+')

genefile.write('zfin_id\tgene_name\tensemble_id\tuberon_annotations\n')


for i in newgenedic:
    # getting the functions
    v = newgenedic[i]
    # removing duplicates
    v = list(set(v))
    genefile.write('%s\t%s\t%s\t'%(geneids[i]['zfin_id'],i,geneids[i]['ens_id']))

    # iterating through the annotations
    for j in v:
        if j.startswith('uberon'):
            if j == v[-1]:
                genefile.write('%s\n' % (j))
            else:
                genefile.write('%s,' % (j))

print 'Number of genes in the only randomly picked function file:', len(newgenedic)






