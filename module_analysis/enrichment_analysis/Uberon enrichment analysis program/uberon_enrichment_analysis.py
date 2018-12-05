#!/usr/bin/python
import collections
from scipy import stats
from operator import itemgetter
from collections import OrderedDict

#u = raw_input('Enter anatomy profile file name: ')
#g = raw_input('Enter your gene list: ')
#ou = raw_input('Enter output file: ')

out = open('uberon_enrichmentout.txt', 'wb+')
p = open('newzebrafish_anatomy_profiles.txt', 'r')  # open the output file

# storing in a multiple hash: the gene symbol and uberon id
pheno = collections.defaultdict(list)
uid = []
glist = []
names = {}

# c = 0
# for line in p:
#     x = line.split()
#     y = x[0]
#     z = x[1]
#     n = x[2]
#     if y != 'gene_symbol':
#         c += 1
#         pheno[y].append(z)
#         #names[z] = n
#         if z != 'null':  # removing any null values
#             uid.append(z)  # storing uberon ids in a list

# populating the uberon name dictionary
unames = open('uberonnames.txt','r')

for line in unames:
    if line!='\n' and 'uneron_id' not in line:
        lsplit =line.strip().split('\t')
        names[lsplit[0].lower().replace(':','_')]=lsplit[1].lower()



# reading the gene file in implementation
for line in p:
    if line != '\n':
        line = line.strip()
        a = line.split('\t')  # splitting the line by tab to separate gene name vs functions
        if a[1] != 'gene_name':  # excluding the header
            a[1] = a[1].lower()  # convert to lower case to avoid mismatches
            # print a[0]
            # here, we are evaluating, so only the proteins with functions will be considered
            if len(a) > 1:  # ths is required because there are proteins without one single function
                b = a[3].split(',')  # seperating the functions by comma
                # store each function in the multi dic by function as the key and genes as values
                for i in b:
                    pheno[a[1]].append(i)  # this dictionary stores genes as keys
                    uid.append(i)
                    #genefunc[i].append(a[1])


# print c
print len(uid)
druid = list(set(uid))  # removing duplicates of uberon id list
print len(druid)
# for ln in druid:
#     print ln
# print pheno['abca12']

##removing the duplicates of multi dict by set theory
for key in pheno:
    k = pheno[key]
    pheno[key] = list(set(k))
# print pheno['abca12']

gtot = len(pheno)

gn = open('mygenelist.txt', 'r')  # open the output file

# storing the gene list in a list
for line in gn:
    x = line.split()
    h= x[0].lower()
    glist.append(h)

ltot = len(glist)
##geting statistics
##total number of genes
print 'the total number of genes:', gtot

#
# oddsratio, pvalue = stats.fisher_exact([[20, ltot], [233, gtot]])
# print (pvalue)
# searching for uberon id (fin limb)
# count = 0
# for key in pheno:
#     k = pheno[key]
#     if 'UBERON:0001003' in k:
#         out.write('%s\n ' % (key))
#         count += 1
#
# out.write(str(count))
dummyid = ['UBERON:0007380']
final = {}
finalsorted = {}

for item in druid:
    gcount = 0
    ggene =[]
    for key in pheno:
        k = pheno[key]
        if item in k:
            gcount += 1
            ggene.append(key)


    intl = list(set(ggene) & set(glist))
    cintl = len(intl)
    # if item == 'uberon_0007380':
    #     print 'gcount',gcount
    #     print ggene
    #     print intl
    #print gcount
    ngcount = gtot - gcount
    nlcount = ltot - cintl
    oddsratio, pvalue = stats.fisher_exact([[cintl, nlcount], [gcount, ngcount]])
    #if want to check the fisher function use the line below
    #print item,pvalue,cintl, nlcount,gcount, ngcount
    final[item] = pvalue

#finalsorted = sorted(final.items(), key=itemgetter(1))
finalsorted = OrderedDict(sorted(final.items(), key=lambda t: t[1]))

print  final
print finalsorted
out.write('Uberon_id\tP-value\tTerm name\n')

for k in finalsorted:
    out.write('%s\t%s\t%s\n'%(k,finalsorted[k],names[k]))
    #print k,finalsorted[k]
