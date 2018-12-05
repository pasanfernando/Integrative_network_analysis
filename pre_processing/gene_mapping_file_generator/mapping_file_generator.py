'''
Generate a mapping file in the following format for a given species
gene name   string id   ensemble id
'''
from collections import defaultdict
import re


# Opening the stringID to ensemble id mapping file
enstostrmapfile = open('mouse_ensemblemappingstostring.txt','r')

# defining a dic to store: ensembleid: stringID mapping
ensemble_to_str ={}

for line in enstostrmapfile:
    if line != '\n':
        b = line.split()
        if 'string_id' not in line:
            # print line
            # storing string id and ensemble id in the  mapping dictionary

            ensemble_to_str[b[1].lower()] = b[0].lower()

print len(ensemble_to_str)

# closing the stringID to ensemble id mapping file
enstostrmapfile.close()

# Opening the gene name to ensemble id mapping file from biomart
biomart_mapfile = open('mouse_mart_export.txt','r')

# defining a dic to store: gene name: ensembleid mapping
gene_to_ensemble ={}

# defining a dic to store:  ensembleid: gene name  mapping
ensemble_to_gene ={}
for line in biomart_mapfile:
    if line != '\n':
        b = line.split('\t')
        #print b
        if 'Gene stable ID' not in line:
            # print line
            # storing gene name and ensemble id in the  mapping dictionary

            gene_to_ensemble[b[2].lower().strip()] = b[0].lower()
            ensemble_to_gene[b[0].lower()] = b[2].lower().strip()

print len(gene_to_ensemble)
#print gene_to_ensemble

# finding duplicate ensemble ids, which has different gene names
ens_duplicates = defaultdict(list)



for k,v in gene_to_ensemble.items():
    ens_duplicates[v].append(k)

# single line statement to check the duplicates; not sure about the accuracy
#print [v for k,v in ens_duplicates.items() if len(v)>1]

# following loop will print gene names that have the same ensemble id
print 'if there are different genes with the same ensemble id, they will be printed below.'
for k,v in ens_duplicates.items():
    if len(v)>1:
        print v

# generating the mapping file
mapfile = open('three_way_mapping.txt', 'wb+')

mapfile.write('gene_name\tstring_ID\tensemble_id\n')

matched_strings =[]
for k,v in gene_to_ensemble.items():
    # checking for the string id
    if v in ensemble_to_str:
        strid = ensemble_to_str[v]
        matched_strings.append(strid)
    else:
        strid = 'na'

    mapfile.write('%s\t%s\t%s\n'%(k,v,strid))

mapfile.close()


# printing the mapping stats
print 'number of gene names and ensemble ids:', len(gene_to_ensemble)
print 'number of matched string ids:', len(matched_strings)
print 'number of matched string ids after removing duplications:', len(set(matched_strings))

########################################################################################################### # check whether all phenoscape gene names are matching with the biomart gene names
#
# # adapt this code for the anatomy network gene list at the begining
#
# # defining a dictionary to store phenoscape gene name to ensemble id mapping
# phenoscape_ensmble ={}
#
# # opening the phenoscape mapping file that prashanti sent
# in1 = open('KBGenes_2016_EnsemblIDs.tsv', 'r')
# total =0
#
# for line in in1:
#     if line != '\n':
#         if '?gene' not in line:
#             a = line.split('\t')
#             gid = a[0]
#             gn = a[1]
#             org = a[2]
#             ensid = a[3].strip().lower()
#             total = total +1
#
#             # extracting gene id
#             gid = gid.strip('<|>')
#             #print gid
#             #print gn
#             # extracting gene name
#             result = re.search('"(.*)"', gn)
#             #print result
#             x = result.group(1).strip().lower()
#             # x = x.strip()
#             # x = x.lower()
#
#             # extracting species name
#             result = re.search('"(.*)"', org)
#             # print result
#             y = result.group(1)
#             y = y.strip()
#
#             # change the species name to match the species you are interested in
#             if y == 'Mus musculus':
#                 phenoscape_ensmble[x]=ensid
#
#
# #print phenoscape_ensmble
# print 'the length of the phenoscape gene list is:',len(phenoscape_ensmble)
#
# # checking the phenoscape genes in the biomart list
#
# for k,v in phenoscape_ensmble.items():
#     if k not in gene_to_ensemble:
#         if v in ensemble_to_gene:
#             updated_name = ensemble_to_gene[v]
#             print k,v,updated_name
#
#         else:
#             print k,v,'no ensemble match'

