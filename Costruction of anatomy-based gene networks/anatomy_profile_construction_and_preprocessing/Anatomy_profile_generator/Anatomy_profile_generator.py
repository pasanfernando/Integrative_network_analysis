'''
Generate a anatomy profiles for mouse and zebrafish using the scala output
'''
from collections import defaultdict
import re


# method to output annotation statistics for the two organisms
def annotation_stats(orgdic, spname):

    #counters to count uberon, gene ontology and other annotations
    ubtotal =0
    gototal =0
    othertotal =0

    for k,v in orgdic.items():
        # we need to only keep the uberon ids otherwise a comma will be printed at the end of some uberon ids (last element)
        # this tempory list only copies the uberon ids to the orgdic
        templist =[]
        for i in v:
            if i.startswith('uberon'):
                ubtotal+=1
                templist.append(i)

            elif i.startswith('go'):
                gototal+=1

            else:
                othertotal+=1
        orgdic[k]= templist

    statout.write('********************************************************************************\n')
    statout.write('Annotation statistics for %s is given below\n'%(spname))
    statout.write('Total number of gene ids: %s\n'%(len(orgdic)))
    statout.write('Total number of Uberon annotations: %s\n'%(ubtotal))
    statout.write('Average Uberon annotations per gene: %f\n'%(float(ubtotal)/len(orgdic)))
    statout.write('Total number of GO annotations: %s\n' % (gototal))
    statout.write('Total number of other annotations: %s\n' % (othertotal))
    statout.write('\n')
    statout.write('\n')
    return

# Given a list, writes the components of the list in the statistical output file; e.g: genes without gene names and ensemble ids
def listwrite(list):
    for i in list:
        statout.write('%s\n'%(i))

    return

###############################################################################################################################
# Opening the scala output

scalafile = open('jim_new_output.txt','r')

# defining a multidic to store zfinid: uberon annotations
zfin_uberon = defaultdict(list)

# defining a multidic to store mgiid: uberon annotations
mgi_uberon = defaultdict(list)
# important: make sure all the names and ids are converted to lowercase before storing them

for line in scalafile:
    if line != '\n':
        a = line.split('\t')
        #print a

        # extracting tha anatomy annotations first because they are common to both zebrafish and mouse
        if 'http://purl.obolibrary.org/' in a[1]:
            b1 = a[1].strip().lower().split('/')
            #print b1[4]

        # pre-processing the urls to extract zfin or mgi id
        # distinguishing between mouse and zebrafish
        # first extract zebra ids
        if 'http://zfin.org' in a[0]:
            # extracting the zfin id at the end of the string url
            # i am using split method, if you want to use regex u can do that
            #print a
            b = a[0].lower().split('/')
            #print b[3]
            if b[3].startswith('zdb'):
                zfin_uberon[b[3]].append(b1[4])
            # printing the ids that do not start with zdb
            else:
                print b[3]

        elif 'http://www.informatics.jax.org' in a[0]:
            # extracting the mgi id at the end of the string url
            m1 = a[0].lower().split('/')
            if m1[4].startswith('mgi'):
                mgi_uberon[m1[4]].append(b1[4])
            # printing the ids that do not start with mgi
            else:
                print m1[4]

# print zfin_uberon
# print mgi_uberon


scalafile.close()

# removing duplicate annotations in zfin dic
# for these steps u can write a general method but there are only two organisms and code snippet is too small
for k,v in zfin_uberon.items():
    tempv= set(v)
    zfin_uberon[k]=list(tempv)


# removing duplicate annotations in mouse dic
for k,v in mgi_uberon.items():
    tempv= set(v)
    mgi_uberon[k]=list(tempv)

# opening an output file to save the statistics
statout = open('annotation_stats.txt', 'wb+')
# calling the annotation stats method

annotation_stats(zfin_uberon,'zebrafish')
annotation_stats(mgi_uberon,'mouse')



# generating the annotation files

# reading biomart downloads and storing them in nested dictionaries

# two nested dictionaries to store mapping data
#zebrafish dic: key: zfin id; values: ensemble and gene ids
zebramap ={}
mousemap ={}

# opening the zebrafish biomart file
zebmart = open('zebra_mart_export_9_7_17.txt', 'r')

for line in zebmart:
    if line != '\n' and 'Gene name'not in line:
        a = line.strip().split('\t')
        #print a
        #print a[6]
        # checking wether list index 6, which corresponds to zfin id exist; some records do not have this index
        # it seems like all ensemble ids and all gene names are present for zfin id
        if len(a)>6:
            # filling the nested dictionary
            zebramap[a[6].lower()]={}
            zebramap[a[6].lower()]['ensemble_id'] = a[0].lower()
            zebramap[a[6].lower()]['gene_name'] = a[2].lower()
        # it seems like all ensemble ids and all gene names are present for zfin id
        # if a[2] == '':
        #     print a[0]

print 'Number of zfin ids in the biomart file:',len(zebramap)

# opening the mouse biomart file
mousemart = open('mouse_mart_export_9_7_17.txt', 'r')

for line in mousemart:
    if line != '\n' and 'Gene name'not in line:
        a = line.strip().split('\t')
        #print a
        #print a[6]
        # checking wether list index 4, which corresponds to mgi id exist; some records do not have this index
        # it seems like all ensemble ids and all gene names are present for mgi id
        if len(a)>4:
            # filling the nested dictionary
            mousemap[a[4].lower()]={}
            mousemap[a[4].lower()]['ensemble_id'] = a[0].lower()
            mousemap[a[4].lower()]['gene_name'] = a[6].lower()

mousemart.close()
print 'Number of mgi ids in the biomart file:',len(mousemap)

# Since some zfin ids are not mapping to the biomart download, i am using the zfin download file to match the remainder
pe = open('phenoGeneCleanData_fish_2017.09.11.txt', 'r')

# dictionary to store zfin id to gene symbol mappings
zfinid_to_name ={}

for line in pe:
    if line != '\n' and  'Gene Symbol' not in line:
        #print line
        a = line.strip().split('\t')
        zfinid_to_name[a[2].lower()]= a[1].lower()

pe.close()
# generating the profiles for zebrafish

#print zfinid_to_name

############################################################################### code for zfin profile generation ##############################
# opening the output file to write zebrafish profiles
zebraprof = open('zebrafish_anatomy_profiles.txt','wb+')

# writing the header line
zebraprof.write('zfin_id\tgene_name\tensemble_id\tuberon_annotations\n')

# a list to store zfinids without uberon annotations (they have GO annotations)]
ids_without_uberon =[]

# a list to store zfin ids that only miss ensemble ids (these have gene names)
zebra_ensemble_missing =[]

# a list to store zfin ids that  miss both ensemble ids and gene names
zebramissing =[]


# iterating through  the zfin ids from the scala output
for k,v in zfin_uberon.items():
    # this if statement is needed to remove the genes without uberon ids (the ones with only GO annotations)
    if v:
        # these are the genes that have both gene name and ensemble ids
        if k in zebramap:
            info = zebramap[k]
            zebraprof.write('%s\t%s\t%s\t'%(k,info['gene_name'],info['ensemble_id']))

        # These genes only have a gene name
        elif k in zfinid_to_name:
            zebra_ensemble_missing.append(k)
            zebraprof.write('%s\t%s\tnot_found\t' % (k,zfinid_to_name[k]))

        # These genes do not have gene names and ensemble ids
        else:
            #print k
            zebraprof.write('%s\tnot_found\tnot_found\t' % (k))
            zebramissing.append(k)

        # iterating through the annotations
        for i in v:
            if i.startswith('uberon'):
                if i == v[-1]:
                    zebraprof.write('%s' % (i))
                else:
                    zebraprof.write('%s,' % (i))

        zebraprof.write('\n')
    else:

        ids_without_uberon.append(k)

zebraprof.close()

print len(zebramissing)

# writing the list of zfin ids that are missing some data in statistics output: gene names, ensemble ids

# the ones that only miss ensemble ids
statout.write('zfin ids that only miss ensemble ids: %s\n'%(len(zebra_ensemble_missing)))
for i in zebra_ensemble_missing:
    statout.write('%s\t%s\n'%(i,zfinid_to_name[i]))

statout.write('\n')
statout.write('***********************************************************************************************************\n')
# the ones that miss both ensemble ids and gene names
statout.write('zfin ids that miss both ensemble ids and gene names: %s\n'%(len(zebramissing)))

listwrite(zebramissing)

statout.write('\n')
statout.write('***********************************************************************************************************\n')
# the zfin ids that do not have any uberon annotations ( they have GO or other annotaions)
statout.write('the zfin ids that do not have any uberon annotations ( they have GO or other annotaions): %s\n'%(len(ids_without_uberon)))
listwrite(ids_without_uberon)



#################################### code for mouse profile generation #################################################

# Since some MGI ids are not mapping to the biomart download, i am using the MGI download file to match the remainder
mgifile = open('MRK_List2_9_18_17.rpt', 'r')

# dictionary to store mgi id to gene symbol mappings
mgiid_to_name ={}

for line in mgifile:
    if line != '\n' and  'MGI Accession ID' not in line:
        #print line
        a = line.strip().split('\t')
        #print a
        mgiid_to_name[a[0].lower()]= a[6].lower()

mgifile.close()

print 'Total number of MGI mappings to gene symbols from MGI download:',len(mgiid_to_name)

# opening the output file to write mouse profiles
mouseprof = open('mouse_anatomy_profiles.txt','wb+')

# writing the header line
mouseprof.write('mgi_id\tgene_name\tensemble_id\tuberon_annotations\n')

# a list to store mgi_ids without uberon annotations (they have GO annotations)]
mgiids_without_uberon =[]

# a list to store mgi ids that only miss ensemble ids (these have gene names)
mouse_ensemble_missing =[]

# a list to store mgi ids that  miss both ensemble ids and gene names
mousemissing =[]


# iterating through  the mgi ids from the scala output
for k,v in mgi_uberon.items():
    # this if statement is needed to remove the genes without uberon ids (the ones with only GO annotations)
    if v:
        # these are the genes that have both gene name and ensemble ids
        if k in mousemap:
            info = mousemap[k]
            mouseprof.write('%s\t%s\t%s\t'%(k,info['gene_name'],info['ensemble_id']))

        # These genes only have a gene name, no ensemble id
        elif k in mgiid_to_name:
            mouse_ensemble_missing.append(k)
            mouseprof.write('%s\t%s\tnot_found\t' % (k,mgiid_to_name[k]))

        # These genes do not have gene names and ensemble ids
        else:
            #print k
            mouseprof.write('%s\tnot_found\tnot_found\t' % (k))
            mousemissing.append(k)

        # iterating through the annotations
        for i in v:
            if i.startswith('uberon'):
                if i == v[-1]:
                    mouseprof.write('%s' % (i))
                else:
                    mouseprof.write('%s,' % (i))

        mouseprof.write('\n')
    else:

        mgiids_without_uberon.append(k)

mouseprof.close()

print 'mgi ids that do not have both ensemble ids and gene names:',len(mousemissing)

# writing the list of MGI ids that are missing some data in statistics output: gene names, ensemble ids

statout.write('********************* Mouse unmapped ********************\n')

# the ones that only miss ensemble ids
statout.write('MGI ids that only miss ensemble ids: %s\n'%(len(mouse_ensemble_missing)))
for i in mouse_ensemble_missing:
    statout.write('%s\t%s\n'%(i,mgiid_to_name[i]))

statout.write('\n')
statout.write('***********************************************************************************************************\n')
# the ones that miss both ensemble ids and gene names
statout.write('MGI ids that miss both ensemble ids and gene names: %s\n'%(len(mousemissing)))

listwrite(mousemissing)

statout.write('\n')
statout.write('***********************************************************************************************************\n')
# the MGI ids that do not have any uberon annotations ( they have GO or other annotaions)
statout.write('the MGI ids that do not have any uberon annotations ( they have GO or other annotaions): %s\n'%(len(mgiids_without_uberon)))
listwrite(mgiids_without_uberon)

statout.close()