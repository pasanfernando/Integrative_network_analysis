import networkx as nx
from collections import defaultdict
from collections import OrderedDict
from operator import itemgetter
from collections import OrderedDict
import numpy as np


#import matplotlib.pyplot as plt

""" Compares the detected modules
    Uses the orthology mapping file between mouse and zebrafish
    
    Assumption: All mouse genes have ids/ Otherwise there will be a dictionary lookup error
    
    Needs mouse gene id to symbol mapping file downloaded from MGI
    
    please change the direct_annotation_to term name with the correct term names
    
    10/5/2018 Update
    Updated the rank to increase the number after duplicated entries
    updated the zfin id mapping to get the ids from the anatomical profiles
    
    
    """

############################################################################# ###########

#define the terms
mouseterm='hindlimb'
zebraterm='pelvic_fin'

# defining multi dictionaries to map orthology information

# a multi dictionary to map zebrafish gene symbols: mouse symbols
zebradic = defaultdict(set)

# a multi dictionary to map mouse gene symbols: zebrafish symbols
mousedic = defaultdict(set)

# storing zfin and mgi ids
# zebra names: zfin id
zid ={}

# mouse names: mgi id
mid ={}


# Opening the orthology mapping file downloaded from the zfin database

orthfile = open('mouse_orthos_2018.06.26.txt', 'r')

for line in orthfile:
    if line != '\n' and line.startswith('ZDB-GENE'):
        sep = line.strip().split('\t')

        # populating the zebrafish and mouse orthology dictionaries
        zebradic[sep[1].lower()].add(sep[3].lower())
        mousedic[sep[3].lower()].add(sep[1].lower())

        # populating the zebrafish and mouse id dictionaries
        zid[sep[1].lower()]=sep[0].lower()

        mid[sep[3].lower()] = sep[5].lower()
#print zebradic['cyp26b1']
#print mousedic
#print zebradic
# print zid
# print mid

# first checking whether there are multiple mouse orthologs for zebrafish; Usually its multiple zebrafish orthologs for mouse
# for i in zebradic:
#     if len(zebradic[i])>1:
#         print i, zebradic[i]



orthfile.close()

#opening the mouse gene data file downloaded from MGI database
midfile = open('mouse_gene_ids.rpt', 'r')

# filling in the missing mouse ids for some genes
for line in midfile:
    if line != '\n' and line.startswith('MGI:'):
        sep = line.strip().split('\t')
        if sep[6].lower() not in mid:
            mid[sep[6].lower()]=sep[0].lower()

midfile.close()

#Opening anatomical profile file
mgenelist = open('newmouse_anatomy_profiles.txt', 'r')


# filling the missing mouse ids from the anatomical profile file
# reading the gene file in implementation
for line in mgenelist:
    if line != '\n':
        line = line.strip()
        a = line.split('\t') # splitting the line by tab to separate gene name vs functions
        if a[1] != 'gene_name': # excluding the header
            a[1] = a[1].lower() # convert to lower case to avoid mismatches
            if a[1] not in mid:
                mid[a[1]]=a[0]


mgenelist.close()

# opening the zebrafish gene data file downloaded from zfin database
zidfile = open('zebraids.txt', 'r')

# filling in the missing mouse ids for some genes
for line in zidfile:
    if line != '\n' and line.startswith('ZDB-'):
        sep = line.strip().split('\t')
        if sep[1].lower() not in zid:
            zid[sep[1].lower()] = sep[0].lower()

zidfile.close()

# Opening zebrafish anatomical profile file
zgenelist = open('newzebrafish_anatomy_profiles.txt', 'r')

# filling the missing mouse ids from the anatomical profile file
# reading the gene file in implementation
for line in zgenelist:
    if line != '\n':
        line = line.strip()
        a = line.split('\t')  # splitting the line by tab to separate gene name vs functions
        if a[1] != 'gene_name':  # excluding the header
            a[1] = a[1].lower()  # convert to lower case to avoid mismatches
            if a[1] not in zid:
                zid[a[1]] = a[0]
zgenelist.close()

# defining dictionaries to store the degree for each gene
mousedeg ={}
zebradeg ={}
mouswdeg ={}
zebrawdeg={}

# a method to read a module network and calculate degree for each gene

def moddegree (filename, dicname,dic2):
    # Reading the mouse module network file
    in1 = open(filename, 'r')
    G = nx.Graph()
    for line in in1:
        if line != '\n':
            #print line
            if 'combined_score' not in line:
                line = line.strip('\n')
                a = line.split()
                #print a[0],a[1],a[2]
                G.add_edge(a[0].lower(),a[1].lower(), weight = a[2])

    # calculating the degree for each gene
    for i in G.nodes():
        tot =0
        wtot=0
        nb = nx.all_neighbors(G, i)
        for j in nb:
            weightdic = G[i][j]
            weight = float(weightdic['weight'])
            tot = tot+1
            wtot=wtot+weight
        dicname[i]=tot
        dic2[i]=wtot

# calculating the degree for the mouse network
moddegree('uberon_0002103_gene_extraction.txt',mousedeg,mouswdeg)

# calculating the degree for the zebrafish network
moddegree('uberon_0000152_gene_extraction.txt',zebradeg,zebrawdeg)


# reading the predicted mouse gene file
predmouselist =[]

# reading the predicted zebrafish gene file
mousepredfile = open('hindlimb_predicted_genes.txt','r')
for line in mousepredfile:
    if line != '\n':
        #a = line.split('\t')
        predmouselist.append(line.strip().lower())

# a list to store predicted zebrafish genes
predzebralist =[]

# reading the predicted zebrafish gene file
zebrapredfile = open('pelvic_predicted_genes.txt','r')
for line in zebrapredfile:
    if line != '\n':
        #a = line.split('\t')
        predzebralist.append(line.strip().lower())


# print mousedeg
# print zebradeg

# a dictionary to store rank
#mouserank ={}
#zebrarank ={}
# mousedegsorted=OrderedDict(sorted(mousedeg.items(),key=lambda t:t[1],reverse=True))
# print mousedegsorted

# # ranking the dictionaries
# r = {key: rank for rank, key in enumerate(sorted(set(mousedeg.values()), reverse=True), 1)}
# mouserank={k: r[v] for k,v in mousedeg.items()}
# #print mouserank
#
# r = {key: rank for rank, key in enumerate(sorted(set(zebradeg.values()), reverse=True), 1)}
# zebrarank={k: r[v] for k,v in zebradeg.items()}

# a method for ranking data
def ranking(data):
    rank, count, previous, result = 0, 0, None, {}
    s_data = sorted(data.items(), key=lambda item: item[1],reverse=True)
    for key, num in s_data:
        count += 1
        if num != previous:
            rank += count
            previous = num
            count = 0
        result[key] = rank

    return result

# # ranking the dictionaries
mouserank = ranking(mousedeg)
zebrarank= ranking(zebradeg)
mousewrank = ranking(mouswdeg)
zebrawrank = ranking(zebrawdeg)

#print zebrarank

# var= None
# for i in mousedegsorted:
#
#     var=mousedegsorted[i]


#lists for not found mouse and zebrafish genes
mousenotfound =[] # mouse genes not found in zebrafish
zebranotfound=[]
commongenes =[]


# performing the module comparison

# writing the output table

out = open('module_comparison_output.txt', 'wb+')

out.write('Mouse_gene_name\tMouse_id\tPredicted_status\tNumber_of_mouse_orthologs\tDegree\tRank\tWeighted-degree\tWeighted_rank\tZebrafish_gene_name\tZfin_id\tPredicted_status\tNumber_of_zebrafish_orthologs\tDegree\tRank\tWeighted-degree\tWeighted-rank\n')

# obtaining the list of genes in the mouse module
mouselist = mousedeg.keys()

# obtaining the list of genes in the zebrafish module
zebralist = zebradeg.keys()

# a set to iterate through the remaining zebrafish genes
zebraremain = set(zebralist)



for gene in mouselist:
    # a variable to store predicted status
    mousepred = 'direct_annotation_to_the_'+mouseterm
    zebrapred = 'direct_annotation_to_the_'+zebraterm

    # checking whether mouse gene is predicted
    if gene in predmouselist:
        mousepred = 'predicted'
    if gene in mousedic:
        # If there are multiple zebrafish orthologs: check them one by one
        if len(mousedic[gene])>1:
            #print gene
            # a variable to check whether gene is there in the zebrafish list
            found = 'no'
            for i in mousedic[gene]:
                if i in zebralist:
                    # checking whether zebrafish gene is predicted
                    if i in predzebralist:
                        zebrapred = 'predicted'
                    found='yes'
                    #removing the gene from the zebra remaining set
                    zebraremain.remove(i)

                    out.write('%s\t%s\t%s\t'%(gene,mid[gene],mousepred))
                    # if there are multiple mouse orthologs
                    if len(zebradic[i])>1:
                        out.write('%s:'%(len(zebradic[i])))
                        #writing the multiple mouse orthologs
                        out.write(','.join(zebradic[i]))
                    # if there are no multiple mouse orthologs
                    elif len(zebradic[i])==1:
                        out.write('1:%s'%(gene))
                    else:
                        out.write('NA')


                    # continuing filling the table
                    out.write('\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t'%(mousedeg[gene],mouserank[gene],mouswdeg[gene],mousewrank[gene],i,zid[i],zebrapred))
                    #going through zebrafish orthologs
                    out.write('%s:' % (len(mousedic[gene])))
                    # writing the multiple zebra orthologs
                    out.write(','.join(mousedic[gene]))
                    #writing the remainder
                    out.write('\t%s\t%s\t%s\t%s\n'%(zebradeg[i],zebrarank[i],zebrawdeg[i],zebrawrank[i]))
            # if mouse genes are not found in the zebrafish list
            if found == 'no':
                if gene in zebralist:
                    print 'matched by name change code', gene
                    commongenes.append(gene)

                else:
                    mousenotfound.append(gene)
                    out.write('%s\t%s\t%s\tNA\t%s\t%s\t%s\t%s\tNA\tNA\tNA\t'%(gene,mid[gene],mousepred,mousedeg[gene],mouserank[gene],mouswdeg[gene],mousewrank[gene]))
                    # going through zebrafish orthologs
                    out.write('%s:' % (len(mousedic[gene])))
                    # writing the multiple zebra orthologs
                    out.write(','.join(mousedic[gene]))
                    # ending the line with NAs
                    out.write('\tNA\tNA\tNA\tNA\n')

            else:
                commongenes.append(gene)

        #if there is only one zebrafish ortholog
        else:
            # extracting the zebrafish ortholog
            #print gene
            # getting the first element in the set
            i = list(mousedic[gene])[0]
            #print i
            if i in zebralist:
                #print gene
                # removing the gene from the zebra remaining set
                commongenes.append(i)
                zebraremain.remove(i)

                # checking whether zebrafish gene is predicted
                if i in predzebralist:
                    zebrapred = 'predicted'

                out.write('%s\t%s\t%s\t' % (gene, mid[gene],mousepred))
                # if there are multiple mouse orthologs
                if len(zebradic[i]) > 1:
                    out.write('%s:' % (len(zebradic[i])))
                    # writing the multiple mouse orthologs
                    out.write(','.join(zebradic[i]))
                # if there are no multiple mouse orthologs
                elif len(zebradic[i]) == 1:
                    out.write('1:%s' % (gene))
                else:
                    out.write('NA')


                out.write('\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t1:%s\t%s\t%s\t%s\t%s\n' % (mousedeg[gene], mouserank[gene],mouswdeg[gene],mousewrank[gene], i, zid[i],zebrapred,i,zebradeg[i],zebrarank[i],zebrawdeg[i],zebrawrank[i]))

            else:
                if gene in zebralist:
                    print 'matched by name change code',gene

                else:
                    mousenotfound.append(gene)
                    out.write('%s\t%s\t%s\tNA\t%s\t%s\t%s\t%s\tNA\tNA\tNA\t' % (gene, mid[gene], mousepred, mousedeg[gene], mouserank[gene],mouswdeg[gene],mousewrank[gene]))
                    # going through zebrafish orthologs
                    out.write('%s:' % (len(mousedic[gene])))
                    # writing the multiple zebra orthologs
                    out.write(','.join(mousedic[gene]))
                    # ending the line with NAs
                    out.write('\tNA\tNA\tNA\tNA\n')

    # if gene has no orthologs
    else:
        #print gene
        if gene in zebralist:
            #print 'matched by name change code', gene

            # checking whether zebrafish gene is predicted
            if gene in predzebralist:
                zebrapred = 'predicted'

            commongenes.append(gene)

            # removing the gene from the zebra remaining set
            zebraremain.remove(gene)

            out.write('%s\t%s\t%s\t' % (gene, mid[gene],mousepred))
            # if there are multiple mouse orthologs
            if len(zebradic[i]) > 1:
                out.write('%s:' % (len(zebradic[gene])))
                # writing the multiple mouse orthologs
                out.write(','.join(zebradic[gene]))
            # if there are no multiple mouse orthologs
            elif len(zebradic[i]) == 1:
                out.write('1:%s' % (gene))
            else:
                out.write('NA')


            out.write('\t%s\t%s\t%s\t%s\t%s\tNA\t%s\t1:%s\t%s\t%s\t%s\t%s\n' % (mousedeg[gene], mouserank[gene],mouswdeg[gene],mousewrank[gene], gene, zebrapred,i,zebradeg[gene], zebrarank[gene],zebrawdeg[gene],zebrawrank[gene]))


        else:
            mousenotfound.append(gene)
            if gene in mid:
                out.write('%s\t%s\t%s\tNA\t%s\t%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n' % (gene, mid[gene], mousepred,mousedeg[gene], mouserank[gene],mouswdeg[gene],mousewrank[gene]))
            else:
                print gene
                out.write('%s\tNA\t%s\tNA\t%s\t%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n' % (gene, mousepred,mousedeg[gene], mouserank[gene],mouswdeg[gene],mousewrank[gene]))




print len(zebraremain)

#iterating through the remaining zebrafish genes
for i in zebraremain:
    zebrapred='direct_annotation_to_the_'+zebraterm
    # double check whether the gene is in mouse
    if i in mouselist:
        print 'this remaining zebra gene is found in mouse, please check',g
        # checking whether zebrafish gene is predicted
    if i in predzebralist:
        zebrapred = 'predicted'
    out.write('NA\tNA\tNA\t')
    # if there are multiple mouse orthologs
    if len(zebradic[i]) > 1:
        out.write('%s:' % (len(zebradic[i])))
        # writing the multiple mouse orthologs
        out.write(','.join(zebradic[i]))
    # if there are no multiple mouse orthologs


    elif len(zebradic[i]) == 1:
        out.write('1:%s' % (list(zebradic[i])[0]))
    else:
        out.write('NA')


    # continuing filling the table
    if i  in zid:
        out.write('\tNA\tNA\tNA\tNA\t%s\t%s\t%s\tNA\t%s\t%s\t%s\t%s\n' % ( i, zid[i],zebrapred,zebradeg[i], zebrarank[i],zebrawdeg[i],zebrawrank[i]))
    else:
        out.write('\tNA\tNA\tNA\tNA\t%s\tNA\t%s\tNA\t%s\t%s\t%s\t%s\n' % (i, zebrapred,zebradeg[i], zebrarank[i],zebrawdeg[i],zebrawrank[i]))
        print 'zebrafish gene without id:',i






out.close()

# statistics

print len(mouselist)

# writing comparison statistics
statfile = open('comparison_stats.txt','wb+')

statfile.write('Number of genes in the mouse module: %s\n'%(len(mouselist)))
statfile.write('Number of genes in the zebrafish module: %s\n'%(len(zebralist)))
statfile.write('Number of genes common to the two modules: %s\n'%(len(commongenes)))
statfile.write('Number of genes only found in zebrafish: %s\n'%(len(zebraremain)))
statfile.write('Number of genes only found in the mouse module: %s\n'%(len(mousenotfound)))
statfile.write('\n')

statfile.write('Number of predicted genes in the mouse module: %s\n'%(len(set(mouselist)&set(predmouselist))))
statfile.write('Number of predicted genes in the zebrafish module: %s\n'%(len(set(zebralist)&set(predzebralist))))
statfile.write('Number of predicted mouse genes common to the two modules: %s\n'%(len(set(commongenes)&set(predmouselist))))
statfile.write('Number of predicted zebradfish genes common to the two modules: %s\n'%(len(set(commongenes)&set(predzebralist))))
statfile.write('Number of predicted genes only found in zebrafish: %s\n'%(len(set(zebraremain)&set(predzebralist))))
statfile.write('Number of predicted genes only found in the mouse module: %s\n'%(len(set(mousenotfound)&set(predmouselist))))


statfile.close()