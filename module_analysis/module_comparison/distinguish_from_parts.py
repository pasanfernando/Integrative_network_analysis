import networkx as nx
from collections import defaultdict
from collections import OrderedDict
from operator import itemgetter
from collections import OrderedDict
import numpy as np


#import matplotlib.pyplot as plt

""" Further modifies the moducle comparison output table to distinguish the orignal annotations or original annotation to part
    
    Also gives the statistics. Such as common genes, only found in either list, etc.
    
    """

############################################################################# ###########

# lists to store genes annotated to parts
mousepartsgenes =[]
zebrapartsgenes =[]

zpartgfile = open('pelvic_partgenes.txt','r')

for line in zpartgfile:

    if line !='\n':
        #print line
        zebrapartsgenes.append(line.strip())

print zebrapartsgenes

mpartgfile = open('hindlimb_partgenes.txt','r')

for line in mpartgfile:

    if line !='\n':
        #print line
        mousepartsgenes.append(line.strip())

print mousepartsgenes

#list to count matching genes
zmatches=[]
mmatches=[]

# lists to store common genes
zebracommon =[]
mousecommon =[]

#lists to confirm the predicted genes
mousepredpred =[]
mousepredori =[]
zebrapredpred =[]
zebrapredori=[]

# genes that are only in one organism
zebraonly =[]
mouseonly =[]

comfile = open('module_comparison_output.txt', 'r')

out = open('pel_hind_comparisontable.txt','wb+')
out.write('Mouse_gene_name\tMouse_id\tPredicted_status\tNumber_of_mouse_orthologs\tDegree\tRank\tWeighted-degree\tWeighted_rank\tZebrafish_gene_name\tZfin_id\tPredicted_status\tNumber_of_zebrafish_orthologs\tDegree\tRank\tWeighted-degree\tWeighted-rank\n')

for line in comfile:
    if line !='\n' and 'Mouse_gene_name' not in line:
        a = line.strip().split('\t')
        if a[0] in mousepartsgenes:
            a[2]='annotated_to_a_part_or_bud'
            mmatches.append(a[0])
        if a[8] in zebrapartsgenes:
            a[10]='annotated_to_a_part_or_bud'
            zmatches.append(a[6])

        if a[0]!='NA' and a[8]!= 'NA':
            zebracommon.append(a[8])
            mousecommon.append(a[0])

            if a[2]=='predicted' and a[10]=='predicted':
                zebrapredpred.append(a[8])
                mousepredpred.append(a[0])

            elif a[2]=='predicted' and a[10]!='predicted':
                mousepredori.append(a[0])

            elif a[2]!='predicted' and a[10]=='predicted':
                zebrapredori.append(a[8])

        if a[0]!='NA' and a[8]== 'NA':
            mouseonly.append(a[0])

        if a[0]== 'NA' and a[8] != 'NA':
            zebraonly.append(a[8])

        out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],a[10],a[11],a[12],a[13],a[14],a[15]))


out.close()
comfile.close()

print 'number of mouse genes annotated to parts that were replaced',len(mmatches)
print 'number of zebrafish genes annotated to parts that were replaced',len(zmatches)

#opening the statistics file

stat = open('comparison_postprocessed_stats.txt','wb+')
stat.write('Number of mouse genes that are common: %s\n'%(len(set(mousecommon))))
stat.write('Number of zebrafish genes that are common: %s\n'%(len(set(zebracommon))))
stat.write('\n')

stat.write('Number of predicted mouse genes that are also predicted in zebrafish: %s\n'%(len(set(mousepredpred))))
stat.write('Number of predicted zebrafish genes that are also predicted in mouse: %s\n'%(len(set(zebrapredpred))))
stat.write('Number of predicted mouse genes that are not predicted in zebrafish: %s\n'%(len(set(mousepredori))))
stat.write('Number of predicted zebrafish genes that are not predicted in mouse: %s\n'%(len(set(zebrapredori))))
stat.write('\n')

stat.write('Number of of genes only found in mouse module: %s\n'%(len(set(mouseonly))))
stat.write('Number of of genes only found in zebrafish module: %s\n'%(len(set(zebraonly))))
stat.write('\n')

stat.write('Number of mouse genes annotated to parts that were replaced: %s\n'%(len(set(mmatches))))
stat.write('Number of zebrafish genes annotated to parts that were replaced: %s\n'%(len(set(zmatches))))

stat.close()