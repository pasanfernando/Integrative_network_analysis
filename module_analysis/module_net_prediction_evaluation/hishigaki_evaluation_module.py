import networkx as nx
from collections import defaultdict
from operator import itemgetter
from collections import OrderedDict
import numpy as np

import curve_plotter_single_function as cp

#import matplotlib.pyplot as plt

""" This code is for module prediction phase
    

    Code to gene function prediction based on the hishigaki methods
    This code is for evaluation using leave one out method
    The evaluation method for multiple classes/ functions is newly proposed
    
    Note: if you feel there is some thing wrong with the code, please uncomment the print statements throughout 
    the code. It will print results for each step. Do this with the prototype.
    
    There are two options to read the gene list
        1. Original annotations for the anatomical profile file
        2. list for custom gene list from a given file
    """

############################################################################# ###########


#methods for calculating evaluation metrics

def evaluation_metrics(p,n,tp,tn,fp,fn):

    # calculating the accuracy
    try:
        accuracy = (tp+tn)/float(p+n)
    except ZeroDivisionError:
        accuracy = 0
    #print 'Accuracy is:',accuracy

    #error rate
    try:
        error_rate = (fp+fn)/float(p+n)
    except ZeroDivisionError:
        error_rate = 0
    #print 'Error rate is:', error_rate

    #sensitivity/ same as recall
    try:
        sensitivity = tp/float(p)
    except ZeroDivisionError:
        sensitivity = 0
    #print 'sensitivity:', sensitivity

    #specificity1
    try:
        specificity = tn/float(n)
    except ZeroDivisionError:
        specificity = 0
    #print 'specificity:', specificity

    #precision
    try:
        precision = tp/float(tp+fp)
    except ZeroDivisionError:
        precision=0
    #print 'precision:', precision

    #true positive rate
    try:
        tpr = tp/float(p)
    except ZeroDivisionError:
        tpr=0

    # false positive rate
    try:
        fpr = fp/float(n)
    except ZeroDivisionError:
        fpr=0

    # false positive rate

    return accuracy,error_rate,sensitivity,specificity,precision,tpr,fpr


########################################################################################

#input1 = raw_input('enter the interaction file:')

# prompting the user for gene list selection method
input1 = raw_input('Extract a custom gene list from a file (y or n):')

if input1 == 'y':
    # a list to store the genes from the custom gene list
    clist =[]

    # if yes please check the name of the file
    cusgenelist = open('hindlimbgenelist_including_parts.txt', 'r')

    for line in cusgenelist:
        if line != '\n':
            line = line.strip()
            clist.append(line)

    cusgenelist.close()

# creating a nested dictionary to store chi square values for different functions for each gene
chigenes = defaultdict(dict)

#Reading the input gene function annotation file

genelist = open('newmouse_anatomy_profiles.txt', 'r')

# defining a multiple dictionary to store functions: genes
genefunc = defaultdict(list)
# defining a multiple dictionary to store gene: functions
genedic = defaultdict(list)

# Then convert the above input file into function: gene format

# this indicator is required for leave one out method

# reading the gene file
# for line in genelist:
#     if line != '\n':
#         line = line.strip()
#         a = line.split('\t') # splitting the line by tab to separate gene name vs functions
#         if a[0] != 'genes': # excluding the header
#             a[0] = a[0].lower() # convert to lower case to avoid mismatches
#             #print a[0]
#             # here, we are evaluating, so only the proteins with functions will be considered
#             if len(a)>1: # ths is required because there are proteins without one single function
#                 b = a[1].split(',') # seperating the functions by comma
#                 # store each function in the multi dic by function as the key and genes as values
#                 for i in b:
#                     genedic[a[0]].append(i) # this dictionary stores genes as keys
#                     genefunc[i].append(a[0])

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

# print genedic
# print len(genefunc)

### reading the network file########################################################
in1 = open('mslick0.3integrated_network.txt', 'r')
G = nx.Graph()
for line in in1:
    if line != '\n':
        #print line
        if 'combined_score' not in line:
            line = line.strip('\n')
            a = line.split()
            #print a[0],a[1],a[2]
            G.add_edge(a[0].lower(),a[1].lower(), weight = a[2])

# visualizing the network; please turn this off when loading large networks
# nx.draw_networkx(G, with_labels=True)
# plt.draw()
# plt.show()

######### calculate the function frequency or phi #########
# total proteins in the network
totps = nx.number_of_nodes(G)

# total proteins for the given function
# funps = len(genelist) # you are assuming all genes matches the ones in the network; need to check one by one

 # customized gene function list
customfunctions = ['uberon_0002103']
#customfunctions = ['uberon_0005724']

genes_in_network = nx.nodes(G)
print 'genes in the network:',genes_in_network

chvals =[]
# gene counter for runtime purposes
genecounter =0

genes_not_found = []
for gene in genedic:
    print 'gene is:',gene
    genecounter +=1
    print genecounter

    #dictionaries to store total and function counts, chi square value
    total={}
    ccount={}
    chisq = {}

    #counting the genes with given function in immediate neighborhood of gene in consideration of leave one out method

    # going through the different functions in the function list
    for i in customfunctions:
    #for i in genefunc:
        #print 'function is:',i
        # reading the genes from a custom list if the option is given
        if input1 == 'y':
            genelist = clist
        else:
            genelist = genefunc[i]
        #print genelist

        # making sure the gene in question is removed from the gene list because we are assuming it does not have any function
        if gene in genelist:
            # removing the gene by list comprehension; otherwise it will mutate the original dictionary
            genelist = [x for x in genelist if x != gene]
        #print genelist
        #print genefunc[i]

        # making sure all the genes in the input data file are in the network

        notfound = set(genelist) - set(genes_in_network)  # genes that are not found in the network
        #print 'genes that are not found in the network:', notfound
        found = set(genelist) & set(genes_in_network)  # genes that are found in the network
        #print 'genes that are found in the network:', found


        funps = len(found)

        fracps = float(funps) / totps

        #print fracps

        # counting the genes with given function in immediate neighborhood of gene in consideration in leave one out method
        totalc =0 # initially, total count is zero
        ccountc=0 # initially, neighborhood count is zero
        if gene in G:
            nb = nx.all_neighbors(G, gene)
            #print i
            for j in nb:
                #print j
                totalc=totalc+1 # counting the total neighbors; could have done this using len(nb) as well
                if j in genelist: # if the neighbor is in the gene list that neghbor has the function
                    ccountc=ccountc+1

            # storing the above calculated values in designated dictionaries
            total[i]=totalc
            ccount[i]=ccountc

            #calculate expected frequency values
            ef = fracps*totalc
            # calculating the chi-square value
            try: # this error handling is required to avoid zero division when ef=0
                chi = ((ccountc-ef)**2)/ef
            except ZeroDivisionError:
                chi = 0
            chisq[i] = chi
            chvals.append(chi)
            chigenes[gene] = chisq
        else:
            genes_not_found.append(gene)

    #print chisq
    #print total
    #print ccount

print chigenes
# print chigenes['p1']
print chigenes

# getting the min max of chvals
minv = min(chvals)
maxv=max(chvals)

# ordering each dictionary based on the highest function

for i in chigenes:
    chigenes[i]= OrderedDict(sorted(chigenes[i].items(),key=lambda t:t[1], reverse=True))
    #print chigenes[i]

# Normalizing the chi values
# for i in chigenes:
#
#     func_chi = chigenes[i]
#
#     #print lastelement[-1]
#     for j in func_chi:
#         #print j
#         #performing the normalization
#         func_chi[j]= (float(func_chi[j])-minv)/float(maxv-minv)
#
#     chigenes[i]=func_chi
#
# print chigenes

# a list to store chivalues; used to get the maximum chi value
chivalues = []
############################################################################################
# a dictionary to find the maximum chi value for each function

chimaxfunc = defaultdict(list)

# writing the predicted gene functions in an output file in sorted order
out = open('predicted_functions.txt', 'wb+')
out.write('genes\tpredicted_functions\n')

for i in chigenes:
    out.write('%s\t'%(i))
    func_chi = chigenes[i]
    lastelement = list(func_chi.keys())
    #print lastelement[-1]
    for j in func_chi:
        #print j
        # loading the chi values for each function to find the maximum chi value for each function
        chimaxfunc[j].append(func_chi[j])
        chivalues.append(float(func_chi[j]))# appending chi values
        if j == lastelement[-1]:# this avoids printing a comma after the last element
            out.write('%s:%s' % (j, func_chi[j]))
        else:
            out.write('%s:%s,'%(j,func_chi[j]))

    out.write('\n')
out.close()

# printing the chi max dictionary
print chimaxfunc

###########################################################################################
# module to calculate FP,TP,TN, FN



# getting total functions as a set
totfunc = set(genefunc.keys())
#print 'all the functions:',totfunc

# defining the dictionary for function annotation statistics
funcstat ={}


# a list to find the last element in the gene function list
func_last_element = list(genefunc.keys())

# a dictionary to draw all roc curves in the same plot
generalroc = defaultdict(dict)

# a dictionary to draw all precision-recall curves in the same plot
general_prec_recall = defaultdict(dict)

# counter for functions
funccounter = 0
# In standard one function evaluation mode, its required to iterate through function list first
for j in customfunctions:
#for j in genefunc:
    funccounter+=1
    print 'function is:',j
    print funccounter
    #It can iterate through all the functions or user can give set of functions to select from
    # in user specified case define a string of functions and iterate through them
    # defining dictionaries to store evaluation metrics for each threshold value
    accuhash ={}
    errorhash ={}
    sensehash ={}
    specifichash ={}
    precisionhash ={}

    # one multi dictionary to store fpr and tpr for thresholods
    # the key should be the threshold
    rochash = defaultdict(list)

    # another hash for prescision recall curve
    precision_recall_hash = defaultdict(list)

    #print 'the function is',j

    # then, for each function, we need to move through different threshold values
    # numpy is used to iterate through float values like 0.5
    # In this method the maximum chi value is the maximum of all functions; for function specifice one use the other line
    #for threshold in np.arange(0,(max(chivalues)+0.2),0.2):

    # function specific maximum chi value
    chi_max = max(chimaxfunc[j])

    print 'maximum chi value for function',j,':', chi_max
    print chimaxfunc[j]
    # THIS segment is required to incriment the threshold slowly at first stage then faster in the second stage
    # This is to get a curvy ROC while not wasting time
    initialrange = np.arange(0, 2, 0.01)
    finalrange = np.arange(2, (chi_max + 0.2), 0.2)
    thresholdrange = np.append(initialrange,finalrange)
    for threshold in thresholdrange:
    #for threshold in np.arange(0,(chi_max+0.2),0.2):
    # following line negates the function specific maximum chi value uncomment to use it
    #for threshold in np.arange(0, (max(chivalues) + 0.1), 1):
        #print 'the threshold is:',threshold
        #threshold use as the cutoff to draw the roc curve
        #threshold = 0

        # defining the variables
        pcountfinal =0
        ncountfinal =0
        tpcountfinal =0
        tncountfinal =0
        fpcountfinal =0
        fncountfinal =0

        # calculating the confusion matrix: tp,fp,tn,fn
        for i in chigenes:
            pred_functions =[] # list to store functions higher than the threshold
            genefunctions = chigenes[i] # we are assuming that all the funtions are in the chigenes dictionary
            #print genefunctions

            # if the function is in actual functions list for the specific gene it is a actual positive
            if j in genedic[i]:
                pcountfinal = pcountfinal + 1
                # positive can be a true positive if it is predicted or false negative if else
                if float(genefunctions[j]) >= threshold: # if the score is higher than the threshold it is predicted
                    tpcountfinal = tpcountfinal + 1

                # In this case, this would be a false negative not a false postive
                else:
                    fncountfinal = fncountfinal + 1


            # if the function is not in the actual functions list for the specific gene it is a actual negative
            else:
                ncountfinal = ncountfinal + 1
                # negative can be a true negative if it is not predicted or false positive otherwise
                if float(genefunctions[j]) >= threshold:
                    # here, it is predicted as a positive but actually it is a negative
                    fpcountfinal = fpcountfinal + 1

                else:
                    tncountfinal = tncountfinal + 1

            # the following statements adds all the variables for the given threshold
             # adding to the parent variable that sums positive count for all the genes







        # print 'P:', pcountfinal
        # print 'N:',ncountfinal
        # print 'TP:',tpcountfinal
        # print 'TN:',tncountfinal
        # print 'FP:',fpcountfinal
        # print 'FN:',fncountfinal

        # getting all the evaluation metrics including the tpr and fpr
        accuracy,error_rate,sensitivity,specificity,precision,tpr,fpr= evaluation_metrics(pcountfinal,ncountfinal,tpcountfinal,tncountfinal,fpcountfinal,fncountfinal)

        # storing metrics in defined dictionaries
        accuhash[threshold]=accuracy
        errorhash[threshold]=error_rate
        sensehash[threshold]=sensitivity
        specifichash[threshold]=specificity
        precisionhash[threshold]=precision

        rochash[threshold].append(fpr)
        rochash[threshold].append(tpr)

        #if precision!=0:
        precision_recall_hash[threshold].append(sensitivity)
        precision_recall_hash[threshold].append(precision)


    #print 'roc hasH:',rochash
    #print 'precision-recall hash:', precision_recall_hash

    # appending roc data for each function in the general ROC hash
    generalroc[j]=rochash

    # appending precision-recall data for each function in the general ROC hash
    general_prec_recall[j] = precision_recall_hash

    #print max(chivalues)



    rauc= cp.ROCplotter(rochash,j,'roc')

    # plotting the precision-recall curve for given function
    pauc = cp.ROCplotter(precision_recall_hash, j, 'precision_recall')

    # creating the nested dictionary for function statistics
    funcstat[j] = {}
    funcstat[j]['rauc']= rauc
    funcstat[j]['pauc'] =pauc
    # selecting based on custum gene list input
    if input1 == 'y':
        funcstat[j]['num_functions'] = len(clist)
        matchnum = set(clist) - set(genes_not_found)
        funcstat[j]['matched_num_functions'] = len(matchnum)
    else:
        funcstat[j]['num_functions'] = len(genefunc[j])
        matchnum = set(genefunc[j]) - set(genes_not_found)
        funcstat[j]['matched_num_functions'] = len(matchnum)


    #print accuhash
    cp.metric_plotter(accuhash,'Accuracy',j)
    cp.metric_plotter(errorhash,'Error_rate',j)
    cp.metric_plotter(sensehash,'Sensitivity',j)
    cp.metric_plotter(specifichash,'Specificity',j)
    cp.metric_plotter(precisionhash,'Precision',j)

    # this code is needed to compare the precision curves of two functions

    preciout = open(j+'precision_hash.txt','wb+')
    for i in precisionhash:
        preciout.write('%s\t%s\n'%(i,precisionhash[i]))
    preciout.close()

    # writing the ROC curve values in a file for each function
    preciout = open(j+'precision_recall_hash.txt', 'wb+')
    rocout = open(j+'roc_hash.txt', 'wb+')

    for i in rochash:
        a = rochash[i]
        rocout.write('%s\t%s\t%s\n' % ( a[0], a[1], i))

    for i in precision_recall_hash:
        a = precision_recall_hash[i]
        preciout.write('%s\t%s\t%s\n' % ( a[0], a[1],i))


# # this general roc plot can be used to find the best function with greatest ROC curve
cp.general_ROCplotter(generalroc,'roc')
#
# # general precision-recall curves for all functions
cp.general_ROCplotter(general_prec_recall,'precision_recall')

print 'number of genes not found in the network:',len(set(genes_not_found))

# writing the function stats
out1 = open('function_stats.txt', 'wb+')
out1.write('function\tnumber_of_genes\tmatched_number_of_genes\tROC_auc\tPR_auc\n')

# print funcstat

for i in funcstat:
    out1.write('%s\t%s\t%s\t%s\t%s\n'%(i,funcstat[i]['num_functions'],funcstat[i]['matched_num_functions'],funcstat[i]['rauc'],funcstat[i]['pauc']))