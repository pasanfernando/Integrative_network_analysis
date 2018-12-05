import networkx as nx
from collections import defaultdict


'''
Extract set of genes for the given function from the network based on a given threshold
Only the new genes higher than the threshold will be selected 
these genes must be interconnected together to be selected 

You can give a custom gene list as the input or the full profiles file 
'''

# the function that the prediction is done
function = 'uberon_0002103'

# the threshold value for the prediction
threshold = 300


# creating a nested dictionary to store chi square values for different functions for each gene
chigenes = defaultdict(dict)
#Opening anatomical profile file
genelist = open('newmouse_anatomy_profiles.txt', 'r')

# defining a multiple dictionary to store functions: genes
genefunc = defaultdict(list)
# defining a multiple dictionary to store gene: functions
genedic = defaultdict(list)

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

print 'number of functions in the profile:',len(genefunc)
print 'number of genes in the original profile:',len(genedic)
# a dictionary to store predicted genes : score value
predgenes ={}
# getting the predicted candidate genes
# reading the chigenes dictionary from predicted functions file
# a list to store chivalues; used to get the maximum chi value
chivalues = []
############################################################################################
# a dictionary to find the maximum chi value for each function

# store all the genes and threshold score for the given function
chimaxfunc = {}
chidicfile = open('newfullhindlimbincludingparts_functions.txt', 'r')
for line in chidicfile:
    if line != '\n':
        line = line.strip()
        firstsplit = line.split('\t')
        if firstsplit[0] != 'genes':
            chigenes[firstsplit[0]]={}
            secondsplit= firstsplit[1].split(',')
            #print secondsplit
            for i in secondsplit:
                #print i
                thirdsplit = i.split(':')
                chigenes[firstsplit[0]][thirdsplit[0]]=float(thirdsplit[1])
                # if the predicted score is higher or equal to the threshold store that gene
                if thirdsplit[0] == function:
                    if float(thirdsplit[1])>= threshold:
                        predgenes[(firstsplit[0])] = float(thirdsplit[1])
                    chimaxfunc[firstsplit[0]]=float(thirdsplit[1])
                #chivalues.append(float(thirdsplit[1]))


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
    # reading the genes from a custom list if the option is given
    orgenelist = clist
# Otherwise read from the original profile
else:
    orgenelist = genefunc[function]



# geting the final gene list by adding original and final gene list
genelist = set(orgenelist)|set(predgenes.keys())

onlypredicted = set(predgenes.keys()) - set(orgenelist)
# a list to store selected genes
selectedgenes =[]

# Opening the original network

# storing genes in the network
netgenes =[]

# counter to count the interactions in the extracted network
intcounter=0

in1 = open('mslick0.3integrated_network.txt', 'r')
# opening the preprocessed network output file
pre_net = open(function+'_lin4gene_extraction.txt', 'wb+')
# writing the header line
pre_net.write('protein1 protein2 combined_score\n')
for line in in1:
    if line != '\n':
        # print line
        if 'combined_score' not in line:
            templist = []
            line = line.strip('\n')
            a = line.split()
            # print a[0],a[1],a[2]
            gene1 = a[0].lower()
            # print gene1
            gene2 = a[1].lower()
            netgenes.append(gene1)
            netgenes.append(gene2)


            # checking wether both genes are in the newtork
            if gene1 in genelist and gene2 in genelist:
                pre_net.write('%s %s %s\n' % (gene1, gene2, a[2]))
                selectedgenes.append(gene1)
                selectedgenes.append(gene2)
                intcounter+=1


pre_net.close()
in1.close()

unselected = set(genelist)- set(selectedgenes)

#not in the network
notinnetwork = set(unselected)-set(netgenes)

alonegenes = set(unselected)&set(netgenes)

print 'total number of original genes:',len(orgenelist)
print 'number of selected genes:',len(set(selectedgenes))
print 'number of unselected genes:',len(unselected)
print unselected
print 'number of genes not in the network',len(notinnetwork)
print notinnetwork
print 'number of genes in the network but alone:',len(alonegenes)
print alonegenes

#opening the statistics file
statfile = open('extractionstats.txt','wb+')
statfile.write('Total number of original genes: %s\n'%(len(set(orgenelist))))
counter =0
for i in orgenelist:
    if i in chimaxfunc:
        if chimaxfunc[i] >= threshold:
            counter+=1
        statfile.write('%s\t%s\n'%(i, chimaxfunc[i]))
    else:
        statfile.write('%s\tNA\n' % (i))
statfile.write('Number of original genes higher than the %s threshold: %s\n'%(threshold,counter))
statfile.write('\n')
statfile.write('Total number of only predicted genes: %s\n'%(len(onlypredicted)))
for i in onlypredicted:
    statfile.write('%s\t%s\n'%(i,predgenes[i]))
statfile.write('\n')
statfile.write('Number of selected genes: %s\n'%(len(set(selectedgenes))))
statfile.write('Number of interactions in the extracted network: %s\n'%(intcounter))
statfile.write('Total number of genes: %s\n'%(len(genelist)))
statfile.write('Number of unselected genes: %s\n'%(len(unselected)))
statfile.write('Number of genes not in the network: %s\n'%(len(notinnetwork)))
for i in unselected:
    statfile.write('%s\n'%(i))

statfile.write('\n')
statfile.write('Number of genes in the network but alone: %s\n'%(len(alonegenes)))
for i in alonegenes:
    statfile.write('%s\n'%(i))

statfile.write('\n')

# writing only the selected genes in a separate file for
modfile = open('modulegenes.txt','wb+')

for i in set(selectedgenes):
    modfile.write('%s\n'%(i))

modfile.close()

# writing only the predicted genes in a separate file for comparison with genes extracted from parts
predfile = open('only_predicted_genes.txt','wb+')

for i in onlypredicted:
    predfile.write('%s\n'%(i))

predfile.close()