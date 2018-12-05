import networkx as nx
import matplotlib.pyplot as plt


# network integration code
# also prints the statistics of the integration as well

# load the network and store it in a dictionary

# A nested dictionary to store combined score of the two networks
# (a,b) = {net1: score1, net2: score2}
net_score_dic ={}
net_dic1 ={}
net_dic2 ={}

# this method will populate the nested net_score_dic; parameter is the name of the network. This will be used as the field for the dictionary so give a meaningful name

def net_store(netfield):


    in1 = open(netfield+'.txt', 'r')

    #netfield = net_name.strip('.txt')
    # a list to store all the genes in the network
    netgenes =[]
    for line in in1:
        if line != '\n':
            #print line
            if 'combined_score' not in line:


                templist =[]
                line = line.strip('\n')
                a = line.split()
                #print a[0],a[1],a[2]
                gene1 = a[0].lower()
                print gene1
                gene2 = a[1].lower()

                # storing the genes in the list
                netgenes.append(gene1)
                netgenes.append(gene2)

                # appending the two genes for a temporary list for sorting
                templist.append(a[0].lower())
                templist.append(a[1].lower())

                # sorting the list so always the key would be A,B; not B,A
                templist.sort()

                # if the combination is already there, no need for a new nested dic
                if templist[0]+','+templist[1] in net_score_dic:

                    net_score_dic[templist[0]+','+templist[1]][netfield]= a[2]
                # gene combination is not already there; need to create a new nested dic before storing the combination
                else:
                    net_score_dic[templist[0] + ',' + templist[1]] = {}
                    net_score_dic[templist[0] + ',' + templist[1]][netfield] = a[2]

    #converting netgenes to a set to remove duplicates
    netgenes = set(netgenes)

    return netgenes
# names of the two network files
network1 ='wang_gene_randomizedprofiles.txt'
network2= 'string.txt'

# getting stripped out versions of .txt
net_stripped1 = network1.strip('.txt')
net_stripped2 = network2.strip('.txt')


# calling the net_store function for the network 1
net1genes=net_store(net_stripped1)

#print net_score_dic

# calling the net_store function for the network 2
net2genes=net_store(net_stripped2)
#print net_score_dic

# network score integration segment

# generating the two parameters based on the accuracies/ or other metric of the two networks
# use needs to input the two metrics; metric 1 corresponds to network 1
metric1 = 0.97
metric2= 0.85
tot_metric = metric1+metric2
#generating the parameter values
para1= metric1/tot_metric
print 'para 1:', para1

para2= metric2/tot_metric
print 'para 2:', para2

# generating the integrated network

# opening the integrated network output file
int_net = open('integrated_network.txt', 'wb+')

# writing the header line
int_net.write('protein1 protein2 combined_score\n')

# integration code

#lists to store interaction pairs
commoninteractions =[]
onlynet1 =[]
onlynet2=[]
counter =0
# iterating through key value pairs of score dictionary
for k in net_score_dic:
    v=net_score_dic[k]
    counter= counter+1
    print counter

    # if combined scores are available for both networks
    if net_stripped1 in v and net_stripped2 in v:
        #print k,v
        # these are common interactions
        commoninteractions.append(k)
        # must be converted to float datatype
        ntscore1 =float(v[net_stripped1])
        ntscore2 = float(v[net_stripped2])
        # calculating the integrated network score using weighted equation
        integrated_score = (ntscore1*para1)+(ntscore2*para2)

    # If combined score is only available for network 1
    elif net_stripped1 in v:
        onlynet1.append(k)
        #print 'only 1',k, v
        ntscore1 = float(v[net_stripped1])
        # since the other network score is not available, weight*2
        integrated_score = (ntscore1 * para1)

    # If combined score is only available for network 2
    elif net_stripped2 in v:
        onlynet2.append(k)
        #print 'only 2',k, v
        ntscore2 = float(v[net_stripped2])
        integrated_score = (ntscore2*para2)

    # storing the integrated score in the integrated score field
    net_score_dic[k]['integrated_score'] = integrated_score
    gene_temp = k.split(',')
    int_net.write('%s %s %s\n'%(gene_temp[0],gene_temp[1],integrated_score))

#print net_score_dic

# clearing the dictionary to save memory
net_score_dic.clear()
# # writing the protein pare and integrated score; if you don't want decimals round the combined score to nearest integer
# for k,v in net_score_dic.items():
#     gene_temp = k.split(',')
#     int_net.write('%s %s %s\n'%(gene_temp[0],gene_temp[1],v['integrated_score']))

# writing the stat file
statfile = open('integrationstats.txt', 'wb+')
statfile.write('network1: %s\n'%(network1))
statfile.write('network2: %s\n'%(network2))
statfile.write('\n')
statfile.write('Number of genes in network1: %s\n'%(len(net1genes)))
statfile.write('Number of genes in network2: %s\n'%(len(net2genes)))
statfile.write('Number of genes only in network1: %s\n'%(len(net1genes-net2genes)))
statfile.write('Number of genes only in network2: %s\n'%(len(net2genes-net1genes)))
statfile.write('Number of common genes in both the networks: %s\n'%(len(net1genes&net2genes)))
statfile.write('\n')
statfile.write('Number of common interactions in both the networks: %s\n'%(len(commoninteractions)))
statfile.write('Number of interactions in network1: %s\n'%(len(set(commoninteractions)|set(onlynet1))))
statfile.write('Number of interactions in network2: %s\n'%(len(set(commoninteractions)|set(onlynet2))))
statfile.write('Number of interactions only in network1: %s\n'%(len(onlynet1)))
statfile.write('Number of interactions only in network2: %s\n'%(len(onlynet2)))