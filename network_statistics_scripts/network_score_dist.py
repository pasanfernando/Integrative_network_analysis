'''Plotting the genesimilarity score distribution for a given network
'''

import collections
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as ss
from collections import defaultdict


netfile = open('linstrweightedintegrated_network.txt','r')

#list to store the score distribution
scordist =[]

for line in netfile:
    if line != '\n' and 'combined_score' not in line:
        a = line.strip().split(' ')
        scordist.append(float(a[2]))


# histogram for the auc of number of gene annotations curve
plt.hist(scordist,1000, color='red')
plt.xlabel('Gene similarity score')
plt.ylabel('Frequency')
plt.title('Histogram for gene similarity scores for the zebrafish network\n')
#plt.show()
plt.savefig('sshist.jpg')
plt.close()

# cumulative histogram
plt.hist(scordist, 1000,  normed =1,histtype='step',cumulative=True)
plt.xlabel('Gene similarity score')
plt.ylabel('Cumulative frequency')
plt.title('Cumulative frequency histogram for gene similarity scores for the integrative network\n')
#plt.show()
plt.savefig('ssfhist.jpg')
plt.close()

#box plot
plt.boxplot(scordist)
plt.title('Frequency histogram for gene similarity scores for the integrative network\n')
plt.savefig('ssbox.jpg')
plt.close()

print 'number of interactions',len(scordist)