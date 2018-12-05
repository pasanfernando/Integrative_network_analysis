import numpy as np
import matplotlib.pyplot as plt
#humanneighbourcount.txt
#zebrafish_neighbourcount.txt
#human_noncandidates.txt
# input1 = raw_input('enter the interaction file:')
# input2 = raw_input('enter the second interaction file name:')
in1 = open('candidatesneighbourcount.txt', 'r')
in2 = open('noncandidatesneighbourcount.txt', 'r')

y =[]
z = []
y1 = []
z1 = []
for line in in1:
    if line != '\n':
        line = line.strip()
        if 'gene_name' not in line:
            a = line.split()
            y.append(int(a[1]))
            z.append(float(a[3]))


for line in in2:
    if line != '\n':
        line = line.strip()
        if 'gene_name' not in line:
            a = line.split()
            y1.append(int(a[1]))
            z1.append(float(a[3]))

# comparison of the direct gene counts

figs, axes =plt.subplots(nrows=1, ncols=2, sharey=True, squeeze=True,figsize=(12,8))
figs.suptitle('The boxplot comparison of the module gene counts\n',fontsize=15)
# setting custom x axis labels
plt.margins(0.01)
#plt.ylabel('Candidate gene count')
axes[0].set_ylabel('Module gene count', fontsize=15)
axes[0].boxplot(y, showmeans=True)
axes[0].set_xlabel('Module', fontsize=10)
plt.setp(axes[0].get_xticklabels(), visible=False)
axes[1].boxplot(y1, showmeans=True)
axes[1].set_xlabel('Network background', fontsize=10)
plt.setp(axes[1].get_xticklabels(), visible=False)

plt.savefig('boxplotcomparison.jpg',dpi=1200)
plt.close()

# plotting the histogram comparison
plt.margins(0.01)
#figs, axes = plt.subplots(nrows=2, ncols=1, sharey=True, squeeze=True,figsize=(9,11),dpi=600)
figs, axes = plt.subplots(2, 1, sharex=True)
plt.margins(0.01)
plt.subplots_adjust(hspace=0.25)
figs.suptitle('The histogram comparison of the candidate gene counts\n',fontsize=15)
#bins = np.linspace(0, 1, 100)
# plotting the candidate distribution
axes[0].hist(y, 100,color='red')
axes[0].set_title('Module',fontsize=10)
#setting a custom y limit because candidate genes are low
axes[0].set_ylim(0,50)

#axes[0].set_xlim(0,1)
# plotting the non-candidate distribution
axes[1].hist(y1, 100,color='blue')
axes[1].set_title('Network background',fontsize=10)
axes[1].set_ylim(0,10000)
#axes[1].set_xlim(0,1)

plt.legend()
plt.savefig('histcomparison.jpg',dpi=1200)
plt.close()

############ Comparison of the normalized frequencies ##############

figs, axes =plt.subplots(nrows=1, ncols=2, sharey=True, squeeze=True,figsize=(12,8))
figs.suptitle('The boxplot comparison of the normalized candidate gene counts\n',fontsize=15)
# setting custom x axis labels
plt.margins(0.01)
#plt.ylabel('Candidate gene count')
axes[0].set_ylabel('Normalized  andidate gene count', fontsize=15)
axes[0].boxplot(z, showmeans=True)
axes[0].set_xlabel('Module', fontsize=10)
plt.setp(axes[0].get_xticklabels(), visible=False)
axes[1].boxplot(z1, showmeans=True)
axes[1].set_xlabel('Network background', fontsize=10)
plt.setp(axes[1].get_xticklabels(), visible=False)

plt.savefig('normboxplotcomparison.jpg',dpi=1200)
plt.close()

# plotting the histogram comparison
plt.margins(0.01)
#figs, axes = plt.subplots(nrows=2, ncols=1, sharey=True, squeeze=True,figsize=(9,11),dpi=600)
figs, axes = plt.subplots(2, 1, sharex=True)
plt.margins(0.01)
plt.subplots_adjust(hspace=0.25)
figs.suptitle('The histogram comparison of the normalized candidate gene counts\n',fontsize=15)
#bins = np.linspace(0, 1, 100)
# plotting the candidate distribution
axes[0].hist(z, 100,color='red')
axes[0].set_title('Module',fontsize=10)
#setting a custom y limit because candidate genes are low
axes[0].set_ylim(0,50)

#axes[0].set_xlim(0,1)
# plotting the non-candidate distribution
axes[1].hist(z1, 100,color='blue')
axes[1].set_title('Network background',fontsize=10)
axes[1].set_ylim(0,10000)
#axes[1].set_xlim(0,1)

plt.legend()
plt.savefig('normhistcomparison.jpg',dpi=1200)
plt.close()