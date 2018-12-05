
import matplotlib.pyplot as plt
from operator import itemgetter
from collections import OrderedDict
import numpy as np
from sklearn import metrics
from collections import defaultdict


deg1=[]
deg2=[]

dfile = open('forelimb_pred.txt','r')

for line in dfile:
    if line !='\n' and 'predicted' not in line:
        a= line.strip().split('\t')
        if a[0]!= 'NA':
            deg1.append(float(a[0]))
        if a[1] != 'NA':
            deg2.append(float(a[1]))
if max(deg1)>max(deg2):
    maxval = max(deg1)

else:
    maxval=max(deg2)

# setting the xlabels
x1='Originally annotated'
x2='Predicted'
tname= 'forelimb'

plt.title('The histogram comparison of weighted degree distributions\n')
plt.ylabel('Frequency')
plt.xlabel('Degree')
#plt.figure(num=1,figsize=(10,100),dpi=300)
plt.margins(0.01)
figs, axes = plt.subplots(nrows=2, ncols=1, sharey=True, squeeze=True,figsize=(12,11),dpi=1200)
# adjusting the space between the subplots
plt.subplots_adjust(hspace=0.25)
figs.suptitle('The histogram comparison of weighted degree distributions for originally annotated vs predicted genes for the '+tname+'\n',fontsize=15)
bins = np.linspace(0, maxval, 50)
#bins=10

axes[0].hist(deg1, bins,color='red')
axes[0].set_title(x1,fontsize=10)
axes[0].set_xlim(0,maxval)
#axes[0].set_ylim(0,12)

axes[1].hist(deg2, bins,color='blue')
axes[1].set_title(x2,fontsize=10)
axes[1].set_xlim(0,maxval)


plt.legend()
plt.savefig('histcomparison.jpg',dpi=1200)
plt.close()

# plotting the boxplots
plt.ylabel('AUC')
figs, axes =plt.subplots(nrows=1, ncols=2, sharey=True, squeeze=True,figsize=(12,8))
figs.suptitle('The boxplot comparison of weighted degree distributions for originally annotated vs predicted genes for the ' +tname+'\n',fontsize=15)
# setting custom x axis labels
plt.margins(0.01)


axes[0].boxplot(deg1, showmeans=True)
axes[0].set_xlabel(x1, fontsize=10)
axes[0].set_ylabel('Weighted degree', fontsize=10)
plt.setp(axes[0].get_xticklabels(), visible=False)

#axes[0].set_title('Default')
#axes[0].boxplot(y, showmeans=True)
# plt.subplot(122)
axes[1].boxplot(deg2, showmeans=True)
axes[1].set_xlabel(x2, fontsize=10)
plt.setp(axes[1].get_xticklabels(), visible=False)

#plt.legend()
plt.savefig('boxplotcomparison.jpg',dpi=1200)
plt.close()