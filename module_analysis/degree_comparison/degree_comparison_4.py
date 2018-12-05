
import matplotlib.pyplot as plt
from operator import itemgetter
from collections import OrderedDict
import numpy as np
from sklearn import metrics
from collections import defaultdict


deg1=[]
deg2=[]
deg3=[]
deg4=[]

dfile = open('pecfore_degrees4.txt','r')

for line in dfile:
    if line !='\n' and 'mouse' not in line:
        a= line.strip().split('\t')
        if a[0]!= 'NA':
            deg1.append(float(a[0]))
        if a[1] != 'NA':
            deg2.append(float(a[1]))
        if a[2] != 'NA':
            deg3.append(float(a[2]))
        if a[3] != 'NA':
            deg4.append(float(a[3]))
mlist = [max(deg1),max(deg2),max(deg3),max(deg4)]
maxval = max(mlist)

# setting the xlabels


x1='Pectoral fin module-specific genes'
x2='Common pectoral fin module genes'
x3='Common forelimb module genes'
x4='Forelimb module-scpecific genes'

plt.title('The histogram comparison of normalized weighted-degree distributions\n')
plt.ylabel('Frequency')
plt.xlabel('Degree')
#plt.figure(num=1,figsize=(10,100),dpi=300)
plt.margins(0.01)
figs, axes = plt.subplots(nrows=4, ncols=1, sharey=True, squeeze=True,figsize=(9,11),dpi=1200)
# adjusting the space between the subplots
plt.subplots_adjust(hspace=0.25)
figs.suptitle('The histogram comparison of normalized weighted-degree distributions for common genes\n',fontsize=15)
bins = np.linspace(0, maxval, 50)
#bins=10

axes[0].hist(deg1, bins,color='red')
axes[0].set_title(x1,fontsize=10)
axes[0].set_xlim(0,maxval)
#axes[0].set_ylim(0,200)

axes[1].hist(deg2, bins,color='blue')
axes[1].set_title(x2,fontsize=10)
axes[1].set_xlim(0,maxval)

axes[2].hist(deg3, bins,color='green')
axes[2].set_title(x4,fontsize=10)
axes[2].set_xlim(0,maxval)

axes[3].hist(deg4, bins,color='orange')
axes[3].set_title(x4,fontsize=10)
axes[3].set_xlim(0,maxval)


plt.legend()
plt.savefig('histcomparison4.jpg',dpi=1200)
plt.close()

# plotting the boxplots
plt.ylabel('AUC')
figs, axes =plt.subplots(nrows=1, ncols=4, sharey=True, squeeze=True,figsize=(12,8))
figs.suptitle('The boxplot comparison of normalized weighted-degree distributions for module-specific and common genes\n',fontsize=15)
# setting custom x axis labels
plt.margins(0.01)


axes[0].boxplot(deg1, showmeans=True)
axes[0].set_xlabel(x1, fontsize=10)
axes[0].set_ylabel('Normalized weighted-degree', fontsize=10)
plt.setp(axes[0].get_xticklabels(), visible=False)

#axes[0].set_title('Default')
#axes[0].boxplot(y, showmeans=True)
# plt.subplot(122)
axes[1].boxplot(deg2, showmeans=True)
axes[1].set_xlabel(x2, fontsize=10)
plt.setp(axes[1].get_xticklabels(), visible=False)

axes[2].boxplot(deg3, showmeans=True)
axes[2].set_xlabel(x3, fontsize=10)
plt.setp(axes[2].get_xticklabels(), visible=False)

axes[3].boxplot(deg4, showmeans=True)
axes[3].set_xlabel(x4, fontsize=10)
plt.setp(axes[3].get_xticklabels(), visible=False)

#plt.legend()
plt.savefig('boxplotcomparison4.jpg',dpi=1200)
plt.close()