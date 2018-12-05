import matplotlib.pyplot as plt
from operator import itemgetter
from collections import OrderedDict
import numpy as np

# plotting the ROC curve

def ROCplotter(rocdic,curve_type):

    if curve_type == 'roc':
        title = 'ROC curve'
        xlab = 'False positive rate'
        ylab = 'True positive rate'

    elif curve_type == 'precision_recall':
        title = 'Precision-recall curve'
        xlab = 'Recall'
        ylab = 'Precision'

    x=[]
    y=[]
    rocdic = OrderedDict(sorted(rocdic.items(), key=lambda t: t[1]))
    for i in rocdic:
        xy = rocdic[i]
        x.append(xy[0])
        y.append(xy[1])
    #print x
    #print y
    plt.title(title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.plot(x, y,'o-')
    #plt.show()
    # saving the figure
    plt.savefig(title+'.jpg')
    # closing the figure
    plt.close()
    # calculating the area under the curve
    auc = np.trapz(y, x)
    print 'area under the curve is:',auc

    return

def metric_plotter(hash1,metric_name):
    x=[]
    y=[]

    hash1 = OrderedDict(sorted(hash1.items(), key=lambda t: t[0]))

    for i in hash1:
        x.append(i)
        y.append(hash1[i])

    #print x
    #print y
    plt.plot(x, y)
    plt.title(metric_name+' vs threshold plot')
    plt.xlabel('Threshold')
    plt.ylabel(metric_name)
    #plt.show()
    plt.savefig(metric_name+'plot.jpg')
    plt.close()

    return