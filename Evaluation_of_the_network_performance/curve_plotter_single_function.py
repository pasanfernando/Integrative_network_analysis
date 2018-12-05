import matplotlib.pyplot as plt
from operator import itemgetter
from collections import OrderedDict
import numpy as np

plt.switch_backend('agg')

# plotting the ROC curve or precision recall curve
# specify the curve type correctly

def ROCplotter(rocdic,func_name, curve_type):

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
    plt.title(title+'for'+func_name)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.plot(x, y)
    plt.margins(0)
    #plt.show()
    # saving the figure
    plt.savefig('single_function/'+func_name+title+'.pdf')
    # closing the figure
    plt.close()
    # calculating the area under the curve
    auc = np.trapz(y, x)
    print 'area under the curve is:',auc

    return

def metric_plotter(hash1,metric_name,func_name):
    x=[]
    y=[]

    hash1 = OrderedDict(sorted(hash1.items(), key=lambda t: t[0]))

    for i in hash1:
        x.append(i)
        y.append(hash1[i])

    #print x
    #print y
    plt.plot(x, y)
    plt.title(metric_name+' vs threshold plot for'+func_name)
    plt.xlabel('Threshold')
    plt.ylabel(metric_name)
    plt.margins(0)
    #plt.show()
    plt.savefig('single_function/'+func_name+metric_name+'plot.pdf')
    plt.close()

    return

# method to plot all ROC curves or precision-recall curves in the same plot

def general_ROCplotter(genrocdic, curve_type):

    if curve_type == 'roc':
        title = 'ROC curve'
        xlab = 'False positive rate'
        ylab = 'True positive rate'

    elif curve_type == 'precision_recall':
        title = 'Precision-recall curve'
        xlab = 'Recall'
        ylab = 'Precision'


    plt.title(title+'for all functions')
    plt.xlabel(xlab)
    plt.ylabel(ylab)

    for i in genrocdic:

        x=[]
        y=[]
        rocdic = genrocdic[i]
        rocdic = OrderedDict(sorted(rocdic.items(), key=lambda t: t[1]))
        for i in rocdic:
            xy = rocdic[i]
            x.append(xy[0])
            y.append(xy[1])
        #print x
        #print y

        plt.plot(x, y)
        #plt.show()

        # calculating the area under the curve
        auc = np.trapz(y, x)
        print 'area under the curve is:',auc
    plt.margins(0)
        # saving the figure
    plt.savefig('gen_'+title+'.pdf')
    # closing the figure
    plt.close()

    return