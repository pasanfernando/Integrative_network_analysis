
'''Compares enriched uberon terms from two different files
Please look into the file format before use
please define the correct names of the files in the middle of the code'''

from collections import defaultdict



#two multi dics to store term id: term name, pvalue
file1dic =defaultdict(list)
file2dic = defaultdict(list)

def dicpopulator(x,dic):
    # opening kellen's candidate genes
    file1 = open(x, 'r')

    for line in file1:
        if line!= '\n' and 'P-value' not in line:
            a = line.strip().split('\t')
            # only keeping the significant terms
            pval = float(a[1])
            if pval < 0.05:

                #print b
                dic[a[0]].append(a[2])
                dic[a[0]].append(pval)

    file1.close()
    return

# populating file1
dicpopulator('hindlimb_ori.txt',file1dic)
file1name = 'hindlimb_ori'

dicpopulator('hindlimb_pred.txt',file2dic)
file2name = 'hindlimb_pred'

print file1dic
print file2dic

# performing teh set comparisons

# common terms
commenterms = set(file1dic.keys()) & set(file2dic.keys())

# writing the common term file
comfile = open('common_terms.txt','wb+')
comfile.write('Uberon term\tUberon term name\tP-value for %s\tP-value for %s\n'%(file1name,file2name))
for i in commenterms:
    comfile.write('%s\t%s\t%s\t%s\n'%(i,file1dic[i][0],file1dic[i][1],file2dic[i][1]))

comfile.close()

# file1 specific terms
file1specific = set(file1dic.keys()) - set(file2dic.keys())


file1s = open(file1name+'specific_terms.txt','wb+')
file1s.write('Uberon term\tUberon term name\tP-value for %s\n'%(file1name))
for i in file1specific:
    file1s.write('%s\t%s\t%s\n'%(i,file1dic[i][0],file1dic[i][1]))

file1s.close()


# file2 specific terms
file2specific = set(file2dic.keys()) - set(file1dic.keys())


file2s = open(file2name+'specific_terms.txt','wb+')
file2s.write('Uberon term\tUberon term name\tP-value for %s\n'%(file2name))
for i in file2specific:
    file2s.write('%s\t%s\t%s\n'%(i,file2dic[i][0],file2dic[i][1]))

file2s.close()