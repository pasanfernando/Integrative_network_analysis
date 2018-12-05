
from collections import defaultdict

'''generates a list of zebrafish ortholog genes for mouse specific genes'''

genelist =[]

sfile = open('forelimb_specific_zebragenes.txt','r')

for line in sfile:
    if line != '\n' and 'zebrafish' not in line:
        a =line.strip().split('\t')
        a[0] = a[0].strip('\"')
        if a[0]!= 'NA':
            #print a[0]
            b = a[0].split(':')
            #print b
            if ','in b[1]:
                c =b[1].split(',')
                for i in c:
                    print i
                    genelist.append(i)
            else:
                genelist.append(b[1])


            #print gene


sfile.close()
#print sdic

#print genelist

out = open('mouse_specific_zebrafishgenes.txt','wb+')
for i in genelist:
    out.write('%s\n'%(i))
out.close()

lostgenes =[]

lfile = open('pectoral_lost_genes.txt','r')

for line in lfile:
    if line != '\n':
        g = line.strip()
        lostgenes.append(g)



print lostgenes

comgenes = set(lostgenes) & set(genelist)
print 'common genes:',comgenes

out1 = open('statsmouse_specific_zebrafishgenesstats.txt','wb+')
out1.write('Number of mouse specific genes in zebrafish: %s\n'%len(genelist))

out1.close()