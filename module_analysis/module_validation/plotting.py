import numpy as np
import matplotlib.pyplot as plt

input1 = raw_input('enter the interaction file:')
input2 = raw_input('enter the species name:')
in1 = open(input1, 'r')

y =[]
ytot =0
ztot =0
z = []
for line in in1:
    if line != '\n':
        line = line.strip()
        if 'gene_name' not in line:
            a = line.split()
            y.append(int(a[1]))
            z.append(float(a[3]))
            ytot = ytot + int(a[1])
            ztot = ztot + float(a[3])
            #print a[1]

print ytot
print ztot

yrel =[]


for i in y:

    x = float(i)/float(ytot)

    yrel.append(x)



f = plt.figure(1)
plt.hist(y, color='c', histtype='bar', bins=50)
plt.title(input2+' candidate gene count distribution')
plt.xlabel(' gene count')
plt.ylabel('absolute frequency')
f.show()
#plt.hist(z, color='g', histtype='bar')

f.savefig(input2+'genecount.png')
g = plt.figure(2)
plt.hist(z, color='c', histtype='bar')
plt.title(input2+'percentage candidate gene count distribution')
plt.xlabel('gene count')
plt.ylabel('frequency')
g.show()

plt.show()
g.savefig(input2+'relativegenecount.png')

# for relative frequency distribution
print yrel
k = plt.figure(3)
plt.hist(y, color='b', histtype='bar', bins=50, normed=True)
plt.title(input2+'relative candidate gene count frequency distribution')
plt.xlabel(' gene count')
plt.ylabel('relative frequency')
k.show()

plt.show()
k.savefig(input2+'relativegenecount.png')