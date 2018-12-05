##this code uses string alias file to generate string iD to gene names mapping file
# it converts gene symbols to lower case for easier further analysis
# please specify species code during the code otherwise code wont work
# Please change the species id accordingly
in1 = raw_input('Enter input file: ')

out = open('namereplaced.txt', 'wb+')
mis = open('missingids.txt', 'wb+')

f= open(in1, 'r')
human =[]
humangn = {}
blproteins = []
uniproteins =[]
ensproteins = []
unisymbol = {}
ensemblesymbol ={}
blastkegg = {}

allprotein =[]
def printlength(dic):
    print 'original proteins number with names', len(dic)
    print 'protein number with names without duplicates', len(set(dic))

for line  in f:
        # change the species id accordingly
        if line.startswith('10090'):

            line= line.strip()
            a = line.split('\t')
            allprotein.append(a[0])

            #print a
            #activate the following line for humans
            if 'BLAST_KEGG_NAME BLAST_UniProt_GN BioMart_HUGO' in line:
                human.append(a[0])
                humangn[a[0]] = a[1]

            # activate the following line for mouse
            if 'Ensembl_EntrezGene ' in line:
                ensproteins.append(a[0])
                if a[0] not in ensemblesymbol:
                    ensemblesymbol[a[0]] = a[1]
            # excluding uniprot ids and only geting the ensemble uniprot name
            if a[2] == 'Ensembl_UniProt' or 'Ensembl_UniProt ' in line and 'Ensembl_UniProt_ID' not in line:
                uniproteins.append(a[0])
                unisymbol[a[0]] = a[1]
            if a[2] == 'BLAST_UniProt_GN'or'BLAST_UniProt_GN ' in line or 'BLAST_KEGG_NAME' in line:
                blproteins.append(a[0])
                blastkegg[a[0]] = a[1]
            #print a

# print 'original proteins number with names',len(proteins)
# print 'protein number with names without duplicates',len(set(proteins))
# printlength(uniproteins)
# printlength(ensproteins)
#printlength(blproteins)
# print len(unisymbol)
# print len(allprotein)

allprotein = set(allprotein)
print len(allprotein)
genes = {}
out.write("string_id name\n")

#########activate except humans##############################
for i in allprotein:
    if i in ensemblesymbol and i in blastkegg:
        if ensemblesymbol[i] == blastkegg[i]:
            if not ensemblesymbol[i].startswith('LOC'):
                genes[i] = ensemblesymbol[i]
        else:
            if not ensemblesymbol[i].startswith('LOC'):
                #out.write('%s %s\n'%(i,ensemblesymbol[i]))
                genes[i] = ensemblesymbol[i]
            else:
                if i in blastkegg:
                    genes[i] = blastkegg[i].lower()

    elif i in blastkegg:
        genes[i] = blastkegg[i].lower()
    elif i in unisymbol:
        genes[i] = unisymbol[i].lower()

    else:
        mis.write('%s\n'%(i))

print 'final mapped gene names count:',len(genes)

for i in genes:
    out.write('%s %s\n' % (i, genes[i].lower()))
###################################################################
#activate the following line for humans
# for i in allprotein:
#     if i in humangn:
#         out.write('%s %s\n' % (i, humangn[i].lower()))
#     else:
#         mis.write('%s\n' % (i))
#
# print 'final mapped gene names count:',len(humangn)