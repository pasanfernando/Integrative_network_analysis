
'''Perform network reconciliation using anatomy profiles and raw network files (name replaced)
Inputs
Get the anatomy profiles for the genes
StringIDs to name mapping file (genenamer.py)
StringIDs to ensemble id mapping file (string_ensemblemapper.py)
Name replaced raw (without cutoff) network file
'''

# prompting for the anatomical profile
#aprofile =raw_input('Input the name of the anatomical profile file:')

#prompting for the network file
#net =raw_input('Input the name of the network file:')

# Opening the anatomical profile file and storing the necessary information


# Dictionary to store anatomy_gene_name : ensembleid
aname_ens ={}

# a list to store gene names from anatomical profiles
angenelist =[]

anprofile = open('mouse_finalfull_profiles.txt', 'r')
for line in anprofile:
    if line != '\n' and 'ensemble_id' not in line:
        a = line.strip().split('\t')
        # removing genes without gene name
        if a[1] != 'not_found':
            aname_ens[a[1]] =a[2]
            angenelist.append(a[1])

        # detecting especial cases where gene name is not found but there is a ensemble id
        else:
            if a[2]!='not_found':
                print 'special gene ensemble without gene name:',a[0],a[1]

anprofile.close()
#print aname_ens

#finding duplicates in the anatomy profile gene list
seen = set()
uniq = []
for x in angenelist:
    if x not in seen:
        seen.add(x)
    else:
        uniq.append(x)
print 'dupliacted gene name stats are given below'
print uniq
print len(uniq)
#print len(angenelist)

# reading the string id:name file
stringname_id ={}

stringid_name ={}

strname = open('mousenames10.5.txt','r')

for line in strname:
    if line != '\n' and 'string_id' not in line:
        a = line.strip().split(' ')
        stringname_id[a[1]] =a[0].lower()
        stringid_name[a[0].lower()] = a[1]

# reading the ensembleid: string id file
stringens_id = {}

strname.close()

ensid = open('mouseensemblemappingstostring.txt', 'r')

for line in ensid:
    if line != '\n' and 'string_id' not in line:
        a = line.strip().split(' ')
        stringens_id[a[1]] = a[0].lower()

# Opening the gene network and storing the genes in the network in a list
ensid.close()

in1 = open('mouse10.5namenetwork.txt', 'r')

genelist =[]
for line in in1:
    if line != '\n':
        #print line
        if 'combined_score' not in line:
            line = line.strip('\n')
            a = line.split()
            #print a[0],a[1],a[2]
            gene1 = a[0].lower()
            # appending the gene 1 into gene list
            genelist.append(gene1)
            #print gene1
            gene2 = a[1].lower()
            # appending the gene 2 into gene list
            genelist.append(gene2)

networkgeneset = set(genelist)

in1.close()

# the number of overlapping genes in the network file
matchedgenes = set(angenelist)& networkgeneset

# anatomy profile genes that didnot match with the network
firstmismatch = set(angenelist)- networkgeneset

# how many of the mismatched genes are in the stringname mapping file
name_mapping_matches = firstmismatch & set(stringname_id.keys())

# A dictionary to store second round matches: anatomy gene name: networkname or id
secondmatches = {}

# a dictionary to reverse the above dictionary
secondmatchesreversed = {}

# how many of the mismatched genes have ensemble ids matching in the string network

# false synonyms in the second round
secondrumbles ={}
for i in firstmismatch:
    # checking their ensemble id in the stringdb
    if aname_ens[i] in stringens_id:
        #print i
        if stringens_id[aname_ens[i]] in stringid_name:
            #print stringid_name[stringens_id[aname_ens[i]]]
            if stringid_name[stringens_id[aname_ens[i]]] not in matchedgenes:
                secondmatches[i]=stringid_name[stringens_id[aname_ens[i]]]
                secondmatchesreversed[stringid_name[stringens_id[aname_ens[i]]]] = i

            else:
                secondrumbles[i]=stringid_name[stringens_id[aname_ens[i]]]

        else:
            #print stringens_id[aname_ens[i]]
            if stringens_id[aname_ens[i]] not in matchedgenes:
                secondmatches[i]=stringens_id[aname_ens[i]]
                secondmatchesreversed[stringens_id[aname_ens[i]]] =i
            else:
                secondrumbles[i]=stringens_id[aname_ens[i]]



# print 'number of genes in the network:',len(networkgeneset)
# print 'number of genes in the anatomy profile', len(aname_ens)

# out of the second matches how many are in the network
secondnetoverlap = set(secondmatches.values()) & set(networkgeneset)

# how many second matches are not in the network
secondmissing = set(secondmatches.values()) - set(networkgeneset)

# converting the secondnetoverlaps to anatomy gene names
anatomysecondnetoverlaps =[]
for i in secondnetoverlap:
    anatomysecondnetoverlaps.append(secondmatchesreversed[i])

# converting the secondmissing to anatomy gene names
anatomysecondmissing=[]
for i in secondmissing:
    anatomysecondmissing.append(secondmatchesreversed[i])

# Mismatches after second round
secondfinalnetmatches = set(anatomysecondnetoverlaps) | matchedgenes

# final mismatches
secondfinalmismatches = set(angenelist) - secondfinalnetmatches
# Third round of mapping
# in this method, the remaning mismatched are searched directly in the alias file

# Opening the alisa file
alias = open('10090.protein.aliases.v10.5.txt', 'r')

# dictionary to store alias: string id
alias_stringid = {}
for line in alias:
    # change the species id accordingly
    if line != '\n' and 'string_protein_id' not in line:
        line = line.strip()
        a = line.split('\t')
        if len(a) == 3:
            alias_stringid[a[1].lower()] = a[0].lower()
# database matches in the third round
thirddbmatches = []

# third round of false synonyms
thirdrumbles={}
# network matches in the third round anatomyname: string name
thirdnetmatches = {}
thirdmatchesreversed = {}
for i in secondfinalmismatches:
    if i in alias_stringid:
        thirddbmatches.append(i)
        # print i,alias_stringid[i]
        # trying to search the name for the string id in the network
        if alias_stringid[i] in stringid_name:
            if stringid_name[alias_stringid[i]] in networkgeneset:
                if stringid_name[alias_stringid[i]] not in matchedgenes and  stringid_name[alias_stringid[i]] not in anatomysecondnetoverlaps and stringid_name[alias_stringid[i]] not in secondmatchesreversed and stringid_name[alias_stringid[i]] not in thirdmatchesreversed:
                    thirdnetmatches[i] = stringid_name[alias_stringid[i]]
                    thirdmatchesreversed[stringid_name[alias_stringid[i]]] = i
                else:
                    thirdrumbles[i] = stringid_name[alias_stringid[i]]

        elif alias_stringid[i] in networkgeneset:
            if alias_stringid[i] not in matchedgenes and alias_stringid[i] not in anatomysecondnetoverlaps and alias_stringid[i] not in secondmatchesreversed and alias_stringid[i] not in thirdmatchesreversed:
                thirdnetmatches[i] = alias_stringid[i]
                thirdmatchesreversed[alias_stringid[i]] = i

            else:
                thirdrumbles[i] = alias_stringid[i]

# The genes that are in string db but not in the network
thirdmissing = set(thirddbmatches) - set(thirdnetmatches.keys())
# final matches to the network



finalnetmatches = set(anatomysecondnetoverlaps) | matchedgenes | set(thirdnetmatches.keys())

# final mismatches
finalmismatches = set(angenelist) - finalnetmatches

# final anatomyprofile genes that are in the database but not in the network

finaldatabaseonlymatches = name_mapping_matches | set(anatomysecondmissing) |thirdmissing






# Re modifying the anatomy profile file to only include reconciled genes with string network gene names
# correctly written gene number in profile
anmodgenenum=[]
anprofile = open('mouse_finalfull_profiles.txt', 'r')

# opening the modified file
anprofilemodified = open('mouse_anatomy_profilesmodified.txt', 'wb+')
for line in anprofile:
    if line != '\n' and 'ensemble_id' not in line:
        a = line.strip().split('\t')
        # if the gene is in the initial match list (direct name matches)


        if a[1] in matchedgenes:
            # replacing the gene name spaces with underscores
            a[1] = a[1].replace(' ', '_')
            anprofilemodified.write('%s\t%s\t%s\t%s\n'%(a[0],a[1],a[2],a[3]))
            anmodgenenum.append(a[1])

        # else if it is in the second matched list using the ensemble id, replace tha name with string name or id
        elif a[1] in anatomysecondnetoverlaps:
            a[1] = a[1].replace(' ', '_')
            anprofilemodified.write('%s\t%s\t%s\t%s\n' % (a[0], a[1], a[2], a[3]))
            anmodgenenum.append(a[1])

        elif a[1] in thirdnetmatches:


            a[1] = a[1].replace(' ', '_')
            anprofilemodified.write('%s\t%s\t%s\t%s\n' % (a[0], a[1], a[2], a[3]))
            anmodgenenum.append(a[1])


    else:
        anprofilemodified.write(line)

anprofile.close()

anprofilemodified.close()


# Because string names are the outdated ones, replacing network names
in1 = open('mouse10.5namenetwork.txt', 'r')

# Opening the modified network
in2 = open('mouse10.5namenetworkmodified.txt', 'wb+')


for line in in1:
    if line != '\n':
        #print line
        if 'combined_score' not in line:
            line = line.strip('\n')
            a = line.split()
            #print a[0],a[1],a[2]
            gene1 = a[0].lower()
            if gene1 in secondmatchesreversed:
                gene1 = secondmatchesreversed[gene1]
                gene1 = gene1.replace(' ', '_')

            if gene1 in thirdmatchesreversed:
                gene1 = thirdmatchesreversed[gene1]
                gene1 = gene1.replace(' ', '_')

            #print gene1
            gene2 = a[1].lower()
            if gene2 in secondmatchesreversed:
                gene2 = secondmatchesreversed[gene2]
                gene2 = gene2.replace(' ', '_')

            if gene2 in thirdmatchesreversed:
                gene2 = thirdmatchesreversed[gene2]
                gene2 = gene2.replace(' ', '_')

            in2.write('%s %s %s\n'%(gene1,gene2,a[2]))

        else:
            in2.write(line)



in1.close()



statfile = open('statfile.txt','wb+')

statfile.write('number of genes in the network: %s\n'%(len(networkgeneset)))
statfile.write('number of genes in the anatomy profile: %s\n'%(len(aname_ens)))
statfile.write('number of genes with string name mappings: %s\n'%(len(stringname_id)))
statfile.write('number of genes with string ensemble mappings: %s\n'%(len(stringens_id)))
statfile.write('\n')
statfile.write('number of anatomy profile genes that matched with the network: %s\n'%(len(matchedgenes)))
statfile.write('number of anatomy profile genes initially mismatched with the network: %s\n'%(len(firstmismatch)))
statfile.write('of the initial mismatched how many are in the stringname mapping file (string database, but not in the network): %s\n'%(len(name_mapping_matches)))
statfile.write('number of matches using ensemble ids (second matches): %s\n'%(len(secondmatches)))
statfile.write('out of second matches how many are in the network %s\n'%(len(secondnetoverlap)))
statfile.write('out of second matches how many are not in the network (these are in the database)%s\n'%(len(secondmissing)))
statfile.write('\n')
statfile.write('Third round database matches: %s\n'%(len(thirddbmatches)))
statfile.write('out of second matches how many are in the network %s\n'%(len(thirdnetmatches)))
statfile.write('out of third database matches how many are not in the network: %s\n'%(len(thirdmissing)))
statfile.write('\n')
statfile.write('number of final matches to the network: %s\n'%(len(finalnetmatches)))
statfile.write('number of final genes in moified anatomy profile file: %s\n'%(len(anmodgenenum)))
statfile.write('number of final mismatches to the network: %s\n'%(len(finalmismatches)))
statfile.write('number of final anatomy profile genes that are in the string database but not in the network: %s\n'%(len(finaldatabaseonlymatches)))
statfile.write('\n')
statfile.write('*************** Gene lists ********************************\n')
statfile.write('number of matches using ensemble ids (second matches): %s\n'%(len(secondmatches)))
for i in secondmatches:
    statfile.write('%s\t%s\n'%(i,secondmatches[i]))
statfile.write('\n')
statfile.write('third round matches: %s\n'%(len(thirdnetmatches)))
for i in thirdnetmatches:
    statfile.write('%s\t%s\n'%(i,thirdnetmatches[i]))
statfile.write('\n')
statfile.write('number of final anatomy profile genes that are in the string database but not in the network: %s\n'%(len(finaldatabaseonlymatches)))
statfile.write('number of  anatomy profile genes that are in the string database but not in the network in first round (name only): %s\n'%(len(name_mapping_matches)))
for i in name_mapping_matches:
    statfile.write('%s\n'%(i))
statfile.write('number of  anatomy profile genes that are in the string database but not in the network in second round (ensemble matches): %s\n'%(len(secondmissing)))
for i in secondmissing:
    statfile.write('%s\t%s\n'%(secondmatchesreversed[i],i))
statfile.write('\n')
statfile.write('number of  anatomy profile genes that are in the string database but not in the network in third round: %s\n'%(len(thirdmissing)))
for i in thirdmissing:
    statfile.write('%s\t%s\n'%(alias_stringid[i],i))
statfile.write('\n')
statfile.write('number of final mismatches to the network: %s\n'%(len(finalmismatches)))
unmgenes =[]
for i in finalmismatches:

    statfile.write('%s\n'%(i))
    if 'unm_' in i:
        unmgenes.append(i)
statfile.write('number of unnamed mismatches (unm_): %s\n'%(len(unmgenes)))
statfile.write('\n')
statfile.write('####################Pseudo synonyms#######################\n')
statfile.write('These are the genes that matched during second and third round, but matched to a original gene in the profile\n')
statfile.write('In other words they are matching to a old name of anatomy profile gene that matched with the PPI network in a previous round\n')
statfile.write('These are not true synonyms\n')
statfile.write('During second round %s\n'%(len(secondrumbles)))
for i in secondrumbles:
    statfile.write('%s\t%s\n'%(secondrumbles[i],i))
statfile.write('\n')
statfile.write('During third round %s\n'%(len(thirdrumbles)))
for i in thirdrumbles:
    statfile.write('%s\t%s\n'%(thirdrumbles[i],i))
statfile.write('\n')
statfile.close()