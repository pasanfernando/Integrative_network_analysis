#!/usr/bin/python

import re
import collections



p = open('uberonext.obo', 'r')
uberon = {}
name = {}
gene = collections.defaultdict(list)  # creating a multidict to store gene name with zfa id

out = open('uberon_names.txt', 'wb+')

out.write('gene_symbol\tzfin_id\tensemble_id\tuberon_ID\tname\n')




# this code only stores uberon ids
for line in p:

	if line.startswith('[Term]'):
		y = None
		m = None
		d = None

	if line.startswith('id: UBERON'):
		x = line.split()
		y = x[1]
	# print y
	elif 'xref: ZFA:' in line:
		z = line.split()
		m = z[1]
		# print m
		uberon[m] = y

	elif line.startswith('name:'):
		k = re.sub('name: ', '', line, 1)
		z = k.strip()
		d = z.replace(" ", "_")
		name[y] = d

	# when typedef approaches stop reading the file: typedef corresponds to relationships between terms
	if line.startswith('[Typedef]'):
		break

