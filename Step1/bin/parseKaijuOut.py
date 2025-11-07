#!/usr/bin/env python 
__author__ = 'luca.cozzuto@crg.eu'
# -*- coding utf-8 -*-

#MODULES
import sys
import re
import optparse
import collections
import pprint
import copy 
import subprocess
import textwrap 
import gzip 
import os 
import math 
from collections import Counter


#BODY FUNTIONS
def options_arg():
	usage = "usage: %prog -n <sample id>  -g <genome size info>  -i <kaiju output table> -o <output table> -s <stat file>"
	parser = optparse.OptionParser(usage=usage)
	parser.add_option('-i', '--input', help='Input table produced by kaiju ', dest="input" )
	parser.add_option('-n', '--name', help='Sample id name', dest="name" )
	parser.add_option('-o', '--output',help='ouput table file', dest="wotus" )
	parser.add_option('-s', '--stats',help='ouput stats file', dest="stats" )
	parser.add_option('-g', '--genome_stats', help='genome info', dest="genome", default=False )
	parser.add_option('-p', '--profile_type', help='taxonomic or functional', dest="prof", default="Taxonomic")
	(opts,args) = parser.parse_args()
	if opts.input and opts.wotus and opts.stats and opts.name:pass
	else: parser.print_help()
	return (opts)
def __main__ ():
	genome_sizes = False
	if (opts.genome != False ):
		genome_sizes = parseGenome(opts.genome)
	
	parseCountFile(opts.input, opts.wotus, opts.stats, opts.name, opts.prof, genome_sizes)
		

#AUXILIARY MODULES
def parseGenome(genomefile):
	pp = pprint.PrettyPrinter(indent=4)
	out = {}
	infile = open(genomefile, 'r')
	for line in infile:
		line = line.rstrip()
		fields = line.split()
		id = fields[0]
		gsize = fields[4]
		out[id] = gsize.replace(",", "" )	
	return out


def parseCountFile(file, ofile, sfile, name, prof, gsizes=False):
	pp = pprint.PrettyPrinter(indent=4)
	unique = {}
	multi = {}
	alns = []
	tot = 0
	aln = 0 
	# Prepare files
	infile = open(file, 'r')
	if (file.endswith('.gz')):
		infile = gzip.open(file, 'rt')
	for line in infile:
		line = line.rstrip()
		tot = tot + 1
		if (prof=="Functional"):
		  # Reads are mapped
		  	aln = aln + 1
		  	readname, genomeid_list  = line.split('\t')
		  	genomes = genomeid_list.rstrip(",").split(",")
		  	# unique genomes
		  	mygenomes = set(genomes)
		  	genomes = list(mygenomes)
	  
		  	# Reads are mapped in a single genome (unique)
		  	if (len(genomes)==1):
		  		if (genomes[0] not in unique):
		  			unique[genomes[0]] = 1
		  		else:
		  			unique[genomes[0]] = unique[genomes[0]] + 1
		  	# Reads are mapped in more genomes (multi)
		  	elif (len(genomes)> 1):
		  		alns.append(genomes)
		else:
		  if (line[:1]=="C"):
		  	# Reads are mapped
		  	aln = aln + 1
		  	mapped, readname, score, genomeid_list, translation_list  = line.split('\t')
		  	genomes = genomeid_list.rstrip(",").split(",")
		  	# unique genomes
		  	mygenomes = set(genomes)
		  	genomes = list(mygenomes)
	  
		  	# Reads are mapped in a single genome (unique)
		  	if (len(genomes)==1):
		  		if (genomes[0] not in unique):
		  			unique[genomes[0]] = 1
		  		else:
		  			unique[genomes[0]] = unique[genomes[0]] + 1
		  	# Reads are mapped in more genomes (multi)
		  	elif (len(genomes)> 1):
		  		alns.append(genomes)
	infile.close()

	# Scan multi alignments for removing genomes not in unique, then re-assign the new unique and remove from multi alignments
	clean_alns = alns.copy()
	for i, aln_genomes in enumerate(alns):
		new_genomes = aln_genomes.copy()
		for o, aln_genome in enumerate(aln_genomes):
			if (aln_genome not in unique):
				new_genomes.remove(aln_genome)
		if (len(new_genomes) == 1):
			new_uni_genome = new_genomes[0]
			unique[new_uni_genome] = unique[new_uni_genome] + 1
			clean_alns[i] = {} 
		else:
			clean_alns[i] = aln_genomes
	
	# Scan multi alignments for weighting them
	for aln_genomes in clean_alns:
		tot_uni = 0	
		weigth = 0	
		for aln_genome in aln_genomes:
			if (aln_genome in unique):
				tot_uni+=unique[aln_genome]
		for aln_genome in aln_genomes:
			if (aln_genome in unique):
				weigth = unique[aln_genome]/tot_uni
				if (aln_genome not in multi):
					multi[aln_genome] = weigth
				else:
					multi[aln_genome] = multi[aln_genome] + weigth
			
	weigthed_totals = {}
	for genome in unique:
		if (genome in multi):
			weigthed_multi = multi[genome]
			weigthed_totals[genome] = weigthed_multi + unique[genome]
		else:
			weigthed_totals[genome] = unique[genome]

	f = open(ofile,'w')
	inter = open(ofile + ".tmp",'w')
	header = "genome_id\tweighted counts"
	
	if(gsizes!=False):
		header = header + "\tnormalized counts"
		
	f.write(header + "\n")	
	wtot_all = 0
	for genome in weigthed_totals:
		wtotal = round(float(weigthed_totals[genome]), 2)
		wtot_all = wtot_all + wtotal
		uniq = unique[genome]
		multic = 0
		if(genome in multi):
			multic = multi[genome]
		inter.write(genome + "\t" + str(uniq) + "\t" + str(multic) + "\n")
		if(gsizes!=False):
			gsize = int(gsizes[genome])/1000000
			ntotal = wtotal/gsize
			f.write(genome + "\t" + str(wtotal) + "\t" + str(ntotal) + "\n")
		else:
			f.write(genome + "\t" + str(wtotal) + "\n")
	f.close()

	sf = open(sfile,'w')
	sf.write("sample\treads\taligned\n")	
	sf.write(name + "\t" + str(tot) + "\t" + str(wtot_all) + "\n")	
	sf.close()
	
	
	return


#Calling
opts = options_arg()
if opts.input and opts.wotus and opts.stats and opts.name:
	__main__()

