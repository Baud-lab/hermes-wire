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
from collections import Counter


#BODY FUNTIONS
def options_arg():
	usage = "usage: %prog -i <multiqc folder data> -o <output table>"
	parser = optparse.OptionParser(usage=usage)
	parser.add_option('-i', '--input', help='Data folder produced by multiqc ', dest="input" )
	parser.add_option('-o', '--output',help='ouput table file', dest="wotus" )
	(opts,args) = parser.parse_args()
	if opts.input and opts.wotus:pass
	else: parser.print_help()
	return (opts)
def __main__ ():
	parseMultiQCFiles(opts.input, opts.wotus)

#AUXILIAR MODULES
def parseMultiQCFiles(folder, ofile):
	pp = pprint.PrettyPrinter(indent=4)
	outresT = {}
	outres = {}
	trimmo = open(folder + "/multiqc_trimmomatic.txt", 'r')
	for trimmoline in trimmo:
		trimmoline = trimmoline.rstrip()
		SampleT, input_read_pairs, surviving, surviving_pct, forward_only_surviving, forward_only_surviving_pct, reverse_only_surviving, reverse_only_surviving_pct, dropped, dropped_pct = trimmoline.split('\t')
		if (input_read_pairs != "input_read_pairs"):
			if (SampleT not in outres):
				outresT[SampleT] = {}
			outresT[SampleT]['raw_reads'] = int(float(input_read_pairs))
			outresT[SampleT]['surviving'] = int(float(surviving))
			outresT[SampleT]['surviving_pct'] = surviving_pct
	trimmo.close()

	bowtie = open(folder + "/multiqc_bowtie2.txt", 'r')
	for bowtieline in bowtie:
		bowtieline = bowtieline.rstrip()
		bowtie_fields = bowtieline.split('\t')
		SampleB = bowtie_fields[0]
		for SampleT in outresT.keys():
			if SampleT.find(SampleB) > -1:
				pair_tot = int(bowtie_fields[2])
				ov_aln = float(bowtie_fields[10])
				outres[SampleB] = outresT[SampleT]
				outres[SampleB]['est_aligned'] = int(pair_tot*ov_aln/100)
				outres[SampleB]['overall_alignment_rate'] = ov_aln
	bowtie.close()

	kaiju = open(folder + "/multiqc_kaiju_alignment_statistics.txt", 'r')
	for kaijuline in kaiju:
		kaijuline = kaijuline.rstrip()
		kaiju_fields = kaijuline.split('\t')
		SampleK = kaiju_fields[0]
		if SampleK != "Sample":
			outres[SampleK]['available_reads'] = int(float(kaiju_fields[1])) 
			outres[SampleK]['available_reads_pct'] = int(100*(outres[SampleK]['available_reads']/outresT[SampleT]['surviving']))
			outres[SampleK]['aligned'] = int(float(kaiju_fields[2]))/2
			outres[SampleK]['mapping_ratio'] = int(100*(outres[SampleK]['aligned']/outres[SampleK]['available_reads']))

	kaiju.close()

	f = open(ofile,'w')
	f.write("Sample" + "\tRaw_Reads\tReads_After_Trimming\tSurviving_Rate_Trimming\tReads_Aligned_to_Host\tAlignment_Rate\tNon_Host_Reads\tSurviving_Rate_Alignment\tMapped_Reads\tRaw_Mapping_Ratio\n")	
	for sample in outres:
		f.write(sample + "\t" + 
			str(outres[sample]['raw_reads']) + "\t" + 
			str(outres[sample]['surviving']) + "\t" +
			str(outres[sample]['surviving_pct']) + "\t" +
			str(outres[sample]['est_aligned']) + "\t" +
			str(outres[sample]['overall_alignment_rate']) + "\t" +
			str(outres[sample]['available_reads']) + "\t" +
			str(outres[sample]['available_reads_pct']) + "\t" +
			str(outres[sample]['aligned']) + "\t" +
			str(outres[sample]['mapping_ratio']) + "\n"
		)
	f.close()
	return





#Calling
opts = options_arg()
if opts.input and opts.wotus:
	__main__()

