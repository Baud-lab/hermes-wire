#!/usr/bin/env python 
__author__ = 'luca.cozzuto@crg.eu'
# -*- coding utf-8 -*-

#MODULES
import sys
import re
import optparse
from collections import defaultdict
import pprint
import copy 
import os 


#BODY FUNTIONS
def options_arg():
	usage = "usage: %prog -f <file list > -o <output table> -n <norm output table (optional)> -k <keep header>"
	parser = optparse.OptionParser(usage=usage)
	parser.add_option('-f', '--files', help='Input file list', dest="input" )
	parser.add_option('-o', '--output',help='ouput table file', dest="wotus" )
	parser.add_option('-n', '--noutput', help='ouput normalised table file', dest="notus", default=False )
	parser.add_option('-k', '--keep_header', help='keep one header', action="store_true", default=False, dest="keep_header" )

	(opts,args) = parser.parse_args()
	if opts.input and opts.wotus:pass
	else: parser.print_help()
	return (opts)
def __main__ ():
	parseFiles(opts.input, opts.wotus, opts.notus, opts.keep_header)

#AUXILIAR MODULES
def parseFiles(filelist, ofile, nfile=False, kheader=False):
	pp = pprint.PrettyPrinter(indent=4)
	outres = defaultdict(dict)
	noutres = defaultdict(dict)

	header = ""
	ifiles = filelist.split(",")
 	# prepare files files
	for ifile in ifiles:
		nrow = 0
		infile = open(ifile, 'r')
		for line in infile:
			line = line.rstrip()
			nrow = nrow + 1
			if (kheader == True):
				if (nrow == 1):
					header = line
				else:
					fields = line.split('\t')
					id = fields[0]
					#pp.pprint(fields)
					outres[id][ifile] = fields[1]
					if (nfile!=False):
						noutres[id][ifile] = fields[2]
			else:
				fields = line.split('\t')
				id = fields[0]
				outres[id][ifile] = fields[1]
				if (nfile!=False):
					noutres[id][ifile] = fields[2]
		infile.close()
	
	outstr = "Genome_id" + "\t" + "\t".join(ifiles) + "\n" 
	
	for id, values in outres.items():
		vals = []
		for ifile in ifiles:
			val = "0"
			if(ifile in values):
				val = str(values[ifile])
			vals.append(val)
		outstr += id + "\t" + "\t".join(vals) + "\n"
			
	f = open(ofile,'w')
	f.write(outstr)
	f.close()
	
	if (nfile!=False):
		outstr = "Genome_id" + "\t" + "\t".join(ifiles) + "\n" 

		for id, nvalues in noutres.items():
			vals = []
			for ifile in ifiles:
				nval = "0"
				if(ifile in nvalues):
					nval = str(nvalues[ifile])
				vals.append(nval)
			outstr += id + "\t" + "\t".join(vals) + "\n"
			
		f = open(nfile,'w')
		f.write(outstr)
		f.close()
	
	return


#Calling
opts = options_arg()
if opts.input and opts.wotus:
	__main__()

