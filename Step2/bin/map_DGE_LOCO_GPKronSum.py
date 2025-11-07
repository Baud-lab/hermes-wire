#!/usr/bin/env python

#-------------------------#
#  A. IMPORTING MODULES	  #
#-------------------------#	 

import sys
import os
import h5py
import numpy as np
import argparse as argp 
#sys.path.insert(0,'/users/abaud/fmorillo/5_heritability_tests/') # path to dir with social_data_wMT.py and dirIndirVD_wMT
from classes.social_data_wMT import SocialData
import pdb
import copy
import gc
import h5py
import re
import csv
import time
import pandas as pd

from limix_core.gp import GP2KronSum
from limix_core.covar import FreeFormCov
from limix_lmm import LMM

#-------------------------#
#	B. DEFINING OPTIONS	  #
#-------------------------#	 
parser = argp.ArgumentParser(description ='Analysis of phenotypic variance with univariate or bivariate model')

## POSITIONAL ARGUMENTS - NEEDED and IN ORDER: input path / pheno version / covs version / cage version / GRM version
parser.add_argument('in', help = 'full path to input hdf5 file')
parser.add_argument('phenos_version', help= 'type of PHENOTYPES to use, i.e. name of subgroup of `phenotypes` in .h5, e.g. simulations')
parser.add_argument('covs_version', help= 'type of COVARIATES to use, i.e. name of subgroup of `covariates` in .h5, e.g. from_file')
parser.add_argument('cage_version', help= 'type of CAGE to use, i.e. name of subgroup of `cages` in .h5, e.g. real')
parser.add_argument('GRM_version', help= 'kinship type, e.g. prunes_dosage')
#parser.add_argument('chr', help= 'indicates if there is a kinship table per chromosome, e.g. yes')

## OPTIONAL ARGUMENTS
parser.add_argument('-o','--out', help = 'path to folder in which create output folders and files [default="."]', default=".")
parser.add_argument('-p','--pheno', type=int, help = "Number of phenotype's column to analyse [default=1]", default=1)
parser.add_argument('-e','--effects',help='effects - DGE,cageEffect - to be included [default= None to all]', default=None) # FOR NOW, IGE not included, eventually will have to change help including it (DGE,IGE,IEE,cageEffect). Remove "cageEffect" to fit a LMM with random DGE and random cage effects
parser.add_argument('-s','--subset', help='to handle optional subset on individuals [default=None]', default=None)
#parser.add_argument('-d','--indir', help='input directory [default="."]', default=".")
#parser.add_argument('-c','--combins_path',nargs="?", help = 'path to file with pairs of phenotypes to compare (.csv), needed for bivariate') # NOT USED IN UNIVARIATE
#parser.add_argument('-m','--model', help='type of model: univariate (uni) or bivariate (bi) [default=uni]', default="uni") # FOR NOW MODEL IS HARD CODED AS UNIVARIATE

args = vars(parser.parse_args())
start_time = time.time()

if __name__=='__main__':
	
	#---------------------------------------------------------#
	#	 1. GETTING OPTIONS FOR INPUT FILE = ..h5 file		  #
	#---------------------------------------------------------#
	in_file=args['in']
	phenos_version=args['phenos_version'] # e.g. phenos_version = 'bacterial_taxa'
	covs_version=args['covs_version'] 
	cage_version=args['cage_version']
	GRM_version = args['GRM_version'] # e.g. GRM_version = 'pruned_dosages'
	#chr = args['chr']

	
	#---------------------------------------------------------#
	#				  2. GETTING MODEL						  #
	#---------------------------------------------------------#
	
	## FOR NOW, model hard coded to univariate
	model = "uni" # Comment this line if implement bivariate
	#m = args['model']


	#---------------------------------------------------------#
	#	 3. PHENOTYPE'S COLUMN NUMBER TO PYTHON ANNOTATION	  #
	#---------------------------------------------------------#
	if model == "uni": 
		col = (args['pheno']-1)

	#---------------------------------------------------------#
	#			  4. GETTING EFFECTS TO INCLUDE				  #
	#---------------------------------------------------------#
 
	## Setting all to None by default
	DGE = None
	IGE = None
	IEE = None
	cageEffect = None
	
	effects = [DGE, IGE, IEE, cageEffect]
	if args['effects'] is not None:
		effects = args['effects'].split(",")
	
	if "DGE" in effects:
		DGE = "DGE"
		
	if "cageEffect" in effects:
		cageEffect = "cageEffect"
	
	## IGE FOR NOW NOT INCLUDED, will have to uncomment here if want to include it
	#if "IGE" in effects: 
	#	 IGE = "IGE"
	#	 IEE = "IEE"
	
	print("Effects included in model: ", DGE, cageEffect) # When include IGE, change here with ("Effects included in model: ", DGE, IGE, IEE, cageEffect)
	
	
	#-----------------------------------------------------------------#
	#  5. GETTING SUBSET - to handle optional subset on individuals	  #
	#-----------------------------------------------------------------#	 
	subset=args['subset']
	print("Subset: ",  subset)

	
	#-----------------------------------------------------------------#
	#  6. PARSING INPUT .h5 with SocialData() function		  #
	#  Also PARSING INPUT covar and genotypes         		  #
	#-----------------------------------------------------------------#	 

	chr=1 #chr below could be any - only used to get trait, final sampleID etc in covar file 

	## will now use code in social_data_wMT.py	
	data = SocialData(in_file, phenos_version, covs_version, cage_version, GRM_version, subset, chr)
	## Arguments' order in SocialData(self, in_file=None, phenos_version = None,covs_version=None, cage_version=None, GRM_version = None, subset = None, chr = None)
	# for univar
	if model == "uni":
		doto = data.get_data(col)
	# for bivar
	if model == "bi":
		doto = data.get_data(col,col_MT)
	print("Different covariates actually used")
	
	trait=doto['trait']
	print('trait is ' + trait)


	#-----------------------------------------------------------------#
	#    7. DEFINING OUTPUT DIR					  #
	#-----------------------------------------------------------------#	 
   
	## NB: if not existing, CREATE A NEW DIR in the one given as option
	## Try opening pvalues file early so that lmm doesnt run if file opening is going to fail...
	out_dir=args['out']
	#in_dir=args['indir']
	
	if out_dir is None:
		out_dir = os.getcwd()
	pvalues_file_dir = "".join([out_dir,"/",model,'variate/',phenos_version,"/pvalues_LOCO/",GRM_version,"_","_".join(filter(None, ['',subset, DGE, IGE, IEE, cageEffect]))])
	if not os.path.exists(pvalues_file_dir):
		os.makedirs(pvalues_file_dir, exist_ok=True)
	pvalues_file_name = "".join([pvalues_file_dir,'/',trait,'.h5'])

	covar_outfile_dir = "".join([out_dir,"/",model,'variate/',phenos_version,"/null_covars_LOCO/",GRM_version,"_","_".join(filter(None, ['',subset, DGE, IGE, IEE, cageEffect])),'/',trait])
	covar_outfile_name = "".join([covar_outfile_dir,"/",trait,"_chr",str(chr),".h5"])
	covar_outfile = h5py.File(covar_outfile_name,'r')
	saved = {}
	saved['sampleID'] = covar_outfile['sampleID'][:]
	saved['pheno'] = covar_outfile['pheno'][:]
	saved['covs'] = covar_outfile['covs'][:]
	saved['covar_mat'] = covar_outfile['covar_mat'][:]
	covar_outfile.close()

	# Import genotypes. direct and social geno already matched on rows and cols
	input_file = h5py.File(in_file,'r')
	geno  = input_file['direct'] #Change
	geno_matrix = geno['matrix'][:].T
	geno_sample_ID = geno['row_header']['sample_ID'][:]
	position = {
    "chr": np.asarray(geno['col_header']['chr'][:]),
    "pos": np.asarray(geno['col_header']['pos'][:]),
  }
	input_file.close()
	
	# Match genotypes with the rest
	Imatch = np.nonzero(saved['sampleID'][:,np.newaxis]==geno_sample_ID)
	print("Number of individuals in sampleID and with genotypes: " + str(len(Imatch[0])))
	saved['sampleID'] = saved['sampleID'][Imatch[0]]
	saved['pheno'] = saved['pheno'][Imatch[0]]
	saved['covs'] = saved['covs'][Imatch[0],:]
	geno_matrix = geno_matrix[Imatch[1],:]
	saved_geno_matrix = copy.copy(geno_matrix)
	
	pvalues_file = h5py.File(pvalues_file_name,'w')
	pvalues_file.create_dataset(name = 'chr',data = position['chr'])
	pvalues_file.create_dataset(name = 'pos',data = position['pos'])
	
	#-----------------------------------------------------------------#
	#						   8. MAPPING							  #
	#-----------------------------------------------------------------#	 

	# Will define one LMM per chr as covar changes for each cov
	for chr in range(1, 21):
		print('chr is ' + str(chr))
		
		covar_outfile_name = "".join([covar_outfile_dir,"/",trait,"_chr",str(chr),".h5"])
		covar_outfile = h5py.File(covar_outfile_name,'r')
		if "covar_mat" in covar_outfile.keys(): 			
			saved['covar_mat'] = covar_outfile['covar_mat'][:][Imatch[0],:][:,Imatch[0]]
			covar_outfile.close()
		else:
			print("Warning: mapping impossible for chr",chr,"as no null covar available'P values of -999 will be returned, to be understood as missing values")		
			nb = (position['chr'] == chr).sum()
			missing_pvalues = np.full(nb, -999)
			pvalues_file.create_dataset(name = "".join(['pvalues_chr',str(chr)]),data = missing_pvalues)
			continue
		
		#Notation from GP2KronSum:
		#N = number of samples
		#P = number of traits
		#Y = [N, P] phenotype matrix
		#F_i = sample fixed effect design for term i
		#A_i = trait fixed effect design for term i
		#B_i = effect sizes of fixed effect term i
		#Cg = column covariance matrix for signal term respectively
		#Cn = column covariance matrix for noise term respectively
		#R = row covariance matrix for signal term respectively
	
		K = saved['covar_mat']
		Cg = FreeFormCov(1)
		Cn = FreeFormCov(1)
		A = np.eye(1)
	 
		gp = GP2KronSum(Y=saved['pheno'], F=saved['covs'], A=A, Cg=Cg, Cn=Cn, R=K)	  
		gp.covar.Cr.setCovariance(0.5 * np.ones((1, 1)))
		gp.covar.Cn.setCovariance(0.000001 * np.ones((1, 1)))
		info_opt = gp.optimize(verbose=False)
	
		lmm = LMM(saved['pheno'], saved['covs'], gp.covar.solve)
	
		#print('additional noise VC is ')
		#print(Cn.K())
	
		geno_matrix = saved_geno_matrix[:,position['chr']==chr]
		lmm.process(geno_matrix)
		pvalues = lmm.getPv()
		#pdb.set_trace()
		#beta = lmm.getBetaSNP()
		#pvalues = pvalues[0,:]
		pvalues_file.create_dataset(name = "".join(['pvalues_chr',str(chr)]),data = pvalues)
		gc.collect()
	
print("Files saved to", pvalues_file_name)
print("My program took", time.time() - start_time, "sec to run")
	
