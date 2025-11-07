#!/usr/bin/env python

#-------------------------#
#  A. IMPORTING MODULES	  #
#-------------------------#	 

import os
import sys
import argparse as argp 
from classes.social_data_wMT import SocialData	 
from classes.dirIndirVD_wSS import DirIndirVD 
import pdb
import h5py
import re
import gc
import csv
import time
import pandas as pd
import numpy as np

#-------------------------#
#	B. DEFINING OPTIONS	  #
#-------------------------#	 
parser = argp.ArgumentParser(description ='Analysis of phenotypic variance with univariate or bivariate model')

## POSITIONAL ARGUMENTS - NEEDED and IN ORDER: input path / pheno version / covs version / cage version / GRM version
parser.add_argument('in', help = 'path to input hdf5 file')
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
	#		  6. PARSING INPUT .h5 with SocialData() function		  #
	#-----------------------------------------------------------------#	 
	## will now use code in social_data_wMT.py 
	
	for chr in range(1,21) :

		print('chr is ' + str(chr))
		data = SocialData(in_file, phenos_version, covs_version, cage_version, GRM_version, subset, chr)
		## Arguments' order in SocialData(self, in_file=None, phenos_version = None,covs_version=None, cage_version=None, GRM_version = None, subset = None, chr = None)
		# for univar
		if model == "uni":
			doto = data.get_data(col)
		# for bivar
		if model == "bi":
			doto = data.get_data(col,col_MT)
	
		trait=doto['trait']
		print('trait is ' + trait)
		
	#-----------------------------------------------------------------#
	#					 7. DEFINING OUTPUT DIR						  #
	#-----------------------------------------------------------------#	 
   
	## NB: if not existing, CREATE A NEW DIR in the one given as option
		out_dir=args['out']
		if out_dir is None:
			out_dir = os.getcwd()
		if chr == 1:
			covar_outfile_dir = "".join([out_dir,"/",model,'variate/',phenos_version,"/null_covars_LOCO/",GRM_version,"_","_".join(filter(None, ['',subset, DGE, IGE, IEE, cageEffect])),'/',trait])
			if not os.path.exists(covar_outfile_dir):
				os.makedirs(covar_outfile_dir, exist_ok = True) # exist_ok=True should avoid problem in case it already exists
		covar_outfile_name="".join([covar_outfile_dir,'/',trait,'_chr',str(chr),'.h5'])
		if os.path.exists(covar_outfile_name):
			covar_outfile = h5py.File(covar_outfile_name,'r')
			if "covar_mat" in covar_outfile.keys():	#meaning the pipeline has already been run successfully	
				continue

		#if the pipeline has never yet been run or if it has not been run successfully, we'll start trying to run VD up to 5 times
		try_nb_VD = 0
		#sid = 255 # This with -p 380 (trait "g__NSJ-51_MI"); gives ValueError with chr 7 and LinAlgError with chr 14
		while try_nb_VD < 6:
			
	#-----------------------------------------------------------------#
	#				8. VARIANCE DECOMPOSITION ANALYSIS				  #
	#-----------------------------------------------------------------#	 
	
			try_nb_VD = try_nb_VD + 1
		
			## TODO: change here options to vary analysis 
			#pdb.set_trace()
			#print("seed =", sid); np.random.seed(sid); sid = sid+1
			vc = DirIndirVD(pheno = doto['pheno'], pheno_ID = doto['pheno_ID'], covs = doto['covs'], covs_ID = doto['covs_ID'], covariates_names = doto['covariates_names'], kinship_all = doto['kinship_full'], kinship_all_ID = doto['kinship_full_ID'], cage_all = doto['cage_full'], cage_all_ID = doto['cage_full_ID'], DGE = (DGE is not None), IGE = (IGE is not None), IEE = (IEE is not None), cageEffect = (cageEffect is not None), calc_ste=False, subset_IDs = doto['subset_IDs'], SimplifNonIdableEnvs = False, vc_init_type = None, vc_init = None)
			dirIndir_out = vc.getOutput()

			if "optim_error" in dirIndir_out.keys():
			    assert dirIndir_out['optim_error'] is not None, "optimization error is None, when expected to be something"
			    print("Warning: optimization was not successful, have still saved sampleID, pheno and covs")		
			else:
			    break
				
		#will save Infos whether optimization successful (after 5 tries) or not
		toSave_file = h5py.File(covar_outfile_name,'w')
		toSaveInfos = vc.getToSaveInfos()
		toSave_file.create_dataset(name = 'sampleID',data = toSaveInfos['sampleID'])
		toSave_file.create_dataset(name = 'pheno',data = toSaveInfos['pheno'])
		toSave_file.create_dataset(name = 'covs',data = toSaveInfos['covs'])
		#will save Covar if optimization was successful
		if "optim_error" in dirIndir_out.keys():
			#here need to return a warning to let the user know but nothing more as we don't want the pipeline to fail in case optim was not successful
			print("Warning: after 5 tries optimization was still not successful; have still saved sampleID, pheno and covs")				
		else:
			toSaveCovar = vc.getToSaveCovar()
			toSave_file.create_dataset(name = 'covar_mat',data = toSaveCovar['covar_mat'])

		toSave_file.close()
		gc.collect()

print("Files saved to", covar_outfile_dir)
print("My program took", time.time() - start_time, "sec to run")
