#!/usr/bin/env python

#-------------------------#
#  A. IMPORTING MODULES   #
#-------------------------#  
import os
import sys
import argparse as argp 
#sys.path.insert(0,'/users/abaud/fmorillo/Microbiome_Profiling/code/LIMIX/') # path to dir with social_data_wMT.py and dirIndirVD_wSS.py“
from classes.social_data_wMT import SocialData   
from classes.dirIndirVD_wSS import DirIndirVD 
import pdb
import h5py
import re
import csv
import time
import pandas as pd


#-------------------------#
#   B. DEFINING OPTIONS   #
#-------------------------#  
parser = argp.ArgumentParser(description ='Analysis of phenotypic variance with univariate or bivariate model')

## POSITIONAL ARGUMENTS - NEEDED and IN ORDER: input path / pheno version / covs version / cage version / dam version / GRM version
parser.add_argument('in', help = 'path to input hdf5 file')
parser.add_argument('phenos_version', help= 'type of PHENOTYPES to use, i.e. name of subgroup of `phenotypes` in .h5, e.g. simulations')
parser.add_argument('covs_version', help= 'type of COVARIATES to use, i.e. name of subgroup of `covariates` in .h5, e.g. from_file')
parser.add_argument('cage_version', help= 'type of CAGE to use, i.e. name of subgroup of `cages` in .h5, e.g. real')
parser.add_argument('dam_version', help= 'type of DAM to use, i.e. name of subgroup of `dams` in .h5, e.g. real')
parser.add_argument('GRM_version', help= 'kinship type, e.g. Andres_kinship')

## OPTIONAL ARGUMENTS
parser.add_argument('-o','--out', help = 'path to folder in which create output folders and files [default="."]', default=".")
parser.add_argument('-p','--pheno', type=int, help = "Number of phenotype's column to analyse [default=1]", default=1)
parser.add_argument('-e','--effects',help='effects - DGE,cageEffect,maternalEffect - to be included [default= None to all]', default=None) # FOR NOW, IGE not included, eventually will have to change help including it (DGE,IGE,IEE,cageEffect,maternalEffect)
parser.add_argument('-s','--subset', help='to handle optional subset on individuals [default=None]', default=None)
#parser.add_argument('-c','--combins_path',nargs="?", help = 'path to file with pairs of phenotypes to compare (.csv), needed for bivariate') # NOT USED IN UNIVARIATE
#parser.add_argument('-m','--model', help='type of model: univariate (uni) or bivariate (bi) [default=uni]', default="uni") # FOR NOW MODEL IS HARD CODED AS UNIVARIATE


args = vars(parser.parse_args())



start_time = time.time()

if __name__=='__main__':
    
    #---------------------------------------------------------#
    #    1. GETTING OPTIONS FOR INPUT FILE = ..h5 file        #
    #---------------------------------------------------------#
    in_file=args['in']
    phenos_version=args['phenos_version'] # e.g. phenos_version = 'bacterial_taxa'
    covs_version=args['covs_version'] 
    cage_version=args['cage_version']
    dam_version=args['dam_version']
    GRM_version = args['GRM_version'] # e.g. GRM_version = 'exp_DGA_7_Dec_2021'

    
    #---------------------------------------------------------#
    #                 2. GETTING MODEL                        #
    #---------------------------------------------------------#
    
    ## FOR NOW, model hard coded to univariate
    model = "uni" # Comment this line if implement bivariate
    #m = args['model']
    
    ## Will have to uncomment following for BIVARIATE
    #if m == "bivariate" or m == "bi":
    #    model = "bi"
    #    print("model is BIVARIATE: ", model)
    #
    #elif m == "univariate" or m == "uni":
    #    model = "uni"
    #    print("model is UNIVARIATE: ", model)
    #    
    #else:
    #    sys.exit("Need to know which model to use: " + "'"+ m +"' " + "is not valid" + "\n==> Use 'uni' / 'univariate' OR 'bi' / 'bivariate'")
    
    
    #---------------------------------------------------------#
    #    3. PHENOTYPE'S COLUMN NUMBER TO PYTHON ANNOTATION    #
    #---------------------------------------------------------#
    if model == "uni": 
        col = (args['pheno']-1)
    
    ## BIVARIATE MODEL - will have to uncomment 'elif' at line 76
    ##   few lines below go through a .csv file that has on each row PAIRS OF PHENOTYPES (INDXS as in the H5 phenotype matrix) saying which pairs of phenotypes you want to consider (doesnt have to be all combinations; 
    ##   for example for you it may look like
    ##       1 2
    ##       3 4
    ##       5 6 
    ##   if you put paired phenotypes next to each other in the phenotype matrix
    
    #elif model == "bi":
    #    assert args['combins_path'] is not None, 'Need combins file'
    #    # Reading combins file
    #    dfcombins = pd.read_csv(args['combins_path'], header=None, sep=' ')
    #    col, col_MT = dfcombins.iloc[args['pheno'] - 1].values - 1 # -1 here to go from base 1 in H5 file to base 0 in python
    

    #---------------------------------------------------------#
    #             4. GETTING EFFECTS TO INCLUDE               #
    #---------------------------------------------------------#
 
    ## Setting all to None by default
    DGE = None
    IGE = None
    IEE = None
    cageEffect = None
    maternalEffect = None
    
    effects = [DGE, IGE, IEE, cageEffect, maternalEffect]
    if args['effects'] is not None:
        effects = args['effects'].split(",")
    
    if "DGE" in effects:
        DGE = "DGE"
        
    if "cageEffect" in effects:
        cageEffect = "cageEffect"
    
    if "maternalEffect" in effects:
        maternalEffect = "maternalEffect"
    
    ## IGE FOR NOW NOT INCLUDED, will have to uncomment here if want to include it
    #if "IGE" in effects: 
    #    IGE = "IGE"
    #    IEE = "IEE"
    
    print("Effects included in model: ", DGE, cageEffect, maternalEffect) # When include IGE, change here with ("Effects included in model: ", DGE, IGE, IEE, cageEffect, maternalEffect)
    
    
    #-----------------------------------------------------------------#
    #  5. GETTING SUBSET - to handle optional subset on individuals   #
    #-----------------------------------------------------------------#  
    subset=args['subset']
    print("Subset: ",  subset)
    
    
    #-----------------------------------------------------------------#
    #         6. PARSING INPUT .h5 with SocialData() function         #
    #-----------------------------------------------------------------#  
    ## will now use code in social_data_wMT.py 
    
    ## Arguments' order in SocialData(self, in_file=None, phenos_version = None,covs_version=None, cage_version=None, dam_version=None, GRM_version = None, subset = None, chr = None)
    data = SocialData(in_file, phenos_version, covs_version, cage_version, dam_version, GRM_version, subset) 
    # for univar
    if model == "uni":
        doto = data.get_data(col)
    # for bivar
    if model == "bi":
        doto = data.get_data(col,col_MT)

    trait1=doto['trait']
    trait2=doto['trait_MT'] # This is for bivar, in univar will just be None
    print("traits are " + " and ".join(filter(None, [trait1,trait2])))

    
    #-----------------------------------------------------------------#
    #                    7. DEFINING OUTPUT DIR                       #
    #-----------------------------------------------------------------#  
   
    ## NB: if not existing, CREATE A NEW DIR in the one given as option
    out_dir=args['out']
    if out_dir is None:
        out_dir = os.getcwd()
    VD_outfile_dir = "".join([out_dir,"/",model,'variate/',phenos_version,"/",GRM_version,"_","_".join(filter(None, ['',subset, DGE, IGE, cageEffect, maternalEffect])),'/'])
    if not os.path.exists(VD_outfile_dir):
        os.makedirs(VD_outfile_dir, exist_ok=True) # exist_ok=True should avoid problem in case it already exists
    VD_outfile_name="".join([VD_outfile_dir,"_".join(filter(None, [trait1,trait2])),'.txt'])

    #pdb.set_trace()
    
    
    #-----------------------------------------------------------------#
    #               8. VARIANCE DECOMPOSITION ANALYSIS                #
    #-----------------------------------------------------------------#  
    
    ## TODO: change here options to vary analysis 
    vc = DirIndirVD(pheno = doto['pheno'], pheno_ID = doto['pheno_ID'], kinship_all = doto['kinship_full'], kinship_all_ID = doto['kinship_full_ID'], cage_all = doto['cage_full'], cage_all_ID = doto['cage_full_ID'], dam_all = doto['dam_full'], dam_all_ID = doto['dam_full_ID'], DGE = (DGE is not None), cageEffect = (cageEffect is not None), maternalEffect = (maternalEffect is not None), calc_ste=True, subset_IDs = doto['subset_IDs'], vc_init_type = None, vc_init = None)
    rv=vc.getOutput()

    if "optim_error" in rv.keys():
        assert rv['optim_error'] is not None, "optimization error is None, when expected to be something"
        print("Warning: optimization was not successful, have still saved sampleID, pheno and covs")
        toWrite=(trait1,'NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA')
    
    else:    
    #-----------------------------------------------------------------#
    #                    9. GETTING THE OUTPUT                        #
    #-----------------------------------------------------------------#  
        toWrite=(trait1,rv['sample_size'],rv['sample_size_all'],rv['var_Ad'],
            rv['var_As'],rv['STE_Ad'],rv['STE_As'],rv['total_var'],rv['corr_Ads'],
            rv['STE_corr_Ads'],rv['corr_params'],rv['conv'],rv['LML'],rv['var_Ed'],rv['var_Es'],
            rv['corr_Eds'],rv['var_C'],rv['var_M'],rv['env_rho'],rv['env_sigma2'])
  
        print(toWrite)
      
      
    #-----------------------------------------------------------------#
    #                 10. WRITING OUTPUT TO FILE                      #
    #-----------------------------------------------------------------#  
    VD_outfile=open(VD_outfile_name,'w')
    VD_outfile.write("\t".join(str(e) for e in toWrite)+'\n')
    VD_outfile.close()

print("Files saved to", VD_outfile_dir)
print("My program took", time.time() - start_time, "sec to run")

