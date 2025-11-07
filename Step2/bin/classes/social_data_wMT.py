
import numpy as np
import h5py
#re for regular expressions
import re
import pdb

class SocialData():
    
    def __init__(self, in_file=None, phenos_version = None,covs_version=None, cage_version=None, dam_version=None, GRM_version = None, subset = None, chr = None):
        #data = SocialData(in_file,task,covs_task,cage_task,dam_task,kinship_type,subset) # order in exp_bivar impo!
        assert in_file is not None, 'Give me an input!'
        assert phenos_version is not None, 'Specify phenotype type! (phenos_version)'
        assert covs_version is not None, 'Specify covariate type! (covs_version)'
        assert cage_version is not None, 'Specify cage type! (cage_version)'
        assert dam_version is not None, 'Specify dam type! (dam_version)'
        assert GRM_version is not None, 'Specify kinship type!' 
        
        self.phenos_version=phenos_version #this is tied to the type of phenotypes used and the normalisation procedure. 
        
        self.covs_version=covs_version #this is tied to the type of covs used and the normalisation procedure.
        self.cage_version=cage_version
        self.dam_version=dam_version
        self.in_file=in_file 
        self.GRM_version=GRM_version #GRM_version is based on the set of genotypes used (e.g. pruned round 8, pruned round 2 or unpruned round 1). this determines the GRM, the GRM_LOCO (and also the genotype table used in GWAS although for some reason this is read separately in teh GWAS code...hm hm)
        
        self.subset = subset # default is None -> i.e. all 
        self.chr = chr 
        self.load()
    
    def load(self):
        
        print("Reading data from ", self.in_file)
        
        f = h5py.File(self.in_file,'r')
        
        # get phenotypes
        self.all_pheno = f['phenotypes'][self.phenos_version]['matrix'][:].T
        self.measures = f['phenotypes'][self.phenos_version]['col_header']['phenotype_ID'].asstr()[:] # as string
        self.pheno_ID = f['phenotypes'][self.phenos_version]['row_header']['sample_ID'].asstr()[:]

        # get covariates - NB: saved in 'covariates' element in h5
        if 'covariates' not in f.keys():
            print("No covariates in .h5, all_covs2use to None")
            self.all_covs2use = None
        else:
            self.all_covs2use = f['phenotypes'][self.phenos_version]['col_header']['covariatesUsed'].asstr()[:] 
            self.all_covs = f['covariates'][self.covs_version]['matrix'][:].T
            self.covariates = f['covariates'][self.covs_version]['col_header']['covariate_ID'].asstr()[:]
            self.covs_ID = f['covariates'][self.covs_version]['row_header']['sample_ID'].asstr()[:]
        
        # get cages - NB: we decided to dissociate information about the cage from information about the phenotypes
        self.cage_full = f['cages'][self.cage_version]['array'].asstr()[:] #self
        self.cage_full_ID = f['cages'][self.cage_version]['sample_ID'].asstr()[:]
        print('theres cage info in HDF5')
        
        # get dams - NB: we decided to dissociate information about the dam from information about the phenotypes
        self.dam_full = f['dams'][self.dam_version]['array'].asstr()[:] #self
        self.dam_full_ID = f['dams'][self.dam_version]['sample_ID'].asstr()[:]
        print('theres dam info in HDF5')


        if len(self.all_pheno.shape)==1: #very unlikely given highly dimensional datasets we use but still...
            self.all_pheno = self.all_pheno[:,np.newaxis]

        if self.all_covs2use is not None:
            if len(self.all_covs.shape)==1: #very unlikely given highly dimensional datasets we use but still...
                self.all_covs = self.all_covs[:,np.newaxis]


        if self.chr is not None:
            print('social data in LOCO')
            #DO NOT swap [self.GRM_version] to before any GRM thing like below actually # FOR FELIPE AND AMELIE: this have no idea!!
            self.kinship_full = f['GRMs_LOCO'][self.GRM_version][''.join(['chr',str(self.chr)])]['matrix'][:]
            self.kinship_full_ID = f['GRMs_LOCO'][self.GRM_version][''.join(['chr',str(self.chr)])]['row_header']['sample_ID'].asstr()[:]
        else:    
            self.kinship_full = f['GRM'][self.GRM_version]['matrix'][:]
            self.kinship_full_ID = f['GRM'][self.GRM_version]['row_header']['sample_ID'].asstr()[:]

        if self.subset is None:
            self.subset_IDs = self.kinship_full_ID
        else:
            self.subset_IDs = f['subsets'][self.subset].asstr()[:]


    def get_data(self,col,col_MT = None):
        
        self.trait = self.measures[col]

        if col_MT is None:
            self.pheno = self.all_pheno[:,[col]] #using array as index permits keeping the dimensionality kind of
            self.trait_MT = None
        else:
            self.trait_MT = self.measures[col_MT]
            self.pheno = np.concatenate([self.all_pheno[:,[col]], self.all_pheno[:,[col_MT]]],0)

        #that's if no covs in entire study
        if self.all_covs2use is None:
            self.covs_ID = None # This is defined above if there's 'covariates' group - here set to None 
            self.covs = None
            covariates_names = None
          
        else: #at this point it is possible trait1 and/or trait2 have empty covs2use
            covs2use = self.all_covs2use[col].split(',')
            Ic = np.zeros(self.covariates.shape[0],dtype=bool) 
            for cov in covs2use:
                if cov != '':
                    assert any(self.covariates==cov), 'covariate not in cov table'
                    Ic = np.logical_or(Ic,self.covariates==cov)
            
            covariates_names = self.covariates[Ic] # covariates_names will be empty list rather than None if no cov for that phenotype (col)
            print('Initial covs (for first phenotype) in social_data are ' + str(covariates_names))
            if any(Ic): #Ic has only 1 value as len is 1 
                self.covs = self.all_covs[:,Ic] #self.covs is covs for phenotype in column col; since Ic created as an array no need for [Ic]
            else:
                self.covs = None
            
            if col_MT is not None: #should only go there if there are covs in the study
                #at this point covs is None or sthg; covs_ID is sthg; covariates_names is (potentially null) array        

                if self.covs is not None:
                    self.covs = np.append(self.covs,np.zeros(self.covs.shape),0)
                #so for case where phenotype in col has no covs2use we now have a column of 1 then 0 of length 2xN
                #self.covs_ID = np.concatenate([self.covs_ID,self.covs_ID])
                
                # now doing 3 for second phenotype
                covs2use_MT = self.all_covs2use[col_MT].split(',')
                Ic_MT = np.zeros(self.covariates.shape[0],dtype=bool)
                for cov_MT in covs2use_MT:
                    if cov_MT != '': 
                        assert any(self.covariates==cov_MT), 'covariate not in cov table'
                        Ic_MT = np.logical_or(Ic_MT,self.covariates==cov_MT)
                covariates_names_MT = self.covariates[Ic_MT]
                print('Initial covs in social_data for second phenotype are ' + str(covariates_names_MT))
                if any(Ic_MT): #Ic has only 1 value as len is 1 
                    covs_MT = self.all_covs[:,Ic_MT] #self.covs is covs for phenotype in column col
                    covs_MT = np.append(np.zeros(covs_MT.shape),covs_MT,0)
                else:
                    covs_MT = None
                if self.covs is not None and covs_MT is not None:
                    self.covs = np.append(self.covs,covs_MT,1)
                    covariates_names = np.append(covariates_names,covariates_names_MT)
                elif self.covs is None and covs_MT is not None:
                    self.covs = covs_MT
                    covariates_names = covariates_names_MT
            
            
        if self.covs is None:
            self.covs_ID = None
            assert covariates_names is None, "pb in covariates_names"
        
        return {'trait' : self.trait,
                'trait_MT' : self.trait_MT,
                'pheno' : self.pheno,
                'pheno_ID' : self.pheno_ID,
                'covs' : self.covs, 
                'covs_ID' : self.covs_ID, 
                'covariates_names' : covariates_names,
                'kinship_type' : self.GRM_version, 
                'kinship_full' : self.kinship_full,
                'kinship_full_ID' : self.kinship_full_ID, 
                'cage_full' : self.cage_full, 
                'cage_full_ID' : self.cage_full_ID,
                'dam_full' : self.dam_full, 
                'dam_full_ID' : self.dam_full_ID,
                'subset_IDs' : self.subset_IDs}






