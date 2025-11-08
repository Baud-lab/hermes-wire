library(rhdf5, lib.loc = "/users/abaud/fmorillo/R/x86_64-pc-linux-gnu-library/4.2/")
pheno="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/P50_rats_round8_PhenosGRMsDSdosages.h5"
#h5ls(pheno)
H5Fopen(pheno)
phenotypes=h5read(pheno,"/phenotypes/matrix")
rownames(phenotypes)=h5read(pheno,"/phenotypes/row_header/sample_ID")
colnames(phenotypes)=h5read(pheno,"/phenotypes/col_header/phenotype_ID")
H5Fclose(pheno)
save(phenotypes,file="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/phenotypes.RData")

#00                                          /     phenotypes   H5I_GROUP
#101                                /phenotypes     col_header   H5I_GROUP
#102                     /phenotypes/col_header   phenotype_ID H5I_DATASET
#103                                /phenotypes         matrix H5I_DATASET
#104                                /phenotypes     row_header   H5I_GROUP
#105                     /phenotypes/row_header           cage H5I_DATASET
#106                     /phenotypes/row_header      sample_ID H5I_DATASET


########

# Read Phenotypes
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/phenotypes.RData")

# Read Enterotypes
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Norm_nocluster_30k_full/Enterotypes/enterotypes.RData")
number=3
ent_samples=enterotypes_samples_objs[[number]]
ent_samples$sample=sapply(strsplit(ent_samples$sample, "_"),"[",1)

phenotypes_filt=phenotypes[rownames(phenotypes) %in% ent_samples$sample,]
pheno_list=data.frame(colnames(phenotypes_filt))
