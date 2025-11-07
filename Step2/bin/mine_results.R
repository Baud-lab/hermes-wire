#!/usr/bin/env Rscript

# Libraries
suppressMessages(library(argparse))

#Create parser object
parser <- ArgumentParser()

#Define desired outputs:
#GLOBAL FEATURES:
parser$add_argument("-input", "--input_cat", type="character", help="Input file")
parser$add_argument("-kin", "--kinship_version", type="character", help="Name of the kinship model",default="exp_DGA_7_Dec_2021")
parser$add_argument("-model", "--model_type", type="character", help="Name of the model",default="uni")
parser$add_argument("-pheno", "--phenotype_version", type="character", help="Phenotype version",default="bacterial_taxa")
parser$add_argument("-eff", "--effect", type="character", help="Effect to be checked")


#Get command line options, if help option encountered - print help and exit:
args <- parser$parse_args()


## Arguments
MODEL=args$model_type
EFF=args$effect
GRM_V=args$kinship_version     # `GRM_version` as saved in .h5 file
INPUT=args$input_cat
PHENO_V=args$phenotype_version   # `phenos_version` as saved in .h5 file 

#MODEL="uni"
#EFF="DGE,cageEffect,maternalEffect" # or just "DGE" or "cageEffect"; Note: if passing them to qsub, have to put "\" before the comma, EFF="DGE\,cageEffect"
#GRM_V="P50_Rn7_pruned"     # `GRM_version` as saved in .h5 file
#INPUT="~/output_full.txt"

## Read in all the output files coming from the heritability analysis
# irrelebant file=paste0(OUT_DIR,"/",MODEL,'variate/',PHENO_V,"/",GRM_V,"_",effect,"/input_concat.txt")
all_VCs = read.delim(INPUT, header = F, as.is = T)

## Give column and row names
coolnames = c('trait','sample_size','sample_size_all','var_Ad','var_As','STE_Ad','STE_As','total_var','corr_Ads','STE_corr_Ads','corr_params','conv','LML','var_Ed','var_Es','corr_Eds','var_C','var_M','env_rho','env_sigma2')  
if (length(coolnames) != dim(all_VCs)[2]) stop('pb')
colnames(all_VCs)=coolnames
rownames(all_VCs)=all_VCs[,'trait']

## The previous python scripts encode missing data as -999; R needs NA; so we change here
for (col in colnames(all_VCs)) {
  all_VCs[which(all_VCs[,col]==(-999)),col]=NA
}

## Should always be true
print(paste('All models converged: ',all(all_VCs[,'conv']=='True'),sep=''))

## Normalisation procedure - essential
all_VCs[,c('var_Ad','var_As','var_Ed','var_Es','var_C','var_M')] = all_VCs[,c('var_Ad','var_As','var_Ed','var_Es','var_C','var_M')] / all_VCs[,'total_var']

## Save output as R object for further analyses such as plotting etc
all_VCs=all_VCs[order(all_VCs[,'var_Ad'],decreasing=T),]

#EFF="DGE,cageEffect" # or just "DGE" or "cageEffect"; Note: if passing them to qsub, have to put "\" before the comma, EFF="DGE\,cageEffect"
if (EFF == 'DGE,cageEffect,maternalEffect') effect = 'DGE_cageEffect_maternalEffect' else if (EFF == 'DGE,maternalEffect') effect = 'DGE_maternalEffect' else if (EFF == 'DGE,cageEffect') effect = 'DGE_cageEffect' else if (EFF == 'cageEffect,maternalEffect') effect = 'cageEffect_maternalEffect' else if (EFF == 'cageEffect') effect = 'cageEffect' else if (EFF == 'maternalEffect') effect = 'maternalEffect'

save(all_VCs,file=paste(MODEL,effect, PHENO_V, GRM_V,'all_VCs.RData', sep='_')) #have the values of the params in the file name.

