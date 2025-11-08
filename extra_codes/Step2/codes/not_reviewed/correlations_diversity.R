#Load libraries
library("stringr")
library("dplyr")


# -----------------------------------------------------------------------------
# Functions

# 1. Read gzipped or flat files
read_gzipped <- function(input_file) {
  input <- input_file
  if (file_ext(input_file)=="gz") {
    input = gzfile(input_file)  
  }
  return(input)
}

read_tab_file <- function(input_file, header=FALSE) {
  file_content <- read.delim(read_gzipped(input_file), header=header, check.names=FALSE)		
  return(file_content)
}

read_csv_file <- function(input_file) {
  file_content <- read.csv(read_gzipped(input_file), stringsAsFactors = FALSE)	
  return(file_content)
}

# -----------------------------------------------------------------------------

# Extract metadata info
sample_identifier="host_subject_id"
input_meta="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/metadata_Shallow.txt" ### Shallow
if (input_meta != "") {
  print("Extracting metadata info")
  metadata<-read_tab_file(input_meta, TRUE)
  if (ncol(metadata)==1){
    metadata<-read_csv_file(input_meta)
  }
  if (length(which(colnames(metadata)==sample_identifier)) == 0) {
    stop("The metadata file has no column with the same name of the sample identifier informed")
  }
} else {
  stop("No metadata file provided")
}

#Shalow taxa
load("/nfs/users/abaud/fmorillo/P50/microbiome/output/microbiome_profiles/taxonomic_profile/gtdb_profiles/207/step2/abundance_mean/not_norm_comp_16S_207_E1/new_29-01-24/Matrix_Processing/filtered_clr_counts.RData")
ranks="Phylum,Class,Order,Family,Genus,Species"
ranks<-na.omit(c((sapply(strsplit(ranks, ","),"[",1)),
                 (sapply(strsplit(ranks, ","), "[",2)),
                 (sapply(strsplit(ranks, ","), "[",3)),
                 (sapply(strsplit(ranks, ","), "[",4)),
                 (sapply(strsplit(ranks, ","), "[",5)),
                 (sapply(strsplit(ranks, ","), "[",6))))
subset_ids=c("MI","NY","ALL")
permutations <- expand.grid(subset_ids,ranks)
colnames(permutations)=c("subset_ids","ranks")
concatenated <- paste(permutations$ranks, permutations$subset_ids, sep = "_")
for (i in 1:length(concatenated)){
  filtered_clr_counts=filtered_clr_counts_objs[[i]]
  assign(paste('filtered_clr_counts',concatenated[i],sep='_'), filtered_clr_counts)
}
all_traits_ALL=rbind(filtered_clr_counts_Phylum_ALL,
                     filtered_clr_counts_Class_ALL,
                     filtered_clr_counts_Order_ALL,
                     filtered_clr_counts_Family_ALL,
                     filtered_clr_counts_Genus_ALL,
                     filtered_clr_counts_Species_ALL)
all_traits_MI=rbind(filtered_clr_counts_Phylum_MI,
                    filtered_clr_counts_Class_MI,
                    filtered_clr_counts_Order_MI,
                    filtered_clr_counts_Family_MI,
                    filtered_clr_counts_Genus_MI,
                    filtered_clr_counts_Species_MI)
all_traits_NY=rbind(filtered_clr_counts_Phylum_NY,
                    filtered_clr_counts_Class_NY,
                    filtered_clr_counts_Order_NY,
                    filtered_clr_counts_Family_NY,
                    filtered_clr_counts_Genus_NY,
                    filtered_clr_counts_Species_NY)

#Shallow clusters
load("/nfs/users/abaud/fmorillo/P50/microbiome/output/microbiome_profiles/taxonomic_profile/gtdb_profiles/207/step2/abundance_mean/not_norm_comp_16S_207_E1/MicNet_cluster_new_30-01-24/with_outliers/Matrix_Processing/filtered_clr_counts.RData")
all_traits_MI=filtered_clr_counts_objs[[1]]
all_traits_NY=filtered_clr_counts_objs[[2]]
all_traits_ALL=filtered_clr_counts_objs[[3]]

#Include alpha diversity for each level
prefix <- "/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/"
load(paste(prefix,'Output_Alpha_TD/Matrix_Processing/alpha_diversity.RData',sep=''))
alpha_table_TD_MI=alpha_table_objs[[1]]
alpha_table_TD_NY=alpha_table_objs[[2]]
alpha_table_TD_ALL=alpha_table_objs[[3]]
load(paste(prefix,'Output_Alpha_PD/Matrix_Processing/alpha_diversity.RData',sep=''))
alpha_table_PD_MI=alpha_table_objs[[1]]
alpha_table_PD_NY=alpha_table_objs[[2]]
alpha_table_PD_ALL=alpha_table_objs[[3]]

# first remove from Species_matrix blanks (which are not in metadata), 
#All
motch = match(colnames(all_traits_ALL),metadata[,'host_subject_id'])
any(is.na(motch))
metadata_ALL = metadata[motch,]
#MI
motch = match(colnames(all_traits_MI),metadata[,'host_subject_id'])
any(is.na(motch))
metadata_MI = metadata[motch,]
#NY
motch = match(colnames(all_traits_NY),metadata[,'host_subject_id'])
any(is.na(motch))
metadata_NY = metadata[motch,]

for (study in c('MI','NY', 'ALL')) {
  correlations=data.frame()
  this_ALL_traits = get(paste('all_traits',study,sep='_'))
  this_matadata = get(paste('metadata',study,sep='_'))
  this_alpha_pd = get(paste('alpha_table_PD',study,sep='_'))
  #this_alpha_td = get(paste('alpha_table_TD',study,sep='_'))
  for(i in 1:length(rownames(this_ALL_traits))){
    test_spearman_wat = cor.test(this_ALL_traits[i,],this_matadata$weight_at_weaning, method = "spearman", exact=FALSE)
    test_spearman_gaw = cor.test(this_ALL_traits[i,],this_matadata$glucose_at_weaning, method = "spearman", exact=FALSE)
    test_spearman_hw = cor.test(this_ALL_traits[i,],this_matadata$host_weight, method = "spearman", exact=FALSE)
    test_spearman_gad = cor.test(this_ALL_traits[i,],this_matadata$glucose_at_dissection, method = "spearman", exact=FALSE)
    test_spearman_alpha_pd = cor.test(this_ALL_traits[i,],this_alpha_pd$qD, method = "spearman", exact=FALSE)
    #test_spearman_alpha_td = cor.test(this_ALL_traits[i,],this_alpha_td$qD, method = "spearman", exact=FALSE)
    new_row=data.frame(traits=rownames(this_ALL_traits)[i],
                       cor_weight_at_weaning=test_spearman_wat$estimate*100,
                       logP_weight_at_weaning=-log10(test_spearman_wat$p.value),
                       cor_glucose_at_weaning=test_spearman_gaw$estimate*100,
                       logP_glucose_at_weaning=-log10(test_spearman_gaw$p.value),
                       cor_host_weight=test_spearman_hw$estimate*100,
                       logP_host_weight=-log10(test_spearman_hw$p.value),
                       cor_glucose_at_dissection=test_spearman_gad$estimate*100,
                       logP_glucose_at_dissection=-log10(test_spearman_gad$p.value),
                       cor_alpha_pd=test_spearman_alpha_pd$estimate*100,
                       logP_alpha_pd=-log10(test_spearman_alpha_pd$p.value))
                       #cor_alpha_td=test_spearman_alpha_td$estimate*100,
                       #logP_alpha_td=-log10(test_spearman_alpha_td$p.value))
    correlations=rbind(correlations,new_row)
  }
  assign(paste('correlations',study,sep='_'), correlations)
  print(study)
}

load("/nfs/users/abaud/fmorillo/P50/microbiome/output/microbiome_profiles/taxonomic_profile/gtdb_profiles/207/step2/abundance_mean/not_norm_comp_16S_207_E1/new_29-01-24/Heritability/heritability.RData")
load("/nfs/users/abaud/fmorillo/P50/microbiome/output/microbiome_profiles/taxonomic_profile/gtdb_profiles/207/step2/abundance_mean/not_norm_comp_16S_207_E1/new_29-01-24/GWAS/All_QTLs.RData")

x= filtered_clr_counts_Species_ALL["s__Prevotella_sp002251245",]
y= filtered_clr_counts_Species_ALL["s__Acetatifactor_sp910589035",]

cor.test(x,y, method="spearman")
plot(x,y,xlab="s__Prevotella_sp002251245",ylab="s__Acetatifactor_sp910589035",xlim=c(-16,6),ylim=c(-16,6),main="Relative abundances (After CLR transformation)")

x= all_traits_ALL["cluster__13",]
y= all_traits_ALL["cluster__8",]

cor.test(x,y, method="spearman")
plot(x,y,xlab="Cluster 13",ylab="Cluster 8",
     #xlim=c(-11,9),ylim=c(-11,9),
     main="Relative abundances (After CLR transformation)")




######################

#Alpha diversity and host phenotypes
for (study in c('ALL','MI','NY')) {
  this_matadata = get(paste('metadata',study,sep='_'))
  #this_alpha_pd = get(paste('alpha_matrix_PD_q2',study,sep='_'))
  this_alpha_td = get(paste('alpha_table',study,sep='_'))
  #test_spearman_wat_pd = cor.test(this_alpha_pd$qPD,this_matadata$weight_at_weaning, method = "spearman", exact=FALSE)
  #test_spearman_gaw_pd = cor.test(this_alpha_pd$qPD,this_matadata$glucose_at_weaning, method = "spearman", exact=FALSE)
  #test_spearman_hw_pd = cor.test(this_alpha_pd$qPD,this_matadata$host_weight, method = "spearman", exact=FALSE)
  #test_spearman_gad_pd = cor.test(this_alpha_pd$qPD,this_matadata$glucose_at_dissection, method = "spearman", exact=FALSE)
  test_spearman_wat_td = cor.test(this_alpha_td$qD,this_matadata$weight_at_weaning, method = "spearman", exact=FALSE)
  test_spearman_gaw_td = cor.test(this_alpha_td$qD,this_matadata$glucose_at_weaning, method = "spearman", exact=FALSE)
  test_spearman_hw_td = cor.test(this_alpha_td$qD,this_matadata$host_weight, method = "spearman", exact=FALSE)
  test_spearman_gad_td = cor.test(this_alpha_td$qD,this_matadata$glucose_at_dissection, method = "spearman", exact=FALSE)
  #test_spearman_alphas = cor.test(this_alpha_td$qD,this_alpha_pd$qPD, method = "spearman", exact=FALSE)
  correlations_alpha=data.frame(#cor_weight_at_weaning_pd=test_spearman_wat_pd$estimate*100,
                                #logP_weight_at_weaning_pd=-log10(test_spearman_wat_pd$p.value),
                                #cor_glucose_at_weaning_pd=test_spearman_gaw_pd$estimate*100,
                                #logP_glucose_at_weaning_pd=-log10(test_spearman_gaw_pd$p.value),
                                #cor_host_weight_pd=test_spearman_hw_pd$estimate*100,
                                #logP_host_weight_pd=-log10(test_spearman_hw_pd$p.value),
                                #cor_glucose_at_dissection_pd=test_spearman_gad_pd$estimate*100,
                                #logP_glucose_at_dissection_pd=-log10(test_spearman_gad_pd$p.value),
                                #cor_weight_at_weaning_pd=test_spearman_wat_pd$estimate*100,
                                logP_weight_at_weaning_td=-log10(test_spearman_wat_td$p.value),
                                cor_glucose_at_weaning_td=test_spearman_gaw_td$estimate*100,
                                logP_glucose_at_weaning_td=-log10(test_spearman_gaw_td$p.value),
                                cor_host_weight_td=test_spearman_hw_td$estimate*100,
                                logP_host_weight_td=-log10(test_spearman_hw_td$p.value),
                                cor_glucose_at_dissection_td=test_spearman_gad_td$estimate*100,
                                logP_glucose_at_dissection_td=-log10(test_spearman_gad_td$p.value),
                                #cor_alphas=test_spearman_alphas$estimate*100,
                                #logP_alphas=-log10(test_spearman_alphas$p.value)
                                )
  assign(paste('correlations_alpha',study,sep='_'), correlations_alpha)
  print(study)
}

#prefix="/users/abaud/fmorillo/P50/microbiome/output/analyses/correlations/"
save(correlations_ALL,
     correlations_MI,
     correlations_NY,
     correlations_alpha_ALL,
     correlations_alpha_MI,
     correlations_alpha_NY,
     file = paste0(prefix,"correlations_phenotypes_gtdb_new_filter.RData")) 

#Individual Tests
plot(x = filtered_clr_counts_Species_ALL["s__Helicobacter_D_rodentium",], y = alpha_ALL$Species, frame = FALSE, xlab = "s__Helicobacter_D_rodentium (clr)", ylab = "Alpha diversity (after normalization)")
abline(lm(alpha_ALL$Species ~ filtered_clr_counts_Species_ALL["s__Helicobacter_D_rodentium",])) # abline(y ~ x)
cor.test(filtered_clr_counts_Species_ALL["s__Helicobacter_D_rodentium",],alpha_ALL$Species, method="spearman")


####Correlations s__Helicobacter_D_rodentium and other Species

for (study in c('ALL','MI','NY')) {
  correlations=data.frame()
  this_ALL_traits = get(paste('all_traits',study,sep='_'))
  for(i in 1:length(rownames(this_ALL_traits))){
    test_spearman_helic = cor.test(this_ALL_traits[i,],this_ALL_traits["s__Helicobacter_D_rodentium",], method = "spearman", exact=FALSE)
    new_row=data.frame(traits=rownames(this_ALL_traits)[i],
                       cor=test_spearman_helic$estimate*100,
                       logP=-log10(test_spearman_helic$p.value))
    correlations=rbind(correlations,new_row)
  }
  assign(paste('correlations_helico',study,sep='_'), correlations)
  print(study)
}

pvalues_dir = '/users/abaud/fmorillo/P50/microbiome/output/LOCO_files/taxonomic/LOCO_gtdb/pvalues/phenotypes/pruned_dosagesDGE_cageEffect/'
load(paste0(pvalues_dir,'QTLs_alpha1e-04.RData'))
load("/users/abaud/fmorillo/P50/microbiome/output/VD_files/taxonomic/VD_gtdb/heritability_gtdb.RData")
load("/users/abaud/fmorillo/P50/microbiome/output/analyses/correlations/correlations_phenotypes_gtdb.RData")
hist()


####

alpha_MI$study="MI"
alpha_NY$study="NY"
alphas=rbind(alpha_MI,alpha_NY)

helicobacter_MI=data.frame(abundance=all_traits_MI["s__Helicobacter_D_rodentium",],study="MI")
helicobacter_NY=data.frame(abundance=all_traits_NY["s__Helicobacter_D_rodentium",],study="NY")
helicobacters=rbind(helicobacter_MI,helicobacter_NY)

#Compare helicobacter
ggplot(alphas, aes(x=study, y=Species)) + 
  geom_boxplot() +
  xlab("") +
  ylab("Relative abundance of Helicobacter_D_rodentium") +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.bOrder = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=10)) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  geom_signif(comparisons = list(c("MI", "NY")), map_signif_level=TRUE)
#dev.off()

#Compare alpha diversity
ggplot(alphas, aes(x=study, y=Species)) + 
  geom_boxplot() +
  xlab("") +
  ylab("Alpha diversity (Species level)") +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.bOrder = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=10)) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  geom_signif(comparisons = list(c("MI", "NY")), map_signif_level=TRUE)
#dev.off()

boxplot(all_traits_MI["s__Helicobacter_D_rodentium",],all_traits_NY["s__Helicobacter_D_rodentium",])
boxplot(alpha_MI$Species, alpha_NY$Species)

x= filtered_clr_counts_Species_ALL["s__Paraprevotella_sp910584955",]
y= filtered_clr_counts_Species_ALL["s__Acetatifactor_sp910589035",]

cor.test(x,y, method="spearman")
plot(x,y,xlab="s__Paraprevotella_sp910584955",ylab="s__Acetatifactor_sp910589035",xlim=c(-16,6),ylim=c(-16,6))



myplot<-ggplot(metadata_ALL, aes(x=reorder(Study,host_weight,FUN = median), y=host_weight)) +
  geom_boxplot() +
  xlab("") +
  ylab("Weigth (g)") +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=10))+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  geom_signif(comparisons = list(unique(metadata_ALL$Study)),
              map_signif_level = TRUE)
print(myplot)
