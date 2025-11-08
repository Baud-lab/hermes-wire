library(xfun)
library(stats)
library(dplyr)
library(tidyverse)
subset_id="ALL"
# "ALL" "MI" "NY"
method="Shallow"

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

metadata <- read.delim("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/metadata_Shallow.txt")
#metadata=metadata[metadata$host_subject_id %in% colnames(matrix),]
#count(metadata,Study,Sex)

# Calculate relative abundances and prevalence
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Func_Mat_Net/Matrix_Processing/filtered_prev.RData")
number=4
if (subset_id=="MI"){
  number=number-2
} else {
  if (subset_id=="NY"){
    number=number-1
  }
}
filtered_prev_ec4=filtered_prev_objs[[number]]
ab_vs_prev_ec4=data.frame(trait=rownames(filtered_prev_ec4),
                      abundance=100*apply(sweep(filtered_prev_ec4,2,colSums(filtered_prev_ec4),"/"), 1, mean, na.rm=TRUE),
                      prevalence=100*rowSums(filtered_prev_ec4>0)/ncol(filtered_prev_ec4))

#Genus relative abundances
number=3
if (subset_id=="MI"){
  number=number-2
} else {
  if (subset_id=="NY"){
    number=number-1
  }
}
filtered_prev_ec3=filtered_prev_objs[[number]]
ab_vs_prev_ec3=data.frame(trait=rownames(filtered_prev_ec3),
                            abundance=100*apply(sweep(filtered_prev_ec3,2,colSums(filtered_prev_ec3),"/"), 1, mean, na.rm=TRUE),
                            prevalence=100*rowSums(filtered_prev_ec3>0)/ncol(filtered_prev_ec3))

# Family relative abundances
number=2
if (subset_id=="MI"){
  number=number-2
} else {
  if (subset_id=="NY"){
    number=number-1
  }
}
filtered_prev_ec2=filtered_prev_objs[[number]]
ab_vs_prev_ec2=data.frame(trait=rownames(filtered_prev_ec2),
                             abundance=100*apply(sweep(filtered_prev_ec2,2,colSums(filtered_prev_ec2),"/"), 1, mean, na.rm=TRUE),
                             prevalence=100*rowSums(filtered_prev_ec2>0)/ncol(filtered_prev_ec2))

# Family relative abundances
number=1
if (subset_id=="MI"){
  number=number-2
} else {
  if (subset_id=="NY"){
    number=number-1
  }
}
filtered_prev_ec1=filtered_prev_objs[[number]]
ab_vs_prev_ec1=data.frame(trait=rownames(filtered_prev_ec1),
                             abundance=100*apply(sweep(filtered_prev_ec1,2,colSums(filtered_prev_ec1),"/"), 1, mean, na.rm=TRUE),
                             prevalence=100*rowSums(filtered_prev_ec1>0)/ncol(filtered_prev_ec1))


# Correlations with beta and alpha diversity

## Load residuals
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Func_Mat_Net/Matrix_Processing/residuals_qned_counts.RData")
number=4
if (subset_id=="MI"){
  number=number-2
} else {
  if (subset_id=="NY"){
    number=number-1
  }
}
all_traits_ec4=residuals_qned_counts_objs[[number]]

metadata=metadata[metadata$host_subject_id %in% colnames(all_traits_ec4),]
# row‐wise variances (1 = rows, 2 = columns)
variances <- apply(all_traits_ec4, 1, var, na.rm = TRUE)
variances_ec4 <- data.frame(
  EC    = rownames(all_traits_ec4),
  Variance = variances,
  row.names = NULL
)

ec_numbers <- sapply(strsplit(rownames(all_traits_ec4), "_"),"[",3)
ec_components <- data.frame(do.call(rbind, strsplit(ec_numbers, "\\.")))
colnames(ec_components) <- c("L1", "L2", "L3", "L4")
ec_components$Level1 <- ec_components$L1
ec_components$Level2 <- with(ec_components, paste(Level1, L2, sep = "."))
ec_components$Level3 <- with(ec_components, paste(Level2, L3, sep = "."))
ec_components$Level4 <- with(ec_components, paste(Level3, L4, sep = "."))
ec_components=ec_components[,5:8]
rownames(ec_components)=paste0("EC4__",ec_components$Level4)

variances_ec4$Level1=ec_components$Level1[match(variances_ec4$EC,rownames(ec_components),)]
#variances_ec4$Level2=ec_components$Level2[match(variances_ec4$EC,rownames(ec_components),)]
#variances_ec4$Level3=ec_components$Level3[match(variances_ec4$EC,rownames(ec_components),)]

variances_ec4$Level1[which(variances_ec4$Level1=="1")]="Oxidoreductases"
variances_ec4$Level1[which(variances_ec4$Level1=="2")]="Transferases"
variances_ec4$Level1[which(variances_ec4$Level1=="3")]="Hydrolases"
variances_ec4$Level1[which(variances_ec4$Level1=="4")]="Lyases"
variances_ec4$Level1[which(variances_ec4$Level1=="5")]="Isomerases"
variances_ec4$Level1[which(variances_ec4$Level1=="6")]="Ligases"
variances_ec4$Level1[which(variances_ec4$Level1=="7")]="Translocases"


myplot<-ggplot(variances_ec4, aes(x=reorder(Level1,Variance,FUN = median), y=Variance)) +
  geom_boxplot() +
  xlab("") +
  ylab("Variance") +
  #ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Translocases","Oxidoreductases")),
              map_signif_level = TRUE)
print(myplot)

#na.omit(metadata$Age[metadata$Study=="MI"])
#na.omit(metadata$Age[metadata$Study=="NY"])

number=1
if (subset_id=="MI"){
  number=number-2
} else {
  if (subset_id=="NY"){
    number=number-1
  }
}
all_traits_ec1=residuals_qned_counts_objs[[number]]
# row‐wise variances (1 = rows, 2 = columns)
variances <- apply(all_traits_ec1, 1, var, na.rm = TRUE)
variances_ec1 <- data.frame(
  EC    = rownames(all_traits_ec1),
  Variance = variances,
  row.names = NULL
)


number=2
if (subset_id=="MI"){
  number=number-2
} else {
  if (subset_id=="NY"){
    number=number-1
  }
}
all_traits_class=residuals_qned_counts_objs[[number]]
variances <- apply(all_traits_species, 1, var, na.rm = TRUE)
variances_species <- data.frame(
  Species    = rownames(all_traits_species),
  Variance = variances,
  row.names = NULL
)
number=3
if (subset_id=="MI"){
  number=number-2
} else {
  if (subset_id=="NY"){
    number=number-1
  }
}
all_traits_order=residuals_qned_counts_objs[[number]]
number=4
if (subset_id=="MI"){
  number=number-2
} else {
  if (subset_id=="NY"){
    number=number-1
  }
}
all_traits_family=residuals_qned_counts_objs[[number]]
number=5
if (subset_id=="MI"){
  number=number-2
} else {
  if (subset_id=="NY"){
    number=number-1
  }
}
all_traits_genus=residuals_qned_counts_objs[[number]]

## Load cluster residuals
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Func_Mat_Net/Cluster_Analyses/residuals_qned_counts_clusters.RData")
number=1
if (subset_id=="MI"){
  number=number-2
} else {
  if (subset_id=="NY"){
    number=number-1
  }
}
all_traits_cluster=residuals_qned_counts_clusters_objs[[number]]

correlations_cluster3_clusters=data.frame()
for(i in 1:length(rownames(all_traits_cluster))){
  test_pearson = cor.test(all_traits_cluster[i,],all_traits_cluster[which(rownames(all_traits_cluster)=="cluster__3"),], method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_cluster)[i],correlations=test_pearson$estimate,pval=-log(test_pearson$p.value),pvalue=test_pearson$p.value)
  correlations_cluster3_clusters=rbind(correlations_cluster3_clusters,new_row)
}

## Load beta
#load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Func_Mat_Net/Diversity/beta_diversity.RData")
#beta=gp.dist_objs[[3]]
#betaoa <- ape::pcoa(beta)
#beta <- data.frame(betaoa$vectors)

load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Func_Mat_Net/Diversity/residuals_qned_counts_beta.RData")
number=1
if (subset_id=="MI"){
  number=number-2
} else {
  if (subset_id=="NY"){
    number=number-1
  }
}
beta=residuals_qned_counts_beta_objs[[number]]
beta=as.matrix(beta)
rownames(beta)
#[1] "beta__TD_dispersion" "beta__TD_PC1"        "beta__TD_PC2"       
#[4] "beta__PD_dispersion" "beta__PD_PC1"        "beta__PD_PC2"   
beta_df=data.frame(Axis.1=beta[2,],Axis.2=beta[3,])
beta=beta[2,]

## Load alpha
#load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Func_Mat_Net/Diversity/residuals_qned_counts_alpha.RData")
load("~/NextFlow/Microbiome_Profiling/step2/Output_Func_Mat_Net/Diversity/residuals_qned_counts_alpha.RData")
number=1
if (subset_id=="MI"){
  number=number-2
} else {
  if (subset_id=="NY"){
    number=number-1
  }
}
alpha=residuals_qned_counts_alpha_objs[[number]]
rownames(alpha)
#[1] "alpha__TD_q0" "alpha__TD_q1" "alpha__TD_q2" "alpha__PD_q0" "alpha__PD_q1"
#[6] "alpha__PD_q2"
alpha=alpha[6,]

# Glucose ######### Wrong
#glucose=metadata$glucose_at_dissection[match(colnames(all_traits_species),metadata$host_subject_id)]
#names(glucose)=colnames(all_traits_species)
print("Load phenotypes residuals")
load("/users/abaud/fmorillo/paper_figures/pehotypes_residuals.RData")
glucose=df$glucose_at_dissection[match(colnames(all_traits_species),rownames(df))]
names(glucose)=colnames(all_traits_species)

# Weigth
weight=metadata$Weight[match(colnames(all_traits_species),metadata$host_subject_id)]
names(weight)=colnames(all_traits_species)

# Calculate correlation between residuals and beta diversity
correlations_beta_species=data.frame()
for(i in 1:length(rownames(all_traits_species))){
  test_pearson_alpha = cor.test(all_traits_species[i,],beta, method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_species)[i],estimate=test_pearson_alpha$estimate,pval=-log(test_pearson_alpha$p.value),pvalue=test_pearson_alpha$p.value)
  correlations_beta_species=rbind(correlations_beta_species,new_row)
}

correlations_beta_cluster=data.frame()
for(i in 1:length(rownames(all_traits_cluster))){
  test_pearson_alpha = cor.test(all_traits_cluster[i,],beta, method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_cluster)[i],estimate=test_pearson_alpha$estimate,pval=-log(test_pearson_alpha$p.value),pvalue=test_pearson_alpha$p.value)
  correlations_beta_cluster=rbind(correlations_beta_cluster,new_row)
}

correlations_beta_phylum=data.frame()
for(i in 1:length(rownames(all_traits_phylum))){
  test_pearson_alpha = cor.test(all_traits_phylum[i,],beta, method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_phylum)[i],estimate=test_pearson_alpha$estimate,pval=-log(test_pearson_alpha$p.value),pvalue=test_pearson_alpha$p.value)
  correlations_beta_phylum=rbind(correlations_beta_phylum,new_row)
}

correlations_beta_class=data.frame()
for(i in 1:length(rownames(all_traits_phylum))){
  test_pearson_alpha = cor.test(all_traits_class[i,],beta, method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_class)[i],estimate=test_pearson_alpha$estimate,pval=-log(test_pearson_alpha$p.value),pvalue=test_pearson_alpha$p.value)
  correlations_beta_class=rbind(correlations_beta_class,new_row)
}

correlations_beta_order=data.frame()
for(i in 1:length(rownames(all_traits_order))){
  test_pearson_alpha = cor.test(all_traits_order[i,],beta, method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_order)[i],estimate=test_pearson_alpha$estimate,pval=-log(test_pearson_alpha$p.value),pvalue=test_pearson_alpha$p.value)
  correlations_beta_order=rbind(correlations_beta_order,new_row)
}

correlations_beta_family=data.frame()
for(i in 1:length(rownames(all_traits_family))){
  test_pearson_alpha = cor.test(all_traits_family[i,],beta, method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_family)[i],estimate=test_pearson_alpha$estimate,pval=-log(test_pearson_alpha$p.value),pvalue=test_pearson_alpha$p.value)
  correlations_beta_family=rbind(correlations_beta_family,new_row)
}

correlations_beta_genus=data.frame()
for(i in 1:length(rownames(all_traits_genus))){
  test_pearson_alpha = cor.test(all_traits_genus[i,],beta, method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_genus)[i],estimate=test_pearson_alpha$estimate,pval=-log(test_pearson_alpha$p.value),pvalue=test_pearson_alpha$p.value)
  correlations_beta_genus=rbind(correlations_beta_genus,new_row)
}

correlations_beta_all=rbind(
  correlations_beta_species,
  correlations_beta_genus,
  correlations_beta_family,
  correlations_beta_order,
  correlations_beta_class,
  correlations_beta_phylum,
  correlations_beta_cluster
)
correlations_beta_all$bonf_pvalue = correlations_beta_all$pvalue*dim(correlations_beta_all)[1]
correlations_beta_all$qvalue =p.adjust(correlations_beta_all$pvalue, method = "fdr")
correlations_beta_all <- correlations_beta_all %>% 
  mutate( Significance = case_when(
    qvalue < 0.1 & bonf_pvalue < 0.05 ~ "Significant (Bonferroni)",
    qvalue < 0.1 & bonf_pvalue >= 0.05 ~ "Significant (FDR)",
    TRUE ~ "Non-significant"))
correlations_beta_all$Direction <- ifelse(correlations_beta_all$estimate > 0, "Positive", "Negative")
correlations_beta_all = correlations_beta_all[order(correlations_beta_all$pval, decreasing = T),]

# Calculate correlation between residuals and alpha diversity
correlations_alpha_species=data.frame()
for(i in 1:length(rownames(all_traits_species))){
  test_pearson_alpha = cor.test(all_traits_species[i,],alpha, method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_species)[i],estimate=test_pearson_alpha$estimate,pval=-log(test_pearson_alpha$p.value),pvalue=test_pearson_alpha$p.value)
  correlations_alpha_species=rbind(correlations_alpha_species,new_row)
}

correlations_alpha_cluster=data.frame()
for(i in 1:length(rownames(all_traits_cluster))){
  test_pearson_alpha = cor.test(all_traits_cluster[i,],alpha, method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_cluster)[i],estimate=test_pearson_alpha$estimate,pval=-log(test_pearson_alpha$p.value),pvalue=test_pearson_alpha$p.value)
  correlations_alpha_cluster=rbind(correlations_alpha_cluster,new_row)
}

correlations_alpha_phylum=data.frame()
for(i in 1:length(rownames(all_traits_phylum))){
  test_pearson_alpha = cor.test(all_traits_phylum[i,],alpha, method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_phylum)[i],estimate=test_pearson_alpha$estimate,pval=-log(test_pearson_alpha$p.value),pvalue=test_pearson_alpha$p.value)
  correlations_alpha_phylum=rbind(correlations_alpha_phylum,new_row)
}

correlations_alpha_class=data.frame()
for(i in 1:length(rownames(all_traits_phylum))){
  test_pearson_alpha = cor.test(all_traits_class[i,],alpha, method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_class)[i],estimate=test_pearson_alpha$estimate,pval=-log(test_pearson_alpha$p.value),pvalue=test_pearson_alpha$p.value)
  correlations_alpha_class=rbind(correlations_alpha_class,new_row)
}

correlations_alpha_order=data.frame()
for(i in 1:length(rownames(all_traits_order))){
  test_pearson_alpha = cor.test(all_traits_order[i,],alpha, method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_order)[i],estimate=test_pearson_alpha$estimate,pval=-log(test_pearson_alpha$p.value),pvalue=test_pearson_alpha$p.value)
  correlations_alpha_order=rbind(correlations_alpha_order,new_row)
}

correlations_alpha_family=data.frame()
for(i in 1:length(rownames(all_traits_family))){
  test_pearson_alpha = cor.test(all_traits_family[i,],alpha, method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_family)[i],estimate=test_pearson_alpha$estimate,pval=-log(test_pearson_alpha$p.value),pvalue=test_pearson_alpha$p.value)
  correlations_alpha_family=rbind(correlations_alpha_family,new_row)
}

correlations_alpha_genus=data.frame()
for(i in 1:length(rownames(all_traits_genus))){
  test_pearson_alpha = cor.test(all_traits_genus[i,],alpha, method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_genus)[i],estimate=test_pearson_alpha$estimate,pval=-log(test_pearson_alpha$p.value),pvalue=test_pearson_alpha$p.value)
  correlations_alpha_genus=rbind(correlations_alpha_genus,new_row)
}

correlations_alpha_all=rbind(
  correlations_alpha_species,
  correlations_alpha_genus,
  correlations_alpha_family,
  correlations_alpha_order,
  correlations_alpha_class,
  correlations_alpha_phylum,
  correlations_alpha_cluster
)
correlations_alpha_all$bonf_pvalue = correlations_alpha_all$pvalue*dim(correlations_alpha_all)[1]
correlations_alpha_all$qvalue =p.adjust(correlations_alpha_all$pvalue, method = "fdr")
correlations_alpha_all <- correlations_alpha_all %>% 
  mutate( Significance = case_when(
    qvalue < 0.1 & bonf_pvalue < 0.05 ~ "Significant (Bonferroni)",
    qvalue < 0.1 & bonf_pvalue >= 0.05 ~ "Significant (FDR)",
    TRUE ~ "Non-significant"))
correlations_alpha_all$Direction <- ifelse(correlations_alpha_all$estimate > 0, "Positive", "Negative")
correlations_alpha_all = correlations_alpha_all[order(correlations_alpha_all$pval, decreasing = T),]

# Correlations glucose
correlations_glucose_species=data.frame()
for(i in 1:length(rownames(all_traits_species))){
  test_pearson_glucose = cor.test(all_traits_species[i,],glucose, method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_species)[i],estimate=test_pearson_glucose$estimate,pval=-log(test_pearson_glucose$p.value),pvalue=test_pearson_glucose$p.value)
  correlations_glucose_species=rbind(correlations_glucose_species,new_row)
}

correlations_glucose_cluster=data.frame()
for(i in 1:length(rownames(all_traits_cluster))){
  test_pearson_glucose = cor.test(all_traits_cluster[i,],glucose, method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_cluster)[i],estimate=test_pearson_glucose$estimate,pval=-log(test_pearson_glucose$p.value),pvalue=test_pearson_glucose$p.value)
  correlations_glucose_cluster=rbind(correlations_glucose_cluster,new_row)
}

correlations_glucose_phylum=data.frame()
for(i in 1:length(rownames(all_traits_phylum))){
  test_pearson_glucose = cor.test(all_traits_phylum[i,],glucose, method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_phylum)[i],estimate=test_pearson_glucose$estimate,pval=-log(test_pearson_glucose$p.value),pvalue=test_pearson_glucose$p.value)
  correlations_glucose_phylum=rbind(correlations_glucose_phylum,new_row)
}

correlations_glucose_class=data.frame()
for(i in 1:length(rownames(all_traits_phylum))){
  test_pearson_glucose = cor.test(all_traits_class[i,],glucose, method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_class)[i],estimate=test_pearson_glucose$estimate,pval=-log(test_pearson_glucose$p.value),pvalue=test_pearson_glucose$p.value)
  correlations_glucose_class=rbind(correlations_glucose_class,new_row)
}

correlations_glucose_order=data.frame()
for(i in 1:length(rownames(all_traits_order))){
  test_pearson_glucose = cor.test(all_traits_order[i,],glucose, method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_order)[i],estimate=test_pearson_glucose$estimate,pval=-log(test_pearson_glucose$p.value),pvalue=test_pearson_glucose$p.value)
  correlations_glucose_order=rbind(correlations_glucose_order,new_row)
}

correlations_glucose_family=data.frame()
for(i in 1:length(rownames(all_traits_family))){
  test_pearson_glucose = cor.test(all_traits_family[i,],glucose, method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_family)[i],estimate=test_pearson_glucose$estimate,pval=-log(test_pearson_glucose$p.value),pvalue=test_pearson_glucose$p.value)
  correlations_glucose_family=rbind(correlations_glucose_family,new_row)
}

correlations_glucose_genus=data.frame()
for(i in 1:length(rownames(all_traits_genus))){
  test_pearson_glucose = cor.test(all_traits_genus[i,],glucose, method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_genus)[i],estimate=test_pearson_glucose$estimate,pval=-log(test_pearson_glucose$p.value),pvalue=test_pearson_glucose$p.value)
  correlations_glucose_genus=rbind(correlations_glucose_genus,new_row)
}

#################################################

lm_glucose_cluster=data.frame()
for(i in 1:length(rownames(all_traits_cluster))){
  print(rownames(all_traits_cluster)[i])
  my_data=data.frame(glucose=glucose,trait=all_traits_cluster[i,],beta=beta)
  model_beta=summary(lm(trait ~ beta, data=my_data))$coefficients
  model_gluc=summary(lm(glucose ~ trait, data=my_data))$coefficients
  new_row=data.frame(trait=rownames(all_traits_cluster)[i],
                     est_beta_to_trait=model_beta[2,1],
                     est_trait_to_gluc=model_gluc[2,1],
                     final_effect=model_beta[2,1]*model_gluc[2,1],
                     pval_beta_to_trait=model_beta[2,4],
                     pval_trait_to_gluc=model_gluc[2,4])
  lm_glucose_cluster=rbind(lm_glucose_cluster,new_row)
}

lm_glucose_phylum=data.frame()
for(i in 1:length(rownames(all_traits_phylum))){
  print(rownames(all_traits_phylum)[i])
  my_data=data.frame(glucose=glucose,trait=all_traits_phylum[i,],beta=beta)
  model_beta=summary(lm(trait ~ beta, data=my_data))$coefficients
  model_gluc=summary(lm(glucose ~ trait, data=my_data))$coefficients
  new_row=data.frame(trait=rownames(all_traits_phylum)[i],
                     est_beta_to_trait=model_beta[2,1],
                     est_trait_to_gluc=model_gluc[2,1],
                     final_effect=model_beta[2,1]*model_gluc[2,1],
                     pval_beta_to_trait=model_beta[2,4],
                     pval_trait_to_gluc=model_gluc[2,4])
  lm_glucose_phylum=rbind(lm_glucose_phylum,new_row)
}

lm_glucose_class=data.frame()
for(i in 1:length(rownames(all_traits_class))){
  print(rownames(all_traits_class)[i])
  my_data=data.frame(glucose=glucose,trait=all_traits_class[i,],beta=beta)
  model_beta=summary(lm(trait ~ beta, data=my_data))$coefficients
  model_gluc=summary(lm(glucose ~ trait, data=my_data))$coefficients
  new_row=data.frame(trait=rownames(all_traits_class)[i],
                     est_beta_to_trait=model_beta[2,1],
                     est_trait_to_gluc=model_gluc[2,1],
                     final_effect=model_beta[2,1]*model_gluc[2,1],
                     pval_beta_to_trait=model_beta[2,4],
                     pval_trait_to_gluc=model_gluc[2,4])
  lm_glucose_class=rbind(lm_glucose_class,new_row)
}

lm_glucose_order=data.frame()
for(i in 1:length(rownames(all_traits_order))){
  print(rownames(all_traits_order)[i])
  my_data=data.frame(glucose=glucose,trait=all_traits_order[i,],beta=beta)
  model_beta=summary(lm(trait ~ beta, data=my_data))$coefficients
  model_gluc=summary(lm(glucose ~ trait, data=my_data))$coefficients
  new_row=data.frame(trait=rownames(all_traits_order)[i],
                     est_beta_to_trait=model_beta[2,1],
                     est_trait_to_gluc=model_gluc[2,1],
                     final_effect=model_beta[2,1]*model_gluc[2,1],
                     pval_beta_to_trait=model_beta[2,4],
                     pval_trait_to_gluc=model_gluc[2,4])
  lm_glucose_order=rbind(lm_glucose_order,new_row)
}

lm_glucose_family=data.frame()
for(i in 1:length(rownames(all_traits_family))){
  print(rownames(all_traits_family)[i])
  my_data=data.frame(glucose=glucose,trait=all_traits_family[i,],beta=beta)
  model_beta=summary(lm(trait ~ beta, data=my_data))$coefficients
  model_gluc=summary(lm(glucose ~ trait, data=my_data))$coefficients
  new_row=data.frame(trait=rownames(all_traits_family)[i],
                     est_beta_to_trait=model_beta[2,1],
                     est_trait_to_gluc=model_gluc[2,1],
                     final_effect=model_beta[2,1]*model_gluc[2,1],
                     pval_beta_to_trait=model_beta[2,4],
                     pval_trait_to_gluc=model_gluc[2,4])
  lm_glucose_family=rbind(lm_glucose_family,new_row)
}

lm_glucose_genus=data.frame()
for(i in 1:length(rownames(all_traits_genus))){
  print(rownames(all_traits_genus)[i])
  my_data=data.frame(glucose=glucose,trait=all_traits_genus[i,],beta=beta)
  model_beta=summary(lm(trait ~ beta, data=my_data))$coefficients
  model_gluc=summary(lm(glucose ~ trait, data=my_data))$coefficients
  new_row=data.frame(trait=rownames(all_traits_genus)[i],
                     est_beta_to_trait=model_beta[2,1],
                     est_trait_to_gluc=model_gluc[2,1],
                     final_effect=model_beta[2,1]*model_gluc[2,1],
                     pval_beta_to_trait=model_beta[2,4],
                     pval_trait_to_gluc=model_gluc[2,4])
  lm_glucose_genus=rbind(lm_glucose_genus,new_row)
}

lm_glucose_species=data.frame()
for(i in 1:length(rownames(all_traits_species))){
  print(rownames(all_traits_species)[i])
  my_data=data.frame(glucose=glucose,trait=all_traits_species[i,],beta=beta)
  model_beta=summary(lm(trait ~ beta, data=my_data))$coefficients
  model_gluc=summary(lm(glucose ~ trait, data=my_data))$coefficients
  new_row=data.frame(trait=rownames(all_traits_species)[i],
                     est_beta_to_trait=model_beta[2,1],
                     est_trait_to_gluc=model_gluc[2,1],
                     final_effect=model_beta[2,1]*model_gluc[2,1],
                     pval_beta_to_trait=model_beta[2,4],
                     pval_trait_to_gluc=model_gluc[2,4])
  lm_glucose_species=rbind(lm_glucose_species,new_row)
}

lm_glucose_all=rbind(
  lm_glucose_species,
  lm_glucose_genus,
  lm_glucose_family,
  lm_glucose_order,
  lm_glucose_class,
  lm_glucose_phylum,
  lm_glucose_cluster
)

lm_glucose_all$bonf_beta_to_trait = lm_glucose_all$pval_beta_to_trait*dim(lm_glucose_all)[1]
lm_glucose_all$bonf_trait_to_gluc = lm_glucose_all$pval_trait_to_gluc*dim(lm_glucose_all)[1]
lm_glucose_all$qvalue_beta_to_trait =p.adjust(lm_glucose_all$pval_beta_to_trait, method = "fdr")
lm_glucose_all$qvalue_trait_to_gluc =p.adjust(lm_glucose_all$pval_trait_to_gluc, method = "fdr")
lm_glucose_all <- lm_glucose_all %>% 
  mutate( Significance_beta_to_trait = case_when(
    qvalue_beta_to_trait < 0.1 & bonf_beta_to_trait< 0.05 ~ "Significant (Bonferroni)",
    qvalue_beta_to_trait < 0.1 & bonf_beta_to_trait >= 0.05 ~ "Significant (FDR)",
    TRUE ~ "Non-significant"))
lm_glucose_all <- lm_glucose_all %>% 
  mutate( Significance_trait_to_gluc = case_when(
    qvalue_trait_to_gluc < 0.1 & bonf_trait_to_gluc< 0.05 ~ "Significant (Bonferroni)",
    qvalue_trait_to_gluc < 0.1 & bonf_trait_to_gluc >= 0.05 ~ "Significant (FDR)",
    TRUE ~ "Non-significant"))
lm_glucose_all <- lm_glucose_all %>% 
  mutate( All_Sig = case_when(
    Significance_beta_to_trait != "Non-significant" & Significance_trait_to_gluc != "Non-significant" ~ "Significant",
    TRUE ~ "Non-significant"))

lm_glucose_all$Direction <- ifelse(lm_glucose_all$est_trait_to_gluc > 0, "Positive", "Negative")
lm_glucose_all = lm_glucose_all[order(lm_glucose_all$final_effect, decreasing = T),]
lm_glucose_all_filt=lm_glucose_all[lm_glucose_all$All_Sig=="Significant",]

##################################################################

correlations_glucose_all=rbind(
  correlations_glucose_species,
  correlations_glucose_genus,
  correlations_glucose_family,
  correlations_glucose_order,
  correlations_glucose_class,
  correlations_glucose_phylum,
  correlations_glucose_cluster
)
correlations_glucose_all$bonf_pvalue = correlations_glucose_all$pvalue*dim(correlations_glucose_all)[1]
correlations_glucose_all$qvalue =p.adjust(correlations_glucose_all$pvalue, method = "fdr")
correlations_glucose_all <- correlations_glucose_all %>% 
  mutate( Significance = case_when(
    qvalue < 0.1 & bonf_pvalue < 0.05 ~ "Significant (Bonferroni)",
    qvalue < 0.1 & bonf_pvalue >= 0.05 ~ "Significant (FDR)",
    TRUE ~ "Non-significant"))
correlations_glucose_all$Direction <- ifelse(correlations_glucose_all$estimate > 0, "Positive", "Negative")
correlations_glucose_all = correlations_glucose_all[order(correlations_glucose_all$pval, decreasing = T),]

correlations_glucose_alpha=cor.test(alpha,glucose, method = "pearson", exact=FALSE)
#summary(correlations_glucose_alpha)
correlations_glucose_beta=cor.test(beta,glucose, method = "pearson", exact=FALSE)
#summary(correlations_glucose_beta)

# Correlations weight
correlations_weight_species=data.frame()
for(i in 1:length(rownames(all_traits_species))){
  test_pearson_weight = cor.test(all_traits_species[i,],weight, method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_species)[i],estimate=test_pearson_weight$estimate,pval=-log(test_pearson_weight$p.value),pvalue=test_pearson_weight$p.value)
  correlations_weight_species=rbind(correlations_weight_species,new_row)
}

correlations_weight_cluster=data.frame()
for(i in 1:length(rownames(all_traits_cluster))){
  test_pearson_weight = cor.test(all_traits_cluster[i,],weight, method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_cluster)[i],estimate=test_pearson_weight$estimate,pval=-log(test_pearson_weight$p.value),pvalue=test_pearson_weight$p.value)
  correlations_weight_cluster=rbind(correlations_weight_cluster,new_row)
}

correlations_weight_phylum=data.frame()
for(i in 1:length(rownames(all_traits_phylum))){
  test_pearson_weight = cor.test(all_traits_phylum[i,],weight, method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_phylum)[i],estimate=test_pearson_weight$estimate,pval=-log(test_pearson_weight$p.value),pvalue=test_pearson_weight$p.value)
  correlations_weight_phylum=rbind(correlations_weight_phylum,new_row)
}

correlations_weight_class=data.frame()
for(i in 1:length(rownames(all_traits_phylum))){
  test_pearson_weight = cor.test(all_traits_class[i,],weight, method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_class)[i],estimate=test_pearson_weight$estimate,pval=-log(test_pearson_weight$p.value),pvalue=test_pearson_weight$p.value)
  correlations_weight_class=rbind(correlations_weight_class,new_row)
}

correlations_weight_order=data.frame()
for(i in 1:length(rownames(all_traits_order))){
  test_pearson_weight = cor.test(all_traits_order[i,],weight, method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_order)[i],estimate=test_pearson_weight$estimate,pval=-log(test_pearson_weight$p.value),pvalue=test_pearson_weight$p.value)
  correlations_weight_order=rbind(correlations_weight_order,new_row)
}

correlations_weight_family=data.frame()
for(i in 1:length(rownames(all_traits_family))){
  test_pearson_weight = cor.test(all_traits_family[i,],weight, method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_family)[i],estimate=test_pearson_weight$estimate,pval=-log(test_pearson_weight$p.value),pvalue=test_pearson_weight$p.value)
  correlations_weight_family=rbind(correlations_weight_family,new_row)
}

correlations_weight_genus=data.frame()
for(i in 1:length(rownames(all_traits_genus))){
  test_pearson_weight = cor.test(all_traits_genus[i,],weight, method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_genus)[i],estimate=test_pearson_weight$estimate,pval=-log(test_pearson_weight$p.value),pvalue=test_pearson_weight$p.value)
  correlations_weight_genus=rbind(correlations_weight_genus,new_row)
}

correlations_weight_all=rbind(
  correlations_weight_species,
  correlations_weight_genus,
  correlations_weight_family,
  correlations_weight_order,
  correlations_weight_class,
  correlations_weight_phylum,
  correlations_weight_cluster
)
correlations_weight_all$bonf_pvalue = correlations_weight_all$pvalue*dim(correlations_weight_all)[1]
correlations_weight_all$qvalue =p.adjust(correlations_weight_all$pvalue, method = "fdr")
correlations_weight_all <- correlations_weight_all %>% 
  mutate( Significance = case_when(
    qvalue < 0.1 & bonf_pvalue < 0.05 ~ "Significant (Bonferroni)",
    qvalue < 0.1 & bonf_pvalue >= 0.05 ~ "Significant (FDR)",
    TRUE ~ "Non-significant"))
correlations_weight_all$Direction <- ifelse(correlations_weight_all$estimate > 0, "Positive", "Negative")
correlations_weight_all = correlations_weight_all[order(correlations_weight_all$pval, decreasing = T),]


correlations_weight_alpha=cor.test(alpha,weight, method = "pearson", exact=FALSE)
#summary(correlations_weight_alpha)
correlations_weight_beta=cor.test(beta,weight, method = "pearson", exact=FALSE)
#summary(correlations_weight_beta)


correlations_prevotella_phylum=data.frame()
for(i in 1:length(rownames(all_traits_phylum))){
  test_pearson_prevotella_phylum = cor.test(all_traits_phylum[i,],all_traits_genus[which(rownames(all_traits_genus)=="g__Prevotella"),], method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_phylum)[i],estimate=test_pearson_prevotella_phylum$estimate,pval=-log(test_pearson_prevotella_phylum$p.value),pvalue=test_pearson_prevotella_phylum$p.value)
  correlations_prevotella_phylum=rbind(correlations_prevotella_phylum,new_row)
}

correlations_prevotella_genus=data.frame()
for(i in 1:length(rownames(all_traits_genus))){
  test_pearson_prevotella_genus = cor.test(all_traits_genus[i,],all_traits_genus[which(rownames(all_traits_genus)=="g__Prevotella"),], method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_genus)[i],estimate=test_pearson_prevotella_genus$estimate,pval=-log(test_pearson_prevotella_genus$p.value),pvalue=test_pearson_prevotella_genus$p.value)
  correlations_prevotella_genus=rbind(correlations_prevotella_genus,new_row)
}

#correlations_butyribacter_genus=data.frame()
#for(i in 1:length(rownames(all_traits_genus))){
#  test_pearson_butyribacter_genus = cor.test(all_traits_genus[i,],all_traits_genus[which(rownames(all_traits_genus)=="g__Butyribacter"),], method = "pearson", exact=FALSE)
#  new_row=data.frame(trait=rownames(all_traits_genus)[i],estimate=test_pearson_butyribacter_genus$estimate,pval=-log(test_pearson_butyribacter_genus$p.value),pvalue=test_pearson_butyribacter_genus$p.value)
#  correlations_butyribacter_genus=rbind(correlations_butyribacter_genus,new_row)
#}

correlations_prevotella_cluster=data.frame()
for(i in 1:length(rownames(all_traits_cluster))){
  test_pearson_prevotella_cluster = cor.test(all_traits_cluster[i,],all_traits_genus[which(rownames(all_traits_genus)=="g__Prevotella"),], method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_cluster)[i],estimate=test_pearson_prevotella_cluster$estimate,pval=-log(test_pearson_prevotella_cluster$p.value),pvalue=test_pearson_prevotella_cluster$p.value)
  correlations_prevotella_cluster=rbind(correlations_prevotella_cluster,new_row)
}

correlations_prevotella_family=data.frame()
for(i in 1:length(rownames(all_traits_family))){
  test_pearson_prevotella_family = cor.test(all_traits_family[i,],all_traits_genus[which(rownames(all_traits_genus)=="g__Prevotella"),], method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_family)[i],estimate=test_pearson_prevotella_family$estimate,pval=-log(test_pearson_prevotella_family$p.value),pvalue=test_pearson_prevotella_family$p.value)
  correlations_prevotella_family=rbind(correlations_prevotella_family,new_row)
}

correlations_prevotella_class=data.frame()
for(i in 1:length(rownames(all_traits_class))){
  test_pearson_prevotella_class = cor.test(all_traits_class[i,],all_traits_genus[which(rownames(all_traits_genus)=="g__Prevotella"),], method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_class)[i],estimate=test_pearson_prevotella_class$estimate,pval=-log(test_pearson_prevotella_class$p.value),pvalue=test_pearson_prevotella_class$p.value)
  correlations_prevotella_class=rbind(correlations_prevotella_class,new_row)
}

correlations_prevotella_order=data.frame()
for(i in 1:length(rownames(all_traits_order))){
  test_pearson_prevotella_order = cor.test(all_traits_order[i,],all_traits_genus[which(rownames(all_traits_genus)=="g__Prevotella"),], method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_order)[i],estimate=test_pearson_prevotella_order$estimate,pval=-log(test_pearson_prevotella_order$p.value),pvalue=test_pearson_prevotella_order$p.value)
  correlations_prevotella_order=rbind(correlations_prevotella_order,new_row)
}

correlations_prevotella_species=data.frame()
for(i in 1:length(rownames(all_traits_species))){
  test_pearson_prevotella_species = cor.test(all_traits_species[i,],all_traits_genus[which(rownames(all_traits_genus)=="g__Prevotella"),], method = "pearson", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_species)[i],estimate=test_pearson_prevotella_species$estimate,pval=-log(test_pearson_prevotella_species$p.value),pvalue=test_pearson_prevotella_species$p.value)
  correlations_prevotella_species=rbind(correlations_prevotella_species,new_row)
}

correlations_prevotella_all=rbind(
  correlations_prevotella_species,
  correlations_prevotella_genus,
  correlations_prevotella_family,
  correlations_prevotella_order,
  correlations_prevotella_class,
  correlations_prevotella_phylum,
  correlations_prevotella_cluster
)
correlations_prevotella_all$bonf_pvalue = correlations_prevotella_all$pvalue*dim(correlations_prevotella_all)[1]
correlations_prevotella_all$qvalue =p.adjust(correlations_prevotella_all$pvalue, method = "fdr")
correlations_prevotella_all <- correlations_prevotella_all %>% 
  mutate( Significance = case_when(
    qvalue < 0.1 & bonf_pvalue < 0.05 ~ "Significant (Bonferroni)",
    qvalue < 0.1 & bonf_pvalue >= 0.05 ~ "Significant (FDR)",
    TRUE ~ "Non-significant"))
correlations_prevotella_all$Direction <- ifelse(correlations_prevotella_all$estimate > 0, "Positive", "Negative")
correlations_prevotella_all = correlations_prevotella_all[order(correlations_prevotella_all$pval, decreasing = T),]


# Read Clusters
clusters <- read.csv(paste0("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Func_Mat_Net/Cluster_Analyses/Cluster_UMAP_HDBSCAN_",subset_id,".csv"))
colnames(clusters)[1]="Genome"
ALL_RES <- read.csv(paste0("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Func_Mat_Net/Matrix_Processing/matrices_residues/",subset_id,".csv"))
clusters$Genome=ALL_RES$Genome_id
clusters$Cluster=gsub("-1","noise",clusters$Cluster)
clusters$Cluster=paste("cluster",clusters$Cluster,sep="__")
clusters$Species=taxonomy$Species[match(clusters$Genome,taxonomy$Genome)]
clusters$Genus=taxonomy$Genus[match(clusters$Genome,taxonomy$Genome)]
clusters$Family=taxonomy$Family[match(clusters$Genome,taxonomy$Genome)]
clusters$Phylum=taxonomy$Phylum[match(clusters$Genome,taxonomy$Genome)]
stats_clusters=count(clusters,Cluster)


# Read Enterotypes
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Func_Mat_Net/Enterotypes/enterotypes.RData")
number=1
if (subset_id=="MI"){
  number=number-2
} else {
  if (subset_id=="NY"){
    number=number-1
  }
}
ent_samples=enterotypes_samples_objs[[number]]
stats_ent_samples=count(ent_samples,Enterotype,study,sex)
#stats_ent_samples$enterotype=paste0("enterotype__",stats_ent_samples$enterotype)
stats_ent_samples$sex[stats_ent_samples$sex=="F"]="Female"
stats_ent_samples$sex[stats_ent_samples$sex=="M"]="Male"
ent_samples$alpha=alpha[match(ent_samples$sample,names(alpha))]
ent_samples$beta=beta[match(ent_samples$sample,names(beta))]
ent_samples$glucose=glucose[match(ent_samples$sample,names(glucose))]
ent_samples$weight=weight[match(ent_samples$sample,names(weight))]

# Read centrality values
load("/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Func_Mat_Net/Network/networks.RData")
#number=1
#if (subset_id=="MI"){
#  number=number-2
#} else {
#  if (subset_id=="NY"){
#    number=number-1
#  }
#}
centrality=centralities#_objs[[number]]
hubs=hubs_clusters#objs[[number]]
#hubs$hub="Hubs"
#hubs$Cluster=clusters$Cluster[match(hubs$Species,clusters$Species)]
#stats_hubs=stats_hubs_objs[[number]]

## Read Heritability
##heritability="/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Func_Mat_Net/Heritability/heritability.RData"
#heritability="/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Func_Mat_Net/Heritability/heritability.RData"
#load(heritability)
#
####
#
#herit_NY= all_VCs_full[grepl("NY",all_VCs_full$subset_id),]
#herit_MI= all_VCs_full[grepl("MI",all_VCs_full$subset_id),]
#
#herit_NY=herit_NY[herit_NY$trait %in% intersect(herit_NY$trait,herit_MI$trait),]
#herit_MI=herit_MI[herit_MI$trait %in% intersect(herit_NY$trait,herit_MI$trait),]
#
#herit_NY=herit_NY[herit_MI$trait,]
#
#cor.test(herit_NY$var_Ad,herit_MI$var_Ad,method="pearson")
#
####
#
#heritability="/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_16S_Full_Harm/Heritability/heritability.RData"
#heritability="/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Shallow_Full_Harm/Heritability/heritability.RData"
#heritability="/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_16S_Sample_Harm/Heritability/heritability.RData"
#heritability="/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Shallow_Sample_Harm/Heritability/heritability.RData"
heritability="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Func_Mat_Net/Heritability/heritability.RData"
load(heritability)

filtered_VCs=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),]
filtered_VCs=filtered_VCs[grepl("EC4",filtered_VCs$rank),]

ec_numbers <- sapply(strsplit(filtered_VCs$trait, "_"),"[",3)
ec_components <- data.frame(do.call(rbind, strsplit(ec_numbers, "\\.")))
colnames(ec_components) <- c("L1", "L2", "L3", "L4")
ec_components$Level1 <- ec_components$L1
ec_components$Level2 <- with(ec_components, paste(Level1, L2, sep = "."))
ec_components$Level3 <- with(ec_components, paste(Level2, L3, sep = "."))
ec_components$Level4 <- with(ec_components, paste(Level3, L4, sep = "."))
ec_components=ec_components[,5:8]
rownames(ec_components)=paste0("EC4__",ec_components$Level4)


#subset_id="ALL" ### Remove
#sig_f=0.1
#sig_b=0.05
#all_VCs_full=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),] ### Remove
#
#all_VCs_full$bonf_pvalue_DGE = all_VCs_full$pvalue_DGE*dim(all_VCs_full)[1]
#all_VCs_full$qvalue_DGE = qvalue(all_VCs_full$pvalue_DGE)$qvalue
#all_VCs_full$qp_diff = all_VCs_full$qvalue_DGE - all_VCs_full$pvalue_DGE
#all_VCs_full <- all_VCs_full %>% 
#  mutate( Significance = case_when(
#    qvalue_DGE < sig_f & bonf_pvalue_DGE < sig_b ~ "Significant (Bonferroni)",
#    qvalue_DGE < sig_f & bonf_pvalue_DGE >= sig_b ~ "Significant (FDR)",
#    TRUE ~ "Non-significant traits"))
#all_VCs_full = all_VCs_full[order(all_VCs_full$logP, decreasing = T),]
#all_VCs_full$qvalue_DGE = round(all_VCs_full$qvalue_DGE, digits = 2)

# Include columns
# Significance (FDR<10%)
filtered_VCs=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),]
filtered_VCs=filtered_VCs[grepl("Species",filtered_VCs$rank),]
filtered_VCs <- filtered_VCs %>% 
  mutate( All_Sig = case_when(
    Significance ==  "Significant (Bonferroni)" | Significance == "Significant (FDR)" ~ "Heritable",
    TRUE ~ "Non-heritable"))
# Abundane and Prevalence
filtered_VCs$Abundance=ab_vs_prev$abundance[match(filtered_VCs$trait,ab_vs_prev$trait)]
filtered_VCs$Prevalence=ab_vs_prev$prevalence[match(filtered_VCs$trait,ab_vs_prev$trait)]
# Genome Size
filtered_VCs$Genome_size=as.numeric(catalogue_stats$Size[match(filtered_VCs$trait,catalogue_stats$Species)])
filtered_VCs$Number_of_genes=as.numeric(catalogue_stats$Number_of_genes[match(filtered_VCs$trait,catalogue_stats$Species)])
filtered_VCs$Mean_gene_size=as.numeric(catalogue_stats$mean[match(filtered_VCs$trait,catalogue_stats$Species)])
# Correlations with beta
filtered_VCs$correlations_beta=correlations_beta_species$correlations_beta[match(filtered_VCs$trait,correlations_beta_species$trait)]
# Correlations with alpha
filtered_VCs$correlations_alpha=correlations_alpha_species$correlations_alpha[match(filtered_VCs$trait,correlations_alpha_species$trait)]
# Centralities
filtered_VCs$Degree=as.numeric(centralities$degree[match(filtered_VCs$trait,centralities$Traits)])
filtered_VCs$Betweenness=as.numeric(centralities$betweenness[match(filtered_VCs$trait,centralities$Traits)])
filtered_VCs$Closeness=as.numeric(centralities$closeness[match(filtered_VCs$trait,centralities$Traits)])
filtered_VCs$Eigenvector=as.numeric(centralities$eigenvector[match(filtered_VCs$trait,centralities$Traits)])
# Hubs
filtered_VCs$Hubs=hubs_clusters$hub_type[match(filtered_VCs$trait,hubs_clusters$Traits)]
#filtered_VCs$Hub=hubs$hub[match(filtered_VCs$trait,hubs$Species)]
#filtered_VCs$Hub[is.na(filtered_VCs$Hub)]="Others"
# Taxonomy
filtered_VCs$Genome=taxonomy$Genome[match(filtered_VCs$trait,taxonomy$Species)]
filtered_VCs$Phylum=taxonomy$Phylum[match(filtered_VCs$trait,taxonomy$Species)]
filtered_VCs$Class=taxonomy$Class[match(filtered_VCs$trait,taxonomy$Species)]
filtered_VCs$Order=taxonomy$Order[match(filtered_VCs$trait,taxonomy$Species)]
filtered_VCs$Family=taxonomy$Family[match(filtered_VCs$trait,taxonomy$Species)]
filtered_VCs$Genus=taxonomy$Genus[match(filtered_VCs$trait,taxonomy$Species)]
# Cluster
filtered_VCs$Cluster=clusters$Cluster[match(filtered_VCs$trait,clusters$Traits)]
## Enterotype
#filtered_VCs$Enterotype=enterotypes$Enterotype[match(filtered_VCs$trait,enterotypes$Species)]
# Heritability
filtered_VCs$Heritability=as.numeric(filtered_VCs$var_Ad)*100
filtered_VCs$Maternal=as.numeric(filtered_VCs$var_M)*100
filtered_VCs$Cohousing=as.numeric(filtered_VCs$var_C)*100

save(all_VCs_full,
     filtered_VCs,
     hubs_clusters,
     centrality,
     catalogue_stats,
     correlations_alpha_all,
     correlations_beta_all,
     correlations_cluster3_clusters,
     correlations_glucose_all,
     correlations_prevotella_all,
     correlations_weight_all,
     ent_samples,
     taxonomy,
     topology_clusters,
     topology_full_network,
     check_network_res,
     lm_glucose_all,
     file="/users/abaud/fmorillo/paper_figures/full_tax_heritability_analyses.RData")


#no_bacteroidota=filtered_VCs[filtered_VCs$Phylum!="p__Bacteroidota",]
#stats_no_bacteroidota=count(no_bacteroidota,Enterotype,Phylum,All_Sig)

##################################################################### PLOTS

pdf("/users/abaud/fmorillo/paper_figures/heritability_comparisons_ALL.pdf")

##################################################################### ALL EFFECTS

## Boxplot effects

matern=filtered_VCs
matern$Effect=matern$Maternal
matern$Effect_Type="Maternal Effects"
cohous=filtered_VCs
cohous$Effect=cohous$Cohousing
cohous$Effect_Type="Co-housing Effects"
herit=filtered_VCs
herit$Effect=herit$Heritability
herit$Effect_Type="Heritability"
filtered_VCs_effects=rbind(matern,cohous,herit)

#order=c("Co-housing Effects","Maternal Effects","Heritability")
myplot<-ggplot(filtered_VCs_effects, aes(x=reorder(Effect_Type,Effect,FUN = median), y=Effect)) +
  geom_boxplot() +
  xlab("") +
  ylab("Effect Size (%)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Co-housing Effects","Heritability"),c("Co-housing Effects","Heritability"),c("Maternal Effects","Heritability")),
              map_signif_level = TRUE)
print(myplot)

##################################################################### HERITABILITY VALUES

## By rank - Only taxonomic
only_taxonomic=all_VCs_full[(all_VCs_full$subset_id==subset_id) & (all_VCs_full$rank != "Clusters" & all_VCs_full$rank != "Beta Diversity" & all_VCs_full$rank != "Alpha Diversity"),]
#only_taxonomic=all_VCs_full[(all_VCs_full$subset_id==subset_id) & (all_VCs_full$rank == "Species" | all_VCs_full$rank == "Genus"),]
only_taxonomic$Heritability=as.numeric(only_taxonomic$var_Ad)*100
order=c("Phylum", "Class", "Order", "Family", "Genus", "Species")
#order=c("Genus", "Species")
myplot<-ggplot(only_taxonomic, aes(x=factor(rank, levels = order), y=Heritability)) +
  geom_boxplot() +
  xlab("") +
  ylab("Heritability (%)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60))# +
#geom_signif(comparisons = list(c("Species","Genus")),
#          map_signif_level = TRUE)
print(myplot)

## By rank - With clusters
only_taxonomic=all_VCs_full[(all_VCs_full$subset_id==subset_id) & (all_VCs_full$rank != "Beta Diversity" & all_VCs_full$rank != "Alpha Diversity"),]
only_taxonomic$rank[only_taxonomic$rank=="Clusters"]="Guilds"
only_taxonomic$Heritability=as.numeric(only_taxonomic$var_Ad)*100
order=c("Phylum", "Class", "Order", "Family", "Genus", "Species", "Guilds")
myplot<-ggplot(only_taxonomic, aes(x=factor(rank, levels = order), y=Heritability)) +
  geom_boxplot() +
  xlab("") +
  ylab("Heritability (%)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60))# +
#geom_signif(comparisons = list(c("Species","Guilds"),c("Genus","Guilds")),
#            map_signif_level = TRUE)
print(myplot)

# Phylum
myplot<-ggplot(filtered_VCs, aes(x=reorder(Phylum,Heritability,FUN = median), y=Heritability)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species heritability (%)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("p__Firmicutes_A","p__Bacteroidota")),
              map_signif_level = TRUE)
print(myplot)

# Family
myplot<-ggplot(filtered_VCs[filtered_VCs$Phylum=="p__Bacteroidota",], aes(x=reorder(Family,Heritability,FUN = median), y=Heritability)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species heritability (%)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("f__Bacteroidaceae","f__Muribaculaceae")),
              map_signif_level = TRUE)
print(myplot)

# Genus
myplot<-ggplot(filtered_VCs[filtered_VCs$Family=="f__Bacteroidaceae",], aes(x=reorder(Genus,Heritability,FUN = median), y=Heritability)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species heritability (%)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("g__Bacteroides","g__Prevotella")),
              map_signif_level = TRUE)
print(myplot)

# Species heritability by genus (only the ones with more than 10 species)
# Step 1: Count the number of occurrences of each genus
genus_counts <- filtered_VCs %>%
  group_by(Genus) %>%
  summarize(count = n())
# Step 2: Filter out genera with fewer than 10 occurrences
filtered_genus <- genus_counts %>%
  filter(count >= 10) %>%
  pull(Genus)
# Step 3: Filter the original data to include only genera with at least 10 occurrences
filtered_VCs_min_10 <- filtered_VCs %>%
  filter(Genus %in% filtered_genus)
# Genus
myplot<-ggplot(filtered_VCs_min_10, aes(x=reorder(Genus,Heritability,FUN = median), y=Heritability)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species heritability (%)") +
  #ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  #scale_x_discrete(guide = guide_axis(angle = 60))+
  coord_flip()
#geom_signif(comparisons = list(c("f__Bacteroidaceae","f__Muribaculaceae")),
#            map_signif_level = TRUE)
print(myplot)

# Clusters
myplot<-ggplot(filtered_VCs, aes(x=reorder(Cluster,Heritability,FUN = median), y=Heritability)) +
  geom_boxplot() +
  xlab("") +
  ylab("Cluster heritability (%)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60))
#geom_signif(comparisons = list(c("cluster__3","cluster__11"),c("cluster__11","cluster__10"),c("cluster__10","cluster__13")),
#            map_signif_level = TRUE)
print(myplot)


# Hubs
order=c("Non-hub","Cluster Hub","Global Hub","Global and Cluster Hub")
myplot<-ggplot(filtered_VCs, aes(x=factor(Hubs, levels = order), y=Heritability)) +
  geom_boxplot() +
  xlab("") +
  ylab("Heritability (%)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Non-hub","Cluster Hub"),
                                 c("Cluster Hub","Global Hub"),
                                 c("Global Hub","Global and Cluster Hub")
  ),
  map_signif_level = TRUE)
print(myplot)

##################################################################### SIGNIFICANT HERITABILITY

# Porportion of significant by phylum
stats <- filtered_VCs %>%
  group_by(Phylum) %>%
  summarise(count = n()) %>%
  ungroup()
prop_data <- filtered_VCs %>%
  group_by(Phylum, All_Sig) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()
prop_data$proportion=round(100*prop_data$proportion,digit=2)
filtered=prop_data[grepl("Non-heritable",prop_data$All_Sig),]
filtered= filtered[filtered$proportion!=100,]
stats = stats[stats$Phylum %in% filtered$Phylum,]
stats = stats[order(stats$count, decreasing = F),]
myplot <- ggplot(prop_data[prop_data$Phylum %in% stats$Phylum,], aes(x = factor(Phylum, levels = stats$Phylum), y = count, fill = All_Sig)) +
  geom_bar(stat = "identity", position = "stack", linewidth = 2, col = "white", width = 1) +
  #scale_fill_manual(values = colours_phylum) +
  labs(fill = "Species") +
  #ggtitle("Proportion of heritable species (No Bacteroidota)") +
  xlab("") +
  ylab("Number of species") +
  #scale_x_discrete(guide = guide_axis(angle = 60)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),)+
  coord_flip()
print(myplot)

# Porportion of significant by family
stats <- filtered_VCs %>%
  group_by(Family) %>%
  summarise(count = n()) %>%
  ungroup()
prop_data <- filtered_VCs %>%
  group_by(Family, All_Sig) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()
prop_data$proportion=round(100*prop_data$proportion,digit=2)
filtered=prop_data[grepl("Non-heritable",prop_data$All_Sig),]
filtered= filtered[filtered$proportion!=100,]
stats = stats[stats$Family %in% filtered$Family,]
stats = stats[order(stats$count, decreasing = F),]
myplot <- ggplot(prop_data[prop_data$Family %in% stats$Family,], aes(x = factor(Family, levels = stats$Family), y = count, fill = All_Sig)) +
  geom_bar(stat = "identity", position = "stack", linewidth = 2, col = "white", width = 1) +
  #scale_fill_manual(values = colours_phylum) +
  labs(fill = "Species") +
  #ggtitle("Proportion of heritable species (No Bacteroidota)") +
  xlab("") +
  ylab("Number of species") +
  #scale_x_discrete(guide = guide_axis(angle = 60)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),)+
  coord_flip()
print(myplot)

# Porportion of significant by genus
stats <- filtered_VCs %>%
  group_by(Genus) %>%
  summarise(count = n()) %>%
  ungroup()
prop_data <- filtered_VCs %>%
  group_by(Genus, All_Sig) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()
prop_data$proportion=round(100*prop_data$proportion,digit=2)
filtered=prop_data[grepl("Non-heritable",prop_data$All_Sig),]
filtered= filtered[filtered$proportion!=100,]
stats = stats[stats$Genus %in% filtered$Genus,]
stats = stats[order(stats$count, decreasing = F),]
myplot <- ggplot(prop_data[prop_data$Genus %in% stats$Genus,], aes(x = factor(Genus, levels = stats$Genus), y = count, fill = All_Sig)) +
  geom_bar(stat = "identity", position = "stack", linewidth = 2, col = "white", width = 1) +
  #scale_fill_manual(values = colours_phylum) +
  labs(fill = "Species") +
  #ggtitle("Proportion of heritable species (No Bacteroidota)") +
  xlab("") +
  ylab("Number of species") +
  #scale_x_discrete(guide = guide_axis(angle = 60)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),)+
  coord_flip()
print(myplot)

# Porportion of significant in the clusters
prop_data <- filtered_VCs %>%
  group_by(Cluster, All_Sig) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()
prop_data$proportion=round(100*prop_data$proportion,digit=2)
order=prop_data[grepl("Non-heritable",prop_data$All_Sig),]
order = order[order(order$proportion, decreasing = F),]
order=order$Cluster
myplot <- ggplot(prop_data, aes(x = factor(Cluster, levels = order), y = proportion, fill = All_Sig)) +
  geom_bar(stat = "identity", position = "stack", linewidth = 2, col = "white",  width = 1) +
  #scale_fill_manual(values = colours_phylum) +
  labs(fill = "Species") +
  ggtitle("Proportion of species with significant heritability (FDR<10%)") +
  xlab("") +
  ylab("Proportion (%)") +
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),)
print(myplot)

##################################################################### ENTEROTYPES

# Enterotypes statistics (total)
myplot <- ggplot(stats_ent_samples, aes(x = sex, y = n, fill = Enterotype)) +
  geom_bar(stat = "identity", position = "stack", linewidth = 2, col = "white",  width = 1) +
  #scale_fill_manual(values = colours_phylum) +
  labs(fill = "Enterotype") +
  facet_wrap(~ study) +
  #ggtitle("Proportion of species with significant heritability (FDR<10%)") +
  xlab("") +
  ylab("Number of samples") +
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),)
print(myplot)

### Enterotypes relative abundances per phylum
# Make into relative abundance
phylum_matrix <- apply(filtered_prev_new_phylum, 2, function(i) i/sum(i))
colSums(phylum_matrix)[1:5]
# Convert phylum_matrix to a dataframe
phylum_df <- as.data.frame(phylum_matrix)
phylum_df <- tibble::rownames_to_column(phylum_df, var = "phylum")
# Gather the phylum matrix to long format
phylum_long <- phylum_df %>%
  pivot_longer(cols = -phylum, names_to = "sample", values_to = "Abundance")
# Merge with ent_samples to get enterotype information
merged_df <- phylum_long %>%
  inner_join(ent_samples, by = "sample")
# Calculate mean abundance of each phylum per enterotype
mean_abundance <- merged_df %>%
  group_by(Enterotype, phylum) %>%
  summarise(MeanAbundance = mean(Abundance)) %>%
  ungroup()
mean_abundance$MeanAbundance=100*mean_abundance$MeanAbundance
###
reshaped_data_phylum <- mean_abundance %>%
  spread(key = Enterotype, value = MeanAbundance) %>%
  mutate(log_fold_change = log2(`enterotype__2` / `enterotype__1`))
reshaped_data_phylum$diff= reshaped_data_phylum$`enterotype__2`-reshaped_data_phylum$`enterotype__1`
###
# Select top 10 phylum per enterotype
top10_phylum <- mean_abundance %>%
  group_by(Enterotype) %>%
  top_n(nrow(filtered_prev_new_phylum), MeanAbundance) %>%
  ungroup()
top10_phylum$group=paste(top10_phylum$Enterotype,top10_phylum$phylum,sep="_")
# Create a new column to indicate the direction for the bidirectional plot
top10_phylum <- top10_phylum %>%
  mutate(Direction = ifelse(Enterotype == unique(Enterotype)[1], -MeanAbundance, MeanAbundance))
# Order phylum by decreasing mean abundance within each enterotype
top10_phylum <- top10_phylum %>%
  group_by(Enterotype) %>%
  arrange(Enterotype, desc(MeanAbundance)) %>%
  mutate(phylum = factor(phylum, levels = unique(phylum))) %>%
  ungroup()
# Creating the bidirectional barplot with the adjusted theme
max=max(top10_phylum$MeanAbundance)
myplot=ggplot(top10_phylum, aes(x = phylum, y = MeanAbundance, fill = Enterotype)) +
  geom_bar(stat = "identity") +
  #coord_flip() +
  scale_y_continuous(labels = abs) +
  ylim(0,max+1) +
  labs(y = "Mean relative abundance (%)", x = NULL, 
       title = "Phylum by enterotype") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    strip.text = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 15),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 12),
    legend.text=element_text(size=12),
    legend.title=element_text(size=15),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  facet_grid(Enterotype ~ ., scales = "free_y", switch = "y")+
  scale_x_discrete(guide = guide_axis(angle = 60))
plot(myplot)

### Enterotypes relative abundances per class
# Make into relative abundance
class_matrix <- apply(filtered_prev_new_class, 2, function(i) i/sum(i))
colSums(class_matrix)[1:5]
# Convert class_matrix to a dataframe
class_df <- as.data.frame(class_matrix)
class_df <- tibble::rownames_to_column(class_df, var = "Class")
# Gather the class matrix to long format
class_long <- class_df %>%
  pivot_longer(cols = -Class, names_to = "sample", values_to = "Abundance")
# Merge with ent_samples to get enterotype information
merged_df <- class_long %>%
  inner_join(ent_samples, by = "sample")
# Calculate mean abundance of each class per enterotype
mean_abundance <- merged_df %>%
  group_by(Enterotype, Class) %>%
  summarise(MeanAbundance = mean(Abundance)) %>%
  ungroup()
mean_abundance$MeanAbundance=100*mean_abundance$MeanAbundance
###
reshaped_data_class <- mean_abundance %>%
  spread(key = Enterotype, value = MeanAbundance) %>%
  mutate(log_fold_change = log2(`enterotype__2` / `enterotype__1`))
#reshaped_data_class$class=taxonomy$class[match(reshaped_data_class$class,taxonomy$class)]
#reshaped_data_class$Class=taxonomy$Class[match(reshaped_data_class$Class,taxonomy$Class)]
reshaped_data_class$Class=taxonomy$Class[match(reshaped_data_class$Class,taxonomy$Class)]
reshaped_data_class$diff= reshaped_data_class$`enterotype__2`-reshaped_data_class$`enterotype__1`

###
# Select top 10 class per enterotype
top10_class <- mean_abundance %>%
  group_by(Enterotype) %>%
  top_n(20, MeanAbundance) %>%
  ungroup()
top10_class$group=paste(top10_class$Enterotype,top10_class$class,sep="_")
## Create a new column to indicate the direction for the bidirectional plot
#top10_class <- top10_class %>%
#  mutate(Direction = ifelse(Enterotype == unique(Enterotype)[1], -MeanAbundance, MeanAbundance))
# Order class by decreasing mean abundance within each enterotype
top10_class$Phylum=taxonomy$Phylum[match(top10_class$Class,taxonomy$Class)]
top10_class <- top10_class %>%
  group_by(Enterotype) %>%
  arrange(Enterotype, desc(Phylum)) %>%
  arrange(Enterotype, desc(MeanAbundance)) %>%
  mutate(Class = factor(Class, levels = unique(Class))) %>%
  ungroup()
# Creating the bidirectional barplot with the adjusted theme
top10_class$Enterotype=gsub("enterotype__","Enterotype ",top10_class$Enterotype)
#top10_class$LBL[top10_class$class == "f__Bacteroidaceae" & top10_class$Enterotype == "enterotype__2"] <- "*"
#top10_class$LBL[top10_class$class == "f__Oscillospiraceae" & top10_class$Enterotype == "enterotype__1"] <- "*"
#top10_class$LBL[is.na(top10_class$LBL)] <- ""
max=max(top10_class$MeanAbundance)
myplot = ggplot(top10_class, aes(x = Class, y = MeanAbundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  #coord_flip() +
  scale_y_continuous(labels = abs) +
  labs(y = "Mean relative abundance (%)", x = NULL, 
       title = "Top 20 classes by enterotype") +
  theme_minimal() +
  ylim(0, max(top10_class$MeanAbundance) + 1) +  # Adjusted to use max from the data
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    #legend.position = "bottom",
    axis.line = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  facet_grid(Enterotype ~ ., scales = "free_y", switch = "y") +
  scale_x_discrete(guide = guide_axis(angle = 60))
#geom_text(data = top10_class,
#          aes(x = class, y = MeanAbundance,
#              label = format(LBL, nsmall = 0, digits = 1, scientific = FALSE)),
#          color = "black", vjust = -0.5, hjust = 0.5, size = 7)
plot(myplot)

### Enterotypes relative abundances per genus
# Make into relative abundance
genus_matrix <- apply(filtered_prev_new_genus, 2, function(i) i/sum(i))
colSums(genus_matrix)[1:5]
# Convert genus_matrix to a dataframe
genus_df <- as.data.frame(genus_matrix)
genus_df <- tibble::rownames_to_column(genus_df, var = "Genus")
# Gather the genus matrix to long format
genus_long <- genus_df %>%
  pivot_longer(cols = -Genus, names_to = "sample", values_to = "Abundance")
# Merge with ent_samples to get enterotype information
merged_df <- genus_long %>%
  inner_join(ent_samples, by = "sample")
# Calculate mean abundance of each genus per enterotype
mean_abundance <- merged_df %>%
  group_by(Enterotype, Genus) %>%
  summarise(MeanAbundance = mean(Abundance)) %>%
  ungroup()
mean_abundance$MeanAbundance=100*mean_abundance$MeanAbundance
###
reshaped_data_genus <- mean_abundance %>%
  spread(key = Enterotype, value = MeanAbundance) %>%
  mutate(log_fold_change = log2(`enterotype__2` / `enterotype__1`))
reshaped_data_genus$diff=reshaped_data_genus$enterotype__1 - reshaped_data_genus$enterotype__2
reshaped_data_genus$Family=taxonomy$Family[match(reshaped_data_genus$Genus,taxonomy$Genus)]
reshaped_data_genus$Class=taxonomy$Class[match(reshaped_data_genus$Genus,taxonomy$Genus)]
reshaped_data_genus$Phylum=taxonomy$Phylum[match(reshaped_data_genus$Genus,taxonomy$Genus)]
###
# Select top 10 genera per enterotype
top10_genera <- mean_abundance %>%
  group_by(Enterotype) %>%
  top_n(20, MeanAbundance) %>%
  ungroup()
top10_genera$group=paste(top10_genera$Enterotype,top10_genera$Genus,sep="_")
## Create a new column to indicate the direction for the bidirectional plot
#top10_genera <- top10_genera %>%
#  mutate(Direction = ifelse(Enterotype == unique(Enterotype)[1], -MeanAbundance, MeanAbundance))
# Order genera by decreasing mean abundance within each enterotype
top10_genera$Phylum=taxonomy$Phylum[match(top10_genera$Genus,taxonomy$Genus)]
top10_genera <- top10_genera %>%
  group_by(Enterotype) %>%
  arrange(Enterotype, desc(Phylum)) %>%
  arrange(Enterotype, desc(MeanAbundance)) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus))) %>%
  ungroup()
# Creating the bidirectional barplot with the adjusted theme
top10_genera$Enterotype=gsub("enterotype__","Enterotype ",top10_genera$Enterotype)
top10_genera$LBL[top10_genera$Genus == "g__Prevotella" & top10_genera$Enterotype == "enterotype__2"] <- "*"
top10_genera$LBL[top10_genera$Genus == "g__UBA3282" & top10_genera$Enterotype == "enterotype__1"] <- "*"
top10_genera$LBL[is.na(top10_genera$LBL)] <- ""
max=max(top10_genera$MeanAbundance)
myplot = ggplot(top10_genera, aes(x = Genus, y = MeanAbundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  #coord_flip() +
  scale_y_continuous(labels = abs) +
  labs(y = "Mean relative abundance (%)", x = NULL, 
       title = "Top 20 genera by enterotype") +
  theme_minimal() +
  ylim(0, max(top10_genera$MeanAbundance) + 1) +  # Adjusted to use max from the data
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    #legend.position = "bottom",
    axis.line = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  facet_grid(Enterotype ~ ., scales = "free_y", switch = "y") +
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_text(data = top10_genera,
            aes(x = Genus, y = MeanAbundance,
                label = format(LBL, nsmall = 0, digits = 1, scientific = FALSE)),
            color = "black", vjust = -0.5, hjust = 0.5, size = 7)
plot(myplot)

### Enterotypes relative abundances per cluster
# Make into relative abundance
cluster_matrix <- apply(filtered_prev_new_cluster, 2, function(i) i/sum(i))
colSums(cluster_matrix)[1:5]
# Convert cluster_matrix to a dataframe
cluster_df <- as.data.frame(cluster_matrix)
cluster_df <- tibble::rownames_to_column(cluster_df, var = "Cluster")
# Gather the custer matrix to long format
cluster_long <- cluster_df %>%
  pivot_longer(cols = -Cluster, names_to = "sample", values_to = "Abundance")
# Merge with ent_samples to get enterotype information
merged_df <- cluster_long %>%
  inner_join(ent_samples, by = "sample")
# Calculate mean abundance of each cluster per enterotype
mean_abundance <- merged_df %>%
  group_by(Enterotype, Cluster) %>%
  summarise(MeanAbundance = mean(Abundance)) %>%
  ungroup()
mean_abundance$MeanAbundance=100*mean_abundance$MeanAbundance
###
reshaped_data_clusters <- mean_abundance %>%
  spread(key = Enterotype, value = MeanAbundance) %>%
  mutate(log_fold_change = log2(`enterotype__2` / `enterotype__1`))
reshaped_data_clusters$diff=reshaped_data_clusters$`enterotype__2`-reshaped_data_clusters$`enterotype__1`
###
top10_clusters <- mean_abundance %>%
  group_by(Enterotype) %>%
  top_n(nrow(filtered_prev_new_cluster), MeanAbundance) %>%
  ungroup()
top10_clusters$group=paste(top10_clusters$Enterotype,top10_clusters$Cluster,sep="_")
# Create a new column to indicate the direction for the bidirectional plot
top10_clusters <- top10_clusters %>%
  mutate(Direction = ifelse(Enterotype == unique(Enterotype)[1], -MeanAbundance, MeanAbundance))
# Order genera by decreasing mean abundance within each enterotype
top10_clusters <- top10_clusters %>%
  group_by(Enterotype) %>%
  arrange(Enterotype, desc(MeanAbundance)) %>%
  mutate(Cluster = factor(Cluster, levels = unique(Cluster))) %>%
  ungroup()
# Creating the bidirectional barplot with the adjusted theme
max=max(top10_clusters$MeanAbundance)
myplot=ggplot(top10_clusters, aes(x = Cluster, y = MeanAbundance, fill = Enterotype)) +
  geom_bar(stat = "identity") +
  #coord_flip() +
  scale_y_continuous(labels = abs) +
  ylim(0,max+1) +
  labs(y = "Mean Mean relative abundance (%)", x = NULL, 
       title = "Cluster by enterotype") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    strip.text = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 15),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 12),
    legend.text=element_text(size=12),
    legend.title=element_text(size=15),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  facet_grid(Enterotype ~ ., scales = "free_y", switch = "y")+
  scale_x_discrete(guide = guide_axis(angle = 60))
plot(myplot)

## Boxplots

order=c("enterotype__1","enterotype__2")
myplot<-ggplot(ent_samples, aes(x=factor(Enterotype, levels = order), y=alpha)) +
  geom_boxplot() +
  xlab("") +
  ylab("Alpha diversity (Rao's quadratic index)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("enterotype__1","enterotype__2")),
              map_signif_level = TRUE)
print(myplot)

order=c("enterotype__1","enterotype__2")
myplot<-ggplot(ent_samples, aes(x=factor(Enterotype, levels = order), y=beta)) +
  geom_boxplot() +
  xlab("") +
  ylab("Beta diversity (Aitchison distances) PCoA1") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("enterotype__1","enterotype__2")),
              map_signif_level = TRUE)
print(myplot)

order=c("enterotype__1","enterotype__2")
myplot<-ggplot(ent_samples, aes(x=factor(Enterotype, levels = order), y=glucose)) +
  geom_boxplot() +
  xlab("") +
  ylab("Glucose levels (residuals)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("enterotype__1","enterotype__2")),
              map_signif_level = TRUE)
print(myplot)

order=c("enterotype__1","enterotype__2")
myplot<-ggplot(ent_samples, aes(x=factor(Enterotype, levels = order), y=weight)) +
  geom_boxplot() +
  xlab("") +
  ylab("Weight (residuals)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("enterotype__1","enterotype__2")),
              map_signif_level = TRUE)
print(myplot)

############################################################################# CHARACTERIZING SIG VS NON-SIG

# Abundance vs Prevalence

## Scatter plot
plot(filtered_VCs$Prevalence,filtered_VCs$Abundance, xlab="Prevalence (%)", ylab="Mean relative abundance (%)", cex.lab = 1.3)
cor_result=cor.test(filtered_VCs$Prevalence,filtered_VCs$Abundance,method="pearson")
r_value <- cor_result$estimate
p_value <- cor_result$p.value
if (p_value < 0.01){
  P="P < 0.01"
} else {
  if (p_value > 0.01 && p_value < 0.05){
    P="P < 0.05"
  } else {
    P="P > 0.05"
  }
}
legend("topleft", legend = c(paste("rho =", round(r_value, 3)), P), col = c("transparent", "transparent"), pch = 1, bg="white")

# Relative abundances
## Scatter plot
plot(filtered_VCs$Abundance,filtered_VCs$Heritability, xlab="Mean relative abundance (%)", ylab="Species heritability (%)", cex.lab = 1.3)
cor_result=cor.test(filtered_VCs$Abundance,filtered_VCs$Heritability,method="pearson")
r_value <- cor_result$estimate
p_value <- cor_result$p.value
if (p_value < 0.01){
  P="P < 0.01"
} else {
  if (p_value > 0.01 && p_value < 0.05){
    P="P < 0.05"
  } else {
    P="P > 0.05"
  }
}
legend("topright", legend = c(paste("rho =", round(r_value, 3)), P), col = c("transparent", "transparent"), pch = 1, bg="white")

## Boxplot
order=c("Non-heritable","Heritable")
myplot<-ggplot(filtered_VCs, aes(x=factor(All_Sig, levels = order), y=Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  xlab("") +
  ylab("Mean relative abundance (%)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  ylim(0, 0.15) +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Non-heritable","Heritable")),
              map_signif_level = TRUE)
print(myplot)

# Genome Size
## Scatter plot
plot(filtered_VCs$Genome_size,filtered_VCs$Heritability, xlab="Genome Size (AAs)", ylab="Species heritability (%)", cex.lab = 1.3)
cor_result=cor.test(filtered_VCs$Genome_size,filtered_VCs$Heritability,method="pearson")
r_value <- cor_result$estimate
p_value <- cor_result$p.value
if (p_value < 0.01){
  P="P < 0.01"
} else {
  if (p_value > 0.01 && p_value < 0.05){
    P="P < 0.05"
  } else {
    P="P > 0.05"
  }
}
legend("topleft", legend = c(paste("rho =", round(r_value, 3)), P), col = c("transparent", "transparent"), pch = 1, bg="white")
## Boxplot
order=c("Non-heritable","Heritable")
myplot<-ggplot(filtered_VCs, aes(x=factor(All_Sig, levels = order), y=Genome_size)) +
  geom_boxplot() +
  xlab("") +
  ylab("Genome Size (AAs)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Non-heritable","Heritable")),
              map_signif_level = TRUE)
print(myplot)

# Number of Genes
## Scatter plot
plot(filtered_VCs$Number_of_genes,filtered_VCs$Heritability, xlab="Number of Genes", ylab="Species heritability (%)", cex.lab = 1.3)
cor_result=cor.test(filtered_VCs$Number_of_genes,filtered_VCs$Heritability,method="pearson")
r_value <- cor_result$estimate
p_value <- cor_result$p.value
if (p_value < 0.01){
  P="P < 0.01"
} else {
  if (p_value > 0.01 && p_value < 0.05){
    P="P < 0.05"
  } else {
    P="P > 0.05"
  }
}
legend("topleft", legend = c(paste("rho =", round(r_value, 3)), P), col = c("transparent", "transparent"), pch = 1, bg="white")
## Boxplot
order=c("Non-heritable","Heritable")
myplot<-ggplot(filtered_VCs, aes(x=factor(All_Sig, levels = order), y=Number_of_genes)) +
  geom_boxplot() +
  xlab("") +
  ylab("Number of Genes") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Non-heritable","Heritable")),
              map_signif_level = TRUE)
print(myplot)

# Mean Gene Size
## Scatter plot
plot(filtered_VCs$Mean_gene_size,filtered_VCs$Heritability, xlab="Mean Gene Size (AAs)", ylab="Species heritability (%)", cex.lab = 1.3)
cor_result=cor.test(filtered_VCs$Mean_gene_size,filtered_VCs$Heritability,method="pearson")
r_value <- cor_result$estimate
p_value <- cor_result$p.value
if (p_value < 0.01){
  P="P < 0.01"
} else {
  if (p_value > 0.01 && p_value < 0.05){
    P="P < 0.05"
  } else {
    P="P > 0.05"
  }
}
legend("topleft", legend = c(paste("rho =", round(r_value, 3)), P), col = c("transparent", "transparent"), pch = 1, bg="white")
## Boxplot
order=c("Non-heritable","Heritable")
myplot<-ggplot(filtered_VCs, aes(x=factor(All_Sig, levels = order), y=Mean_gene_size)) +
  geom_boxplot() +
  xlab("") +
  ylab("Mean Gene Size (AAs)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Non-heritable","Heritable")),
              map_signif_level = TRUE)
print(myplot)


##################################################################### CHARACTERIZING HUBS

# Genome Size
## Scatter plot
order=c("Non-hub","Cluster Hub","Global Hub","Global and Cluster Hub")
myplot<-ggplot(filtered_VCs, aes(x=factor(Hubs, levels = order), y=Genome_size)) +
  geom_boxplot() +
  xlab("") +
  ylab("Genome Size (AAs)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Non-hub","Cluster Hub"),
                                 c("Cluster Hub","Global Hub"),
                                 c("Global Hub","Global and Cluster Hub")
  ),
  map_signif_level = TRUE)
print(myplot)

# Number of Genes
order=c("Non-hub","Cluster Hub","Global Hub","Global and Cluster Hub")
myplot<-ggplot(filtered_VCs, aes(x=factor(Hubs, levels = order), y=Number_of_genes)) +
  geom_boxplot() +
  xlab("") +
  ylab("Number of Genes") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Non-hub","Cluster Hub"),
                                 c("Cluster Hub","Global Hub"),
                                 c("Global Hub","Global and Cluster Hub")
  ),
  map_signif_level = TRUE)
print(myplot)

# Mean Gene Size
order=c("Non-hub","Cluster Hub","Global Hub","Global and Cluster Hub")
myplot<-ggplot(filtered_VCs, aes(x=factor(Hubs, levels = order), y=Mean_gene_size)) +
  geom_boxplot() +
  xlab("") +
  ylab("Mean Gene Size (AAs)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Non-hub","Cluster Hub"),
                                 c("Cluster Hub","Global Hub"),
                                 c("Global Hub","Global and Cluster Hub")
  ),
  map_signif_level = TRUE)
print(myplot)

##################################################################### CHARACTERIZING TAXA

# Mean Gene Size by Phylum
myplot<-ggplot(filtered_VCs, aes(x=reorder(Phylum,Mean_gene_size,FUN = median), y=Mean_gene_size)) +
  geom_boxplot() +
  xlab("") +
  ylab("Mean Gene Size (AAs)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60))
#geom_signif(comparisons = list(c("p__Bacteroidota","p__Actinobacteriota"),c("p__Actinobacteriota","p__Spirochaetota")),
#            map_signif_level = TRUE)
print(myplot)

# Mean Gene Size by Family
myplot<-ggplot(filtered_VCs[filtered_VCs$Phylum=="p__Bacteroidota",], aes(x=reorder(Family,Mean_gene_size,FUN = median), y=Mean_gene_size)) +
  geom_boxplot() +
  xlab("") +
  ylab("Mean Gene Size (AAs)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60))
#geom_signif(comparisons = list(c("f__Bacteroidaceae","f__Muribaculaceae")),
#            map_signif_level = TRUE)
print(myplot)

# Mean Gene Size by Genus
myplot<-ggplot(filtered_VCs[filtered_VCs$Family=="f__Bacteroidaceae",], aes(x=reorder(Genus,Mean_gene_size,FUN = median), y=Mean_gene_size)) +
  geom_boxplot() +
  xlab("") +
  ylab("Mean Gene Size (AAs)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60))
#geom_signif(comparisons = list(c("f__Bacteroidaceae","f__Muribaculaceae")),
#            map_signif_level = TRUE)
print(myplot)

# Mean Gene Size by Cluster
myplot<-ggplot(filtered_VCs, aes(x=reorder(Cluster,Mean_gene_size,FUN = median), y=Mean_gene_size)) +
  geom_boxplot() +
  xlab("") +
  ylab("Mean Gene Size (AAs)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60))
#geom_signif(comparisons = list(c("cluster__3","cluster__11"),c("cluster__11","cluster__10"),c("cluster__10","cluster__4")),
#            map_signif_level = TRUE)
print(myplot)


######### CENTRALTY

# Degree centrality
## Scatter plot
plot(filtered_VCs$Degree,filtered_VCs$Heritability, xlab="Species degree centrality", ylab="Species heritability (%)", cex.lab = 1.3)
cor_result=cor.test(filtered_VCs$Degree,filtered_VCs$Heritability,method="spearman")
r_value <- cor_result$estimate
p_value <- cor_result$p.value
if (p_value < 0.01){
  P="P < 0.01"
} else {
  if (p_value > 0.01 && p_value < 0.05){
    P="P < 0.05"
  } else {
    P="P > 0.05"
  }
}
legend("topleft", legend = c(paste("rho =", round(r_value, 3)), P), col = c("transparent", "transparent"), pch = 1, bg="white")
## Boxplot
order=c("Non-heritable","Heritable")
myplot<-ggplot(filtered_VCs, aes(x=factor(All_Sig, levels = order), y=Degree)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species degree centrality") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Non-heritable","Heritable")),
              map_signif_level = TRUE)
print(myplot)

# Closeness centrality
## Scatter plot
plot(filtered_VCs$Closeness,filtered_VCs$Heritability, xlab="Species closeness centrality", ylab="Species heritability (%)", cex.lab = 1.3)
cor_result=cor.test(filtered_VCs$Closeness,filtered_VCs$Heritability,method="spearman")
r_value <- cor_result$estimate
p_value <- cor_result$p.value
if (p_value < 0.01){
  P="P < 0.01"
} else {
  if (p_value > 0.01 && p_value < 0.05){
    P="P < 0.05"
  } else {
    P="P > 0.05"
  }
}
legend("topleft", legend = c(paste("rho =", round(r_value, 3)), P), col = c("transparent", "transparent"), pch = 1, bg="white")
## Boxplot
order=c("Non-heritable","Heritable")
myplot<-ggplot(filtered_VCs, aes(x=factor(All_Sig, levels = order), y=Closeness)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species closeness centrality") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Non-heritable","Heritable")),
              map_signif_level = TRUE)
print(myplot)

# Betweeness centrality
## Scatter plot
plot(filtered_VCs$Betweenness,filtered_VCs$Heritability, xlab="Species betweenness centrality", ylab="Species heritability (%)", cex.lab = 1.3)
cor_result=cor.test(filtered_VCs$Betweenness,filtered_VCs$Heritability,method="spearman")
r_value <- cor_result$estimate
p_value <- cor_result$p.value
if (p_value < 0.01){
  P="P < 0.01"
} else {
  if (p_value > 0.01 && p_value < 0.05){
    P="P < 0.05"
  } else {
    P="P > 0.05"
  }
}
legend("topleft", legend = c(paste("rho =", round(r_value, 3)), P), col = c("transparent", "transparent"), pch = 1, bg="white")
## Boxplot
order=c("Non-heritable","Heritable")
myplot<-ggplot(filtered_VCs, aes(x=factor(All_Sig, levels = order), y=Betweenness)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species betweenness centrality") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Non-heritable","Heritable")),
              map_signif_level = TRUE)
print(myplot)

# EIGENVECTOR
## Scatter plot
plot(filtered_VCs$Eigenvector,filtered_VCs$Heritability, xlab="Species eigenvector", ylab="Species heritability (%)", cex.lab = 1.3)
cor_result=cor.test(filtered_VCs$Eigenvector,filtered_VCs$Heritability,method="spearman")
r_value <- cor_result$estimate
p_value <- cor_result$p.value
if (p_value < 0.01){
  P="P < 0.01"
} else {
  if (p_value > 0.01 && p_value < 0.05){
    P="P < 0.05"
  } else {
    P="P > 0.05"
  }
}
legend("topleft", legend = c(paste("rho =", round(r_value, 3)), P), col = c("transparent", "transparent"), pch = 1, bg="white")
## Boxplot
order=c("Non-heritable","Heritable")
myplot<-ggplot(filtered_VCs, aes(x=factor(All_Sig, levels = order), y=Eigenvector)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species eigenvector") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Non-heritable","Heritable")),
              map_signif_level = TRUE)
print(myplot)

###########################



# Correlation with beta diversity
## Scatter plot
plot(filtered_VCs$correlations_beta,filtered_VCs$Heritability, xlab="Correlation with beta diversity (PCoA1)", ylab="Species heritability (%)", cex.lab = 1.3)
cor_result=cor.test(filtered_VCs$correlations_beta,filtered_VCs$Heritability,method="pearson")
r_value <- cor_result$estimate
p_value <- cor_result$p.value
if (p_value < 0.01){
  P="P < 0.01"
} else {
  if (p_value > 0.01 && p_value < 0.05){
    P="P < 0.05"
  } else {
    P="P > 0.05"
  }
}
legend("topleft", legend = c(paste("rho =", round(r_value, 3)), P), col = c("transparent", "transparent"), pch = 1, bg="white")

# Correlation with beta diversity (Bacteroidota)
## Scatter plot
plot(all_traits_phylum[2,],beta, xlab="Bacteroidota relative abundance (CLR residuals)", ylab="Beta diversity (PCoA1)")
cor_result=cor.test(all_traits_phylum[2,],beta,method="pearson")
r_value <- cor_result$estimate
p_value <- cor_result$p.value
if (p_value < 0.01){
  P="P < 0.01"
} else {
  if (p_value > 0.01 && p_value < 0.05){
    P="P < 0.05"
  } else {
    P="P > 0.05"
  }
}
legend("topleft", legend = c(paste("rho =", round(r_value, 3)), P), col = c("transparent", "transparent"), pch = 1, bg="white")

## Boxplot
order=c("Non-heritable","Heritable")
myplot<-ggplot(filtered_VCs, aes(x=factor(All_Sig, levels = order), y=correlations_beta)) +
  geom_boxplot() +
  xlab("") +
  ylab("Correlation with beta diversity (PCoA1)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Non-heritable","Heritable")),
              map_signif_level = TRUE)
print(myplot)

# Correlation with alpha diversity
## Scatter plot
plot(filtered_VCs$correlations_alpha,filtered_VCs$Heritability, xlab="Correlation with alpha diversity (Rao's quadratic index)", ylab="Species heritability (%)", cex.lab = 1.3)
cor_result=cor.test(filtered_VCs$correlations_alpha,filtered_VCs$Heritability,method="pearson")
r_value <- cor_result$estimate
p_value <- cor_result$p.value
if (p_value < 0.01){
  P="P < 0.01"
} else {
  if (p_value > 0.01 && p_value < 0.05){
    P="P < 0.05"
  } else {
    P="P > 0.05"
  }
}
legend("topleft", legend = c(paste("rho =", round(r_value, 3)), P), col = c("transparent", "transparent"), pch = 1, bg="white")
## Boxplot
order=c("Non-heritable","Heritable")
myplot<-ggplot(filtered_VCs, aes(x=factor(All_Sig, levels = order), y=correlations_alpha)) +
  geom_boxplot() +
  xlab("") +
  ylab("Correlation with alpha diversity (Rao's quadratic index)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Non-heritable","Heritable")),
              map_signif_level = TRUE)
print(myplot)

## Correlation Beta and Alpha
plot(beta,alpha, xlab="Beta diversity (PCoA1)", ylab="Alpha diversity (Rao's quadratic index)")
cor_result=cor.test(beta,alpha,method="pearson")
r_value <- cor_result$estimate
p_value <- cor_result$p.value
if (p_value < 0.01){
  P="P < 0.01"
} else {
  if (p_value > 0.01 && p_value < 0.05){
    P="P < 0.05"
  } else {
    P="P > 0.05"
  }
}
legend("topleft", legend = c(paste("rho =", round(r_value, 3)), P), col = c("transparent", "transparent"), pch = 1, bg="white")
plot(filtered_VCs$correlations_beta,filtered_VCs$correlations_alpha, xlab="Correlation with beta diversity (PCoA1)", ylab="Correlation with alpha diversity (Rao's quadratic index)", cex.lab = 1.3)
cor_result=cor.test(filtered_VCs$correlations_beta,filtered_VCs$correlations_alpha,method="pearson")
r_value <- cor_result$estimate
p_value <- cor_result$p.value
if (p_value < 0.01){
  P="P < 0.01"
} else {
  if (p_value > 0.01 && p_value < 0.05){
    P="P < 0.05"
  } else {
    P="P > 0.05"
  }
}
legend("topleft", legend = c(paste("rho =", round(r_value, 3)), P), col = c("transparent", "transparent"), pch = 1, bg="white")

# Correlation with beta by Phylum
myplot<-ggplot(filtered_VCs, aes(x=reorder(Phylum,correlations_beta,FUN = median), y=correlations_beta)) +
  geom_boxplot() +
  xlab("") +
  ylab("Correlation with beta diversity") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("p__Bacteroidota","p__Firmicutes_A")),
              map_signif_level = TRUE)
print(myplot)

# Correlation with beta by Family
myplot<-ggplot(filtered_VCs[filtered_VCs$Phylum=="p__Bacteroidota",], aes(x=reorder(Family,correlations_beta,FUN = median), y=correlations_beta)) +
  geom_boxplot() +
  xlab("") +
  ylab("Correlation with beta diversity") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60))
#geom_signif(comparisons = list(c("f__Bacteroidaceae","f__Muribaculaceae")),
#            map_signif_level = TRUE)
print(myplot)

# Correlation with beta by Genus
myplot<-ggplot(filtered_VCs[filtered_VCs$Family=="f__Bacteroidaceae",], aes(x=reorder(Genus,correlations_beta,FUN = median), y=correlations_beta)) +
  geom_boxplot() +
  xlab("") +
  ylab("Correlation with beta diversity") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60))
#geom_signif(comparisons = list(c("f__Bacteroidaceae","f__Muribaculaceae")),
#            map_signif_level = TRUE)
print(myplot)


# Correlation with beta by Cluster
myplot<-ggplot(filtered_VCs, aes(x=reorder(Cluster,correlations_beta,FUN = median), y=correlations_beta)) +
  geom_boxplot() +
  xlab("") +
  ylab("Correlation with beta diversity") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60))
#geom_signif(comparisons = list(c("cluster__3","cluster__10"),c("cluster__11","cluster__10"),c("cluster__11","cluster__13")),
#            map_signif_level = TRUE)
print(myplot)

# Correlation with alpha by Phylum
myplot<-ggplot(filtered_VCs, aes(x=reorder(Phylum,correlations_alpha,FUN = median), y=correlations_alpha)) +
  geom_boxplot() +
  xlab("") +
  ylab("Correlation with alpha diversity") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("p__Bacteroidota","p__Firmicutes_A")),
              map_signif_level = TRUE)
print(myplot)

# Correlation with alpha by Family
myplot<-ggplot(filtered_VCs[filtered_VCs$Phylum=="p__Bacteroidota",], aes(x=reorder(Family,correlations_alpha,FUN = median), y=correlations_alpha)) +
  geom_boxplot() +
  xlab("") +
  ylab("Correlation with alpha diversity") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60))
#geom_signif(comparisons = list(c("f__Bacteroidaceae","f__Muribaculaceae")),
#            map_signif_level = TRUE)
print(myplot)

# Correlation with alpha by Genus
myplot<-ggplot(filtered_VCs[filtered_VCs$Family=="f__Bacteroidaceae",], aes(x=reorder(Genus,correlations_alpha,FUN = median), y=correlations_alpha)) +
  geom_boxplot() +
  xlab("") +
  ylab("Correlation with alpha diversity") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60))
#geom_signif(comparisons = list(c("f__Bacteroidaceae","f__Muribaculaceae")),
#            map_signif_level = TRUE)
print(myplot)

# Correlation with alpha by Enterotype
myplot<-ggplot(filtered_VCs, aes(x=reorder(Enterotype,correlations_alpha,FUN = median), y=correlations_alpha)) +
  geom_boxplot() +
  xlab("") +
  ylab("Correlation with alpha diversity") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("enterotype__1","enterotype__2")),
              map_signif_level = TRUE)
print(myplot)


# Correlation with alpha by Cluster
myplot<-ggplot(filtered_VCs, aes(x=reorder(Cluster,correlations_alpha,FUN = median), y=correlations_alpha)) +
  geom_boxplot() +
  xlab("") +
  ylab("Correlation with alpha diversity") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60))
#geom_signif(comparisons = list(c("cluster__3","cluster__10"),c("cluster__11","cluster__10"),c("cluster__11","cluster__13")),
#            map_signif_level = TRUE)
print(myplot)



## Heritability by hub
#myplot<-ggplot(filtered_VCs[filtered_VCs$Genus=="g__Prevotella" & filtered_VCs$Cluster=="cluster__10",], aes(x=Hub, y=Heritability)) +
##myplot<-ggplot(filtered_VCs, aes(x=Hub, y=Heritability)) +
#  geom_boxplot() +
#  xlab("") +
#  ylab("Species heritability (%)") +
#  ggtitle(paste0("Subset: ",subset_id)) +
#  theme_bw() +
#  theme(legend.position="bottom",
#        axis.line = element_line(),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        panel.border = element_blank(),
#        panel.background = element_blank(),
#        axis.text = element_text(size = 13),
#        axis.title = element_text(size = 15),
#        legend.text=element_text(size=12),legend.title=element_text(size=15),
#        plot.title = element_text(hjust = 0.5))+
#  scale_x_discrete(guide = guide_axis(angle = 60)) +
#  geom_signif(comparisons = list(c("Hubs","Others")),
#              map_signif_level = TRUE)
#print(myplot)

# Degree Centrality by Phylum
myplot<-ggplot(filtered_VCs, aes(x=reorder(Phylum,Degree,FUN = median), y=Degree)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species degree centrality") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("p__Bacteroidota","p__Firmicutes_A")),
              map_signif_level = TRUE)
print(myplot)

# Degree Centrality by Family
myplot<-ggplot(filtered_VCs[filtered_VCs$Phylum=="p__Bacteroidota",], aes(x=reorder(Family,Degree,FUN = median), y=Degree)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species degree centrality") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60))
#geom_signif(comparisons = list(c("f__Bacteroidaceae","f__Muribaculaceae")),
#            map_signif_level = TRUE)
print(myplot)

# Degree Centrality by Genus
myplot<-ggplot(filtered_VCs[filtered_VCs$Family=="f__Bacteroidaceae",], aes(x=reorder(Genus,Degree,FUN = median), y=Degree)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species degree centrality") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60))
#geom_signif(comparisons = list(c("f__Bacteroidaceae","f__Muribaculaceae")),
#            map_signif_level = TRUE)
print(myplot)

# Degree Centrality by Enterotype
myplot<-ggplot(filtered_VCs, aes(x=reorder(Enterotype,Degree,FUN = median), y=Degree)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species degree centrality") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("enterotype__1","enterotype__2")),
              map_signif_level = TRUE)
print(myplot)

###################################### Cluster

# Degree Centrality by Cluster
myplot<-ggplot(filtered_VCs, aes(x=reorder(Cluster,Degree,FUN = median), y=Degree)) +
  geom_boxplot() +
  xlab("") +
  ylab("Cluster degree centrality") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60))
#geom_signif(comparisons = list(c("cluster__3","cluster__10"),c("cluster__11","cluster__10"),c("cluster__11","cluster__13")),
#            map_signif_level = TRUE)
print(myplot)

# Closeness Centrality by Cluster
myplot<-ggplot(filtered_VCs, aes(x=reorder(Cluster,Closeness,FUN = median), y=Closeness)) +
  geom_boxplot() +
  xlab("") +
  ylab("Cluster closeness centrality") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60))
#geom_signif(comparisons = list(c("cluster__3","cluster__10"),c("cluster__11","cluster__10"),c("cluster__11","cluster__13")),
#            map_signif_level = TRUE)
print(myplot)

# Betweeness Centrality by Cluster
myplot<-ggplot(filtered_VCs, aes(x=reorder(Cluster,Betweeness,FUN = median), y=Betweeness)) +
  geom_boxplot() +
  xlab("") +
  ylab("Cluster betweeness centrality") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60))
#geom_signif(comparisons = list(c("cluster__3","cluster__10"),c("cluster__11","cluster__10"),c("cluster__11","cluster__13")),
#            map_signif_level = TRUE)
print(myplot)

# Heritability by Cluster
myplot<-ggplot(filtered_VCs, aes(x=reorder(Cluster,Heritability,FUN = median), y=Heritability)) +
  geom_boxplot() +
  xlab("") +
  ylab("Cluster Heritability (%)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60))
#geom_signif(comparisons = list(c("cluster__3","cluster__10"),c("cluster__11","cluster__10"),c("cluster__11","cluster__13")),
#            map_signif_level = TRUE)
print(myplot)

# Maternal by Cluster
myplot<-ggplot(filtered_VCs, aes(x=reorder(Cluster,Maternal,FUN = median), y=Maternal)) +
  geom_boxplot() +
  xlab("") +
  ylab("Cluster yMaternal effects (%)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60))
#geom_signif(comparisons = list(c("cluster__3","cluster__10"),c("cluster__11","cluster__10"),c("cluster__11","cluster__13")),
#            map_signif_level = TRUE)
print(myplot)

# Cohousing by Cluster
myplot<-ggplot(filtered_VCs, aes(x=reorder(Cluster,Cohousing,FUN = median), y=Cohousing)) +
  geom_boxplot() +
  xlab("") +
  ylab("Cluster Cohousing effects (%)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60))
#geom_signif(comparisons = list(c("cluster__3","cluster__10"),c("cluster__11","cluster__10"),c("cluster__11","cluster__13")),
#            map_signif_level = TRUE)
print(myplot)




# Closeness Centrality by Phylum
myplot<-ggplot(filtered_VCs, aes(x=reorder(Phylum,Closeness,FUN = median), y=Closeness)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species closeness centrality") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("p__Bacteroidota","p__Firmicutes_A")),
              map_signif_level = TRUE)
print(myplot)

# Closeness Centrality by Family
myplot<-ggplot(filtered_VCs[filtered_VCs$Phylum=="p__Bacteroidota",], aes(x=reorder(Family,Closeness,FUN = median), y=Closeness)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species closeness centrality") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60))
#geom_signif(comparisons = list(c("f__Bacteroidaceae","f__Muribaculaceae")),
#            map_signif_level = TRUE)
print(myplot)



# Closeness Centrality by Genus
myplot<-ggplot(filtered_VCs[filtered_VCs$Family=="f__Bacteroidaceae",], aes(x=reorder(Genus,Closeness,FUN = median), y=Closeness)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species closeness centrality") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60))
#geom_signif(comparisons = list(c("f__Bacteroidaceae","f__Muribaculaceae")),
#            map_signif_level = TRUE)
print(myplot)

# Closeness Centrality by Enterotype
myplot<-ggplot(filtered_VCs, aes(x=reorder(Enterotype,Closeness,FUN = median), y=Closeness)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species closeness centrality") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("enterotype__1","enterotype__2")),
              map_signif_level = TRUE)
print(myplot)

# Closeness Centrality by Cluster
myplot<-ggplot(filtered_VCs, aes(x=reorder(Cluster,Closeness,FUN = median), y=Closeness)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species closeness centrality") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60))
#geom_signif(comparisons = list(c("cluster__3","cluster__10"),c("cluster__11","cluster__10"),c("cluster__11","cluster__13")),
#            map_signif_level = TRUE)
print(myplot)

## Comparisons between Bacteroides and Prevotella
#filtered_VCs_v2=filtered_VCs[filtered_VCs$Genus == "g__Prevotella" | filtered_VCs$Genus == "g__Bacteroides" | filtered_VCs$Genus == "g__CAG-485" | filtered_VCs$Genus == "g__Acetatifactor",]
filtered_VCs_v2=filtered_VCs[filtered_VCs$Genus == "g__Prevotella" | filtered_VCs$Genus == "g__Bacteroides" | filtered_VCs$Genus == "g__CAG-485",]

# Heritability
order=c("g__CAG-485","g__Prevotella","g__Bacteroides")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Genus, levels = order), y=Heritability)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species heritability (%)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("g__Bacteroides","g__Prevotella"),c("g__CAG-485","g__Prevotella")),map_signif_level = TRUE)
#geom_signif(comparisons = list(c("g__Bacteroides","g__Prevotella")),map_signif_level = TRUE)
print(myplot)

# Relative Abundance
order=c("g__CAG-485","g__Prevotella","g__Bacteroides")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Genus, levels = order), y=Abundance)) +
  geom_boxplot() +
  xlab("") +
  ylab("Mean relative abundance (%)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("g__Bacteroides","g__Prevotella"),c("g__CAG-485","g__Prevotella")),map_signif_level = TRUE)
#geom_signif(comparisons = list(c("g__Bacteroides","g__Prevotella")),map_signif_level = TRUE)
print(myplot)

# Prevalence
order=c("g__CAG-485","g__Prevotella","g__Bacteroides")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Genus, levels = order), y=Prevalence)) +
  geom_boxplot() +
  xlab("") +
  ylab("Prevalence (%)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("g__Bacteroides","g__Prevotella"),c("g__CAG-485","g__Prevotella")),map_signif_level = TRUE)
#geom_signif(comparisons = list(c("g__Bacteroides","g__Prevotella")),map_signif_level = TRUE)
print(myplot)

# Genome Size
order=c("g__Prevotella","g__Bacteroides")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Genus, levels = order), y=Genome_size)) +
  geom_boxplot() +
  xlab("") +
  ylab("Genome Size (AAs)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("g__Bacteroides","g__Prevotella"),c("g__CAG-485","g__Prevotella")),map_signif_level = TRUE)
#geom_signif(comparisons = list(c("g__Bacteroides","g__Prevotella")),map_signif_level = TRUE)
print(myplot)

# Number of Genes
order=c("g__CAG-485","g__Prevotella","g__Bacteroides")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Genus, levels = order), y=Number_of_genes)) +
  geom_boxplot() +
  xlab("") +
  ylab("Number of Genes") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("g__Bacteroides","g__Prevotella"),c("g__CAG-485","g__Prevotella")),map_signif_level = TRUE)
#geom_signif(comparisons = list(c("g__Bacteroides","g__Prevotella")),map_signif_level = TRUE)
print(myplot)

# Mean Gene Size
order=c("g__CAG-485","g__Prevotella","g__Bacteroides")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Genus, levels = order), y=Mean_gene_size)) +
  geom_boxplot() +
  xlab("") +
  ylab("Mean Gene Size (AAs)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("g__Bacteroides","g__Prevotella"),c("g__CAG-485","g__Prevotella")),map_signif_level = TRUE)
#geom_signif(comparisons = list(c("g__Bacteroides","g__Prevotella")),map_signif_level = TRUE)
print(myplot)

# Corrrelations Beta
order=c("g__Acetatifactor","g__CAG-485","g__Prevotella","g__Bacteroides")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Genus, levels = order), y=correlations_beta)) +
  geom_boxplot() +
  xlab("") +
  ylab("Correlations with beta diversity (PCoA1)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  #geom_signif(comparisons = list(c("g__Bacteroides","g__Prevotella"),c("g__CAG-485","g__Prevotella"),c("g__CAG-485","g__Acetatifactor")),map_signif_level = TRUE)
  geom_signif(comparisons = list(c("g__Bacteroides","g__Prevotella")),map_signif_level = TRUE)
print(myplot)

# Corrrelations Alpha
order=c("g__Acetatifactor","g__CAG-485","g__Prevotella","g__Bacteroides")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Genus, levels = order), y=correlations_alpha)) +
  geom_boxplot() +
  xlab("") +
  ylab("Correlations with alpha diversity (Rao's quadratic index)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  #geom_signif(comparisons = list(c("g__Bacteroides","g__Prevotella"),c("g__CAG-485","g__Prevotella"),c("g__CAG-485","g__Acetatifactor")),map_signif_level = TRUE)
  geom_signif(comparisons = list(c("g__Bacteroides","g__Prevotella")),map_signif_level = TRUE)
print(myplot)

# Degree
order=c("g__Prevotella","g__Bacteroides")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Genus, levels = order), y=Degree)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species degree centrality") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("g__Bacteroides","g__Prevotella"),c("g__CAG-485","g__Prevotella")),map_signif_level = TRUE)
#geom_signif(comparisons = list(c("g__Bacteroides","g__Prevotella")),map_signif_level = TRUE)
print(myplot)

# Closeness
order=c("g__CAG-485","g__Prevotella","g__Bacteroides")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Genus, levels = order), y=Closeness)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species closeness centrality") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("g__Bacteroides","g__Prevotella"),c("g__CAG-485","g__Prevotella")),map_signif_level = TRUE)
#geom_signif(comparisons = list(c("g__Bacteroides","g__Prevotella")),map_signif_level = TRUE)
print(myplot)

# Betweeness
order=c("g__CAG-485","g__Prevotella","g__Bacteroides")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Genus, levels = order), y=Betweeness)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species betweeness centrality") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("g__Bacteroides","g__Prevotella"),c("g__CAG-485","g__Prevotella")),map_signif_level = TRUE)
#geom_signif(comparisons = list(c("g__Bacteroides","g__Prevotella")),map_signif_level = TRUE)
print(myplot)

## Comparisons between clusters (Prevotella)
genus="g__Prevotella"
filtered_VCs_v2=filtered_VCs[filtered_VCs$Genus == genus,]

# Heritability
order=c("cluster__10","cluster__11")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Cluster, levels = order), y=Heritability)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species heritability (%)") +
  ylim(0,21) +
  ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("cluster__10","cluster__11")),map_signif_level = TRUE)
print(myplot)

# Relative Abundance
order=c("cluster__10","cluster__11")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Cluster, levels = order), y=Abundance)) +
  geom_boxplot() +
  xlab("") +
  ylab("Mean relative abundance (%)") +
  ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("cluster__10","cluster__11")),map_signif_level = TRUE)
print(myplot)

# Prevalence
order=c("cluster__10","cluster__11")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Cluster, levels = order), y=Prevalence)) +
  geom_boxplot() +
  xlab("") +
  ylab("Prevalence (%)") +
  ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("cluster__10","cluster__11")),map_signif_level = TRUE)
print(myplot)

# Genome Size
order=c("cluster__10","cluster__11")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Cluster, levels = order), y=Genome_size)) +
  geom_boxplot() +
  xlab("") +
  ylab("Genome Size (AAs)") +
  ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("cluster__10","cluster__11")),map_signif_level = TRUE)
print(myplot)

# Number of Genes
order=c("cluster__10","cluster__11")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Cluster, levels = order), y=Number_of_genes)) +
  geom_boxplot() +
  xlab("") +
  ylab("Number of Genes") +
  ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("cluster__10","cluster__11")),map_signif_level = TRUE)
print(myplot)

# Mean Gene Size
order=c("cluster__10","cluster__11")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Cluster, levels = order), y=Mean_gene_size)) +
  geom_boxplot() +
  xlab("") +
  ylab("Mean Gene Size (AAs)") +
  ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("cluster__10","cluster__11")),map_signif_level = TRUE)
print(myplot)

# Corrrelations Beta
order=c("cluster__10","cluster__11")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Cluster, levels = order), y=correlations_beta)) +
  geom_boxplot() +
  xlab("") +
  ylab("Correlations with beta diversity (PCoA1)") +
  ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("cluster__10","cluster__11")),map_signif_level = TRUE)
print(myplot)

# Corrrelations Alpha
order=c("cluster__10","cluster__11")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Cluster, levels = order), y=correlations_alpha)) +
  geom_boxplot() +
  xlab("") +
  ylab("Correlations with alpha diversity (Rao's quadratic index)") +
  ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("cluster__10","cluster__11")),map_signif_level = TRUE)
print(myplot)

# Degree
order=c("Heritable","Non-heritable")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(All_Sig, levels = order), y=Degree)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species degree centrality") +
  ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Heritable","Non-heritable")),map_signif_level = TRUE)
print(myplot)

# Closeness
#order=c("cluster__10","cluster__11")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(All_Sig, levels = order), y=Closeness)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species closeness centrality") +
  ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Heritable","Non-heritable")),map_signif_level = TRUE)
print(myplot)

# Betweeness
#order=c("cluster__10","cluster__11")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(All_Sig, levels = order), y=Betweeness)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species betweeness centrality") +
  ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Heritable","Non-heritable")),map_signif_level = TRUE)
print(myplot)


## Comparisons between clusters (Bacteroides)
genus="g__Bacteroides"
filtered_VCs_v2=filtered_VCs[filtered_VCs$Genus == genus,]

# Heritability
order=c("Heritable","Non-heritable")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(All_Sig, levels = order), y=Heritability)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species heritability (%)") +
  ylim(0,21) +
  ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("cluster__10","cluster__3")),map_signif_level = TRUE)
print(myplot)

# Relative Abundance
order=c("cluster__10","cluster__3")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Cluster, levels = order), y=Abundance)) +
  geom_boxplot() +
  xlab("") +
  ylab("Mean relative abundance (%)") +
  ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("cluster__10","cluster__3")),map_signif_level = TRUE)
print(myplot)

# Prevalence
order=c("cluster__10","cluster__3")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Cluster, levels = order), y=Prevalence)) +
  geom_boxplot() +
  xlab("") +
  ylab("Prevalence (%)") +
  ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("cluster__10","cluster__3")),map_signif_level = TRUE)
print(myplot)

# Genome Size
order=c("cluster__10","cluster__3")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Cluster, levels = order), y=Genome_size)) +
  geom_boxplot() +
  xlab("") +
  ylab("Genome Size (AAs)") +
  ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("cluster__10","cluster__3")),map_signif_level = TRUE)
print(myplot)

# Number of Genes
order=c("cluster__10","cluster__3")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Cluster, levels = order), y=Number_of_genes)) +
  geom_boxplot() +
  xlab("") +
  ylab("Number of Genes") +
  ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("cluster__10","cluster__3")),map_signif_level = TRUE)
print(myplot)

# Mean Gene Size
order=c("cluster__10","cluster__3")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Cluster, levels = order), y=Mean_gene_size)) +
  geom_boxplot() +
  xlab("") +
  ylab("Mean Gene Size (AAs)") +
  ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("cluster__10","cluster__3")),map_signif_level = TRUE)
print(myplot)

# Corrrelations Beta
order=c("cluster__10","cluster__3")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Cluster, levels = order), y=correlations_beta)) +
  geom_boxplot() +
  xlab("") +
  ylab("Correlations with beta diversity (PCoA1)") +
  ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("cluster__10","cluster__3")),map_signif_level = TRUE)
print(myplot)

# Corrrelations Alpha
order=c("cluster__10","cluster__3")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Cluster, levels = order), y=correlations_alpha)) +
  geom_boxplot() +
  xlab("") +
  ylab("Correlations with alpha diversity (Rao's quadratic index)") +
  ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("cluster__10","cluster__3")),map_signif_level = TRUE)
print(myplot)










################# "Only Heritable"

# All species


filtered_VCs_v2=filtered_VCs
# Heritable & Degree
order=c("Heritable","Non-heritable")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(All_Sig, levels = order), y=Degree)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species degree centrality") +
  #ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Heritable","Non-heritable")),map_signif_level = TRUE)
print(myplot)

# Heritable & Closeness
order=c("Heritable","Non-heritable")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(All_Sig, levels = order), y=Closeness)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species closeness centrality") +
  #ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Heritable","Non-heritable")),map_signif_level = TRUE)
print(myplot)

# Heritable & Betweeness
order=c("Heritable","Non-heritable")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(All_Sig, levels = order), y=Betweeness)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species betweeness centrality") +
  #ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Heritable","Non-heritable")),map_signif_level = TRUE)
print(myplot)

#order=unique(filtered_VCs_v2$Hubs)
# Hubs & Heritability
order=c("Non-hub",
        "Cluster Hub",
        "Global Hub",
        "Global and Cluster Hub"
)
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Hubs, levels = order), y=Heritability)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species heritability (%)") +
  #ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Non-hub","Cluster Hub"),
                                 c("Cluster Hub","Global Hub"),
                                 c("Global Hub","Global and Cluster Hub")
  ),map_signif_level = TRUE)
print(myplot)

# Hubs & Maternal
order=c("Non-hub",
        "Cluster Hub",
        "Global Hub",
        "Global and Cluster Hub"
)
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Hubs, levels = order), y=Maternal)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species maternal effects (%)") +
  #ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Non-hub","Cluster Hub"),
                                 c("Cluster Hub","Global Hub"),
                                 c("Global Hub","Global and Cluster Hub")
  ),map_signif_level = TRUE)
print(myplot)

# Hubs & Cohousing
order=c("Non-hub",
        "Cluster Hub",
        "Global Hub",
        "Global and Cluster Hub"
)
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Hubs, levels = order), y=Cohousing)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species cohousing effects (%)") +
  #ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Non-hub","Cluster Hub"),
                                 c("Cluster Hub","Global Hub"),
                                 c("Global Hub","Global and Cluster Hub")
  ),map_signif_level = TRUE)
print(myplot)


# Heritablility and Degree
plot(filtered_VCs_v2$Heritability,filtered_VCs_v2$Degree, #main=paste0("Genus: ",gsub("g__","",genus)), 
     xlab="Heritability", ylab="Species degree centrality", cex.lab = 1)
abline(lm(filtered_VCs_v2$Degree ~ filtered_VCs_v2$Heritability),col="red")
cor_result=cor.test(filtered_VCs_v2$Heritability,filtered_VCs_v2$Degree,method="pearson")
r_value <- cor_result$estimate
p_value <- cor_result$p.value
#if (p_value < 0.01){
#  P="P < 0.01"
#} else {
#  if (p_value > 0.01 && p_value < 0.05){
#    P="P < 0.05"
#  } else {
#    P="P > 0.05"
#  }
#}
legend("topleft", legend = c(paste("rho =", round(r_value, 3)), paste0(sprintf("p = %.2e", p_value))), col = c("transparent", "transparent"), pch = 1, bg="white")

# Heritablility and Betweeness
plot(filtered_VCs_v2$Heritability,filtered_VCs_v2$Betweeness, #main=paste0("Genus: ",gsub("g__","",genus)), 
     xlab="Heritability", ylab="Species betweeness centrality", cex.lab = 1)
abline(lm(filtered_VCs_v2$Betweeness ~ filtered_VCs_v2$Heritability),col="red")
cor_result=cor.test(filtered_VCs_v2$Heritability,filtered_VCs_v2$Betweeness,method="pearson")
r_value <- cor_result$estimate
p_value <- cor_result$p.value
legend("topleft", legend = c(paste("rho =", round(r_value, 3)), paste0(sprintf("p = %.2e", p_value))), col = c("transparent", "transparent"), pch = 1, bg="white")



# Maternal and Closeness
plot(filtered_VCs_v2$Maternal,filtered_VCs_v2$Closeness, #main=paste0("Genus: ",gsub("g__","",genus)), 
     xlab="Maternal", ylab="Species closeness centrality", cex.lab = 1)
abline(lm(filtered_VCs_v2$Closeness ~ filtered_VCs_v2$Maternal),col="red")
cor_result=cor.test(filtered_VCs_v2$Maternal,filtered_VCs_v2$Closeness,method="pearson")
r_value <- cor_result$estimate
p_value <- cor_result$p.value
#if (p_value < 0.01){
#  P="P < 0.01"
#} else {
#  if (p_value > 0.01 && p_value < 0.05){
#    P="P < 0.05"
#  } else {
#    P="P > 0.05"
#  }
#}
legend("topleft", legend = c(paste("rho =", round(r_value, 3)), paste0(sprintf("p = %.2e", p_value))), col = c("transparent", "transparent"), pch = 1, bg="white")

genera=c(
  "g__Prevotella",
  "g__Bacteroides"
)
filtered_VCs_v2=filtered_VCs[filtered_VCs$Genus %in% genera & !grepl("g__Bacteroides_F",filtered_VCs$Genus),]

# Genera & Heritability
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Genus, levels = genera), y=Heritability)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species heritability (%)") +
  #ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(genera),map_signif_level = TRUE)
print(myplot)

# Genera & Heritability
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Genus, levels = genera), y=Maternal)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species maternal effects (%)") +
  #ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(genera),map_signif_level = TRUE)
print(myplot)

# Genera & Cohousing
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Genus, levels = genera), y=Cohousing)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species cohousing effects (%)") +
  #ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(genera),map_signif_level = TRUE)
print(myplot)

myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Genus, levels = genera), y=Degree)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species degree centrality") +
  #ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(genera),map_signif_level = TRUE)
print(myplot)

myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Genus, levels = genera), y=Closeness)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species closeness centrality") +
  #ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(genera),map_signif_level = TRUE)
print(myplot)

myplot<-ggplot(filtered_VCs_v2, aes(x=factor(Genus, levels = genera), y=Betweeness)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species betweeness centrality") +
  #ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(genera),map_signif_level = TRUE)
print(myplot)

myplot<-ggplot(filtered_VCs, aes(x=reorder(Cluster,Betweeness,FUN = median), y=Betweeness)) +
  geom_boxplot() +
  xlab("") +
  ylab("Cluster betweeness centrality") +
  #ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60))
#geom_signif(comparisons = list(genera),map_signif_level = TRUE)
print(myplot)



# Heritability vs Maternal vs Cohousing

long_VCs <- filtered_VCs %>%
  pivot_longer(cols = c(Maternal, Heritability, Cohousing),
               names_to = "Component",
               values_to = "Value")

myplot<-ggplot(long_VCs, aes(x=reorder(Component,Value,FUN = median), y=Value)) +
  geom_boxplot() +
  xlab("") +
  ylab("Effect size (%)") +
  #ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Maternal","Heritability"),c("Heritability","Cohousing")),map_signif_level = TRUE)
print(myplot)




################# "Only Heritable"






# Closeness
#order=c("cluster__10","cluster__3")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(All_Sig, levels = order), y=Closeness)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species closeness centrality") +
  ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Heritable","Non-heritable")),map_signif_level = TRUE)
print(myplot)

# Betweeness
#order=c("cluster__10","cluster__3")
myplot<-ggplot(filtered_VCs_v2, aes(x=factor(All_Sig, levels = order), y=Betweeness)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species betweeness centrality") +
  ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Heritable","Non-heritable")),map_signif_level = TRUE)
print(myplot)

dev.off()

myplot<-ggplot(filtered_VCs_v2, aes(x=factor(All_Sig, levels = order), y=Maternal)) +
  geom_boxplot() +
  xlab("") +
  ylab("Species betweeness centrality") +
  ggtitle(paste0("Genus: ",gsub("g__","",genus))) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Heritable","Non-heritable")),map_signif_level = TRUE)
print(myplot)

dev.off()


plot(filtered_VCs_v2$Degree,abs(filtered_VCs_v2$Heritability), xlab="Species degree centrality", ylab="Heritability", cex.lab = 1)
cor_result=cor.test(filtered_VCs_v2$Degree,abs(filtered_VCs_v2$Heritability),method="pearson")
r_value <- cor_result$estimate
p_value <- cor_result$p.value
if (p_value < 0.01){
  P="P < 0.01"
} else {
  if (p_value > 0.01 && p_value < 0.05){
    P="P < 0.05"
  } else {
    P="P > 0.05"
  }
}
legend("topleft", legend = c(paste("rho =", round(r_value, 3)), P), col = c("transparent", "transparent"), pch = 1, bg="white")

plot(filtered_VCs$Closeness,abs(filtered_VCs$correlations_alpha), xlab="Species closeness centrality", ylab="Absolute correlation with alpha diversity (Rao's quadratic index)", cex.lab = 1)
cor_result=cor.test(filtered_VCs$Closeness,abs(filtered_VCs$correlations_alpha),method="pearson")
r_value <- cor_result$estimate
p_value <- cor_result$p.value
if (p_value < 0.01){
  P="P < 0.01"
} else {
  if (p_value > 0.01 && p_value < 0.05){
    P="P < 0.05"
  } else {
    P="P > 0.05"
  }
}
legend("topleft", legend = c(paste("rho =", round(r_value, 3)), P), col = c("transparent", "transparent"), pch = 1, bg="white")




#################################################################################### OLD

# Porportion of phyla in the clusters
prop_data <- filtered_VCs %>%
  group_by(Cluster, Phylum) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()
prop_data$proportion=round(100*prop_data$proportion,digit=2)
myplot <- ggplot(prop_data, aes(x = Cluster, y = proportion, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack", linewidth = 2, col = "white",  width = 1) +
  #scale_fill_manual(values = colours_phylum) +
  labs(fill = "Phylum") +
  ggtitle("Proportion of phyla") +
  xlab("") +
  ylab("Proportion (%)") +
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),)
print(myplot)

# Porportion of bacteroidota families in the clusters
prop_data <- filtered_VCs[filtered_VCs$Phylum=="p__Bacteroidota",] %>%
  group_by(Cluster, Family) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()
prop_data$proportion=round(100*prop_data$proportion,digit=2)
myplot <- ggplot(prop_data, aes(x = Cluster, y = proportion, fill = Family)) +
  geom_bar(stat = "identity", position = "stack", linewidth = 2, col = "white",  width = 1) +
  #scale_fill_manual(values = colours_phylum) +
  labs(fill = "Family") +
  ggtitle("Proportion of Bacteroidota families") +
  xlab("") +
  ylab("Proportion (%)") +
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),)
print(myplot)

# Porportion of genera in clusters 3, 10 and 11
prop_data <- filtered_VCs[filtered_VCs$Cluster=="cluster__3" | filtered_VCs$Cluster=="cluster__10" | filtered_VCs$Cluster=="cluster__11",] %>%
  group_by(Cluster, Genus) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()
prop_data$proportion=round(100*prop_data$proportion,digit=2)
levels_rank <- levels(factor(prop_data$Genus))
set.seed(120)
num_colors <- length(levels_rank)
all_colors <- colors()
random_colors <- sample(all_colors, num_colors, replace = FALSE)
names(random_colors) <- levels_rank
myplot <- ggplot(prop_data, aes(x = Cluster, y = proportion, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack", linewidth = 2, col = "white",  width = 1) +
  #scale_fill_manual(values = colours_phylum) +
  labs(fill = "Genus") +
  ggtitle("") +
  xlab("") +
  ylab("Proportion (%)") +
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  scale_fill_manual(values = random_colors) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=9),)
print(myplot)

# Porportion of significant by phylum
prop_data <- filtered_VCs %>%
  group_by(Phylum, All_Sig) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()
prop_data$proportion=round(100*prop_data$proportion,digit=2)
order=prop_data[grepl("Non-heritable",prop_data$All_Sig),]
order = order[order(order$proportion, decreasing = F),]
order=order$Phylum
myplot <- ggplot(prop_data, aes(x = factor(Phylum, levels = order), y = proportion, fill = All_Sig)) +
  geom_bar(stat = "identity", position = "stack", linewidth = 2, col = "white",  width = 1) +
  #scale_fill_manual(values = colours_phylum) +
  labs(fill = "Species") +
  ggtitle("Proportion of species with significant heritability (FDR<10%)") +
  xlab("") +
  ylab("Proportion (%)") +
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),)
print(myplot)

# Prevalence
## Scatter plot
plot(filtered_VCs$Prevalence,filtered_VCs$Heritability, xlab="Prevalence (%)", ylab="Species heritability (%)", cex.lab = 1.3)
cor_result=cor.test(filtered_VCs$Prevalence,filtered_VCs$Heritability,method="pearson")
r_value <- cor_result$estimate
p_value <- cor_result$p.value
if (p_value < 0.01){
  P="P < 0.01"
} else {
  if (p_value > 0.01 && p_value < 0.05){
    P="P < 0.05"
  } else {
    P="P > 0.05"
  }
}
legend("topleft", legend = c(paste("rho =", round(r_value, 3)), P), col = c("transparent", "transparent"), pch = 1, bg="white")

## Boxplot
order=c("Non-heritable","Heritable")
myplot<-ggplot(filtered_VCs, aes(x=factor(All_Sig, levels = order), y=Prevalence)) +
  geom_boxplot() +
  xlab("") +
  ylab("Prevalence (%)") +
  ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),legend.title=element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  geom_signif(comparisons = list(c("Non-heritable","Heritable")),
              map_signif_level = TRUE)
print(myplot)
