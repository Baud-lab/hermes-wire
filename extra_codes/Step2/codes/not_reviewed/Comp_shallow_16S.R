library(ggsignif)
library(dplyr)
suppressMessages(library(qvalue))

#Load taxonomic table
taxonomy=read.table("/users/abaud/data/secondary/indexes/gtdb/genes/gtdb_ids_taxon.txt", sep="\t")
colnames(taxonomy)=c("genome","taxonomy")
taxonomy$taxonomy<-gsub(" ","_",as.character(taxonomy$taxonomy))
taxonomy$phylum<-sapply(strsplit(taxonomy$taxonomy, ";"), "[",2)
taxonomy$class<-sapply(strsplit(taxonomy$taxonomy, ";"), "[",3)
taxonomy$order<-sapply(strsplit(taxonomy$taxonomy, ";"), "[",4)
taxonomy$family<-sapply(strsplit(taxonomy$taxonomy, ";"), "[",5)
taxonomy$genus<-sapply(strsplit(taxonomy$taxonomy, ";"), "[",6)
taxonomy$species<-sapply(strsplit(taxonomy$taxonomy, ";"), "[",7)

#Load shallow
#Only samples
#load("/users/abaud/fmorillo/P50/microbiome/output/microbiome_profiles/taxonomic_profile/gtdb_profiles/207/step2/abundance_mean/not_norm_comp_16S_207_E1/harmonized_shallow_only_samples/Heritability/heritability.RData")
#Taxa and samples
#load("/nfs/users/abaud/fmorillo/P50/microbiome/output/microbiome_profiles/taxonomic_profile/gtdb_profiles/207/step2/abundance_mean/not_norm_comp_16S_207_E1/harmonized/Heritability/heritability.RData")
#load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Shallow_Full_Harm/Heritability/heritability.RData")
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Shallow_Sample_Harm/Heritability/heritability.RData")

subset_id="ALL" ### Remove
sig_f=0.1
sig_b=0.05
all_VCs_full=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),] ### Remove

all_VCs_full$bonf_pvalue_DGE = all_VCs_full$pvalue_DGE*dim(all_VCs_full)[1]
all_VCs_full$qvalue_DGE = qvalue(all_VCs_full$pvalue_DGE)$qvalue
all_VCs_full$qp_diff = all_VCs_full$qvalue_DGE - all_VCs_full$pvalue_DGE
all_VCs_full <- all_VCs_full %>% 
  mutate( Significance = case_when(
    qvalue_DGE < sig_f & bonf_pvalue_DGE < sig_b ~ "Heritable traits (Bonferroni < 5%)",
    qvalue_DGE < sig_f & bonf_pvalue_DGE >= sig_b ~ "Heritable traits (FDR < 10%)",
    TRUE ~ "Non-heritable traits"))
all_VCs_full = all_VCs_full[order(all_VCs_full$logP, decreasing = T),]
all_VCs_full$qvalue_DGE = round(all_VCs_full$qvalue_DGE, digits = 2)

load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2//Output_Shallow_Full_Harm/Heritability/heritability.RData")

all_shallow=all_VCs_full
all_shallow=all_shallow[grepl("ALL",all_shallow$subset_id),]
all_shallow$method="Shallow Shotgun"
#all_shallow_sig=all_shallow[all_shallow$Significance=="Significant (FDR)",]

#all_shallow_genus=all_shallow[all_shallow$taxa=="Genus",]
#all_shallow_genus_sig=all_shallow_sig[all_shallow_sig$taxa=="Genus",]
#
#all_shallow_NY=all_shallow[grepl("NY",all_shallow$study),]
#all_shallow_MI=all_shallow[grepl("MI",all_shallow$study),]
#
#all_shallow_sig_NY=all_shallow_sig[grepl("NY",all_shallow_sig$study),]
#all_shallow_sig_MI=all_shallow_sig[grepl("MI",all_shallow_sig$study),]

#Load 16S
#Only samples
#load("/users/abaud/fmorillo/P50/microbiome/output/microbiome_profiles/taxonomic_profile/16S/207/step2/mean/harmonized_shallow_only_samples/Heritability/heritability.RData")
#Taxa and samples
#load("~/P50/microbiome/output/microbiome_profiles/taxonomic_profile/16S/207/step2/mean/harmonized_shallow_only_samples/Heritability/heritability.RData")
#load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_16S_Full_Harm/Heritability/heritability.RData")
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_16S_Sample_Harm/Heritability/heritability.RData")

subset_id="ALL" ### Remove
sig_f=0.1
sig_b=0.05
all_VCs_full=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),] ### Remove

all_VCs_full$bonf_pvalue_DGE = all_VCs_full$pvalue_DGE*dim(all_VCs_full)[1]
all_VCs_full$qvalue_DGE = qvalue(all_VCs_full$pvalue_DGE)$qvalue
all_VCs_full$qp_diff = all_VCs_full$qvalue_DGE - all_VCs_full$pvalue_DGE
all_VCs_full <- all_VCs_full %>% 
  mutate( Significance = case_when(
    qvalue_DGE < sig_f & bonf_pvalue_DGE < sig_b ~ "Heritable traits (Bonferroni < 5%)",
    qvalue_DGE < sig_f & bonf_pvalue_DGE >= sig_b ~ "Heritable traits (FDR < 10%)",
    TRUE ~ "Non-heritable traits"))
all_VCs_full = all_VCs_full[order(all_VCs_full$logP, decreasing = T),]
all_VCs_full$qvalue_DGE = round(all_VCs_full$qvalue_DGE, digits = 2)

load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_16S_Full_Harm/Heritability/heritability.RData")

all_16S=all_VCs_full
all_16S=all_16S[grepl("ALL",all_16S$subset_id),]
all_16S$method="16S"


count(all_16S,Significance)
count(all_shallow,Significance)

#all_16S_sig=all_16S[all_16S$Significance=="Significant (FDR)",]

#all_16S_genus=all_shallow[all_shallow$taxa=="Genus",]
#all_16S_genus_sig=all_shallow_sig[all_shallow_sig$taxa=="Genus",]

all_16S=all_16S[intersect(rownames(all_16S),rownames(all_shallow)),]
all_shallow=all_shallow[intersect(rownames(all_16S),rownames(all_shallow)),]

all_methods=rbind(all_shallow,all_16S)
#all_methods_sig=all_methods[all_methods$Significance=="Significant (FDR)",]

#all_shallow_common=all_shallow[all_shallow$trait %in% all_16S$trait,]
#all_16S_common=all_16S[all_16S$trait %in% all_shallow_common$trait,]
#
#all_methods_common=rbind(all_shallow_common,all_16S_common)
#all_methods_common_sig=all_methods[all_methods_common$Significance=="Significant (FDR)",]
#
#### Heritability plot
#myplot<-ggplot(all_methods, aes(x=reorder(method,var_Ad*100,FUN = median), y=var_Ad*100)) + 
#  geom_boxplot() +
#  xlab("Method") +
#  ylab("Heritability (%)") +
#  theme_bw() +
#  theme(legend.position="bottom",
#        axis.line = element_line(),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        panel.border = element_blank(),
#        panel.background = element_blank(),
#        axis.text = element_text(size = 15),
#        axis.title = element_text(size = 15),
#        legend.text=element_text(size=10))+
#  ylim(0,30) +
#  scale_x_discrete(guide = guide_axis(angle = 45)) +
#  geom_signif(comparisons = list(c("Shallow Shotgun","16S")), map_signif_level=TRUE)
#print(myplot)
#
##Volcano plot (just for all)
#max_logP = max(-log10(all_methods[,'pvalue']), na.rm = T)
#max_var_Ad = max(100*(all_methods[,'var_Ad']), na.rm = T)
#pdf('/users/abaud/fmorillo/P50/microbiome/plots/volcano_comp_shallow_16S.pdf',bg='white')
#par(mfrow = c(1,2), mar = c(2,2,2,1))
#all_VCs_plot=all_shallow
#myplot<-ggplot(all_VCs_plot, aes(x=var_Ad*100, y=-log10(pvalue))) + 
#  geom_point(aes(color = Significance), size = 4/5) +
#  ggtitle("Shallow shotgun (GTDB 207)") +
#  xlab(paste("Heritability (%) (",subset_ids,")",sep=' ')) +
#  ylab(expression("-log"[10]*"P-Value Heritability")) +
#  theme_bw() +
#  theme(legend.position="bottom",
#        axis.line = element_line(),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        panel.border = element_blank(),
#        panel.background = element_blank(),
#        axis.text = element_text(size = 22),
#        axis.title = element_text(size = 22)) +
#  coord_cartesian(xlim = c(0, max_var_Ad), ylim = c(0, max_logP)) +
#  scale_y_continuous(breaks = seq(0, max_logP, by = 2)) +
#  scale_x_continuous(breaks = seq(0, max_var_Ad, by = 5)) +
#  scale_color_manual(values = c("black", "deepskyblue4", "firebrick3"),
#                     labels = c("Not-Significant","Significant (Bonferroni)","Significant (FDR)")) +
#  labs(color = "Legend")
#print(myplot)
#dev.off()
#
##Scatterplot
#
##all_16S=all_16S[!grepl("g__Anaerotruncus",all_16S$trait),] #Remove the outlier
#m=intersect(rownames(all_16S),rownames(all_shallow))
#filt_shallow=all_shallow[m,]
#filt_16S=all_16S[m,]
#
#
#plot(x = 100*filt_16S$var_Ad, y = 100*filt_shallow$var_Ad, frame = FALSE, xlab = "Heritability (%) - 16S", ylab = "Heritability (%) - Shallow")
#abline(lm(filt_shallow$var_Ad ~ filt_16S$var_Ad)) # abline(y ~ x)
##dev.off()
#cor.test(filt_16S$var_Ad,filt_shallow$var_Ad, method = "spearman")

### Find common taxa
#x=which(duplicated(all_shallow_sig$trait))
#y=all_shallow_sig$trait[x]
#common=subset(all_shallow_sig, trait %in% y)
#all_filtered=rbind(filt_shallow,filt_16S)
all_filtered=all_methods
all_filtered <- all_filtered %>%
  group_by(trait) %>%
  mutate(paired = group_indices())
all_filtered <- all_filtered %>%
  group_by(paired) %>%
  mutate(Difference = ifelse(
    var_Ad[method == "Shallow Shotgun"] - var_Ad[method == "16S"] > 0,
    "Shallow > 16S",
    "Shallow < 16S"
  ))

#no_species=all_filtered[!grepl("Species",all_filtered$rank),]

ranks=c("All Ranks","Phylum","Class","Order","Family","Genus","Species")
order=c("16S","Shallow Shotgun")
pdf('~/paper_figures/herit_16s_shallow.pdf',bg='white')
for (rank in ranks){
  if (rank=="All Ranks"){
    df=all_filtered
  } else {
    df=all_filtered[all_filtered$rank==rank,]
  }
  p_value <- wilcox.test(df$var_Ad[df$method == "16S"], 
                    df$var_Ad[df$method == "Shallow Shotgun"])$p.value
  maxvar_Ad=max(df$var_Ad)
  myplot<-ggplot(df, aes(x=factor(method, levels = order), y=var_Ad*100)) + 
    geom_point(aes(group=paired),size=2,shape=21, color="lightgrey") +
    xlab("") +
    ylab("Heritability (%)") +
    ggtitle(paste0("Rank: ",rank)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position="bottom",
          axis.line = element_line(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          legend.text=element_text(size=13),legend.title=element_text(size=15))+
    ylim(0,(maxvar_Ad+10)) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    #geom_line(aes(group=paired),color="grey",size = 0.5) +
    geom_line(aes(group = paired, color = Difference), size = 0.5) +  # Use "Difference" for color
    scale_color_manual(values = c("Shallow > 16S" = "lightblue", "Shallow < 16S" = "lightcoral")) +
    geom_boxplot(outlier.shape = NA, fill = "lightgrey", alpha = 0.5) +
    geom_signif(comparisons = list(c("16S", "Shallow Shotgun")), 
                annotations = sprintf("p = %.2e", p_value),
                test = "wilcox.test",
                y_position = maxvar_Ad+9, map_signif_level = FALSE)
  print(myplot)
}
dev.off()
#boxplot(common$var_Ad*100)

com_NY=common[grepl("NY",common$study),]
rownames(com_NY)=com_NY$trait
com_MI=common[grepl("MI",common$study),]
rownames(com_MI)=com_MI$trait
com_MI=com_MI[rownames(com_NY),]
plot(com_NY$var_Ad*100,com_MI$var_Ad*100,xlab="Heritability (%) - NY", ylab="Heritability (%) - MI")
abline()
cor.test

common=common[grepl("species",common$taxa),]
common$genus=taxonomy$genus[match(common$trait, taxonomy$species)]
common$family=taxonomy$family[match(common$trait, taxonomy$species)]
common$phylum=taxonomy$phylum[match(common$trait, taxonomy$species)]

stats_genus=count(common,genus)
stats_family=count(common,family)
stats_phylum=count(common,phylum)
stats_species=count(common,trait)

### Relative abundances
prefix <- "/users/abaud/fmorillo/P50/microbiome/output/microbiome_profiles/taxonomic_profile/gtdb_profiles_10k/"
load(paste(prefix,'taxa_filtered_genomes_matrix_gtdb_noprev.RData',sep=''))
#load(paste(prefix,'full_genomes_matrix_gtdb.RData',sep=''))

#Collapse
for (taxa in c("phylum","class","order","family","genus","species")){
  for (study in c('all','MI','NY')) {
    matrix = get(paste('filtered_genome_',study,sep=''))
    matrix = t(matrix)
    colnames(matrix)=taxonomy[[taxa]][match(colnames(matrix), taxonomy$genome)]
    collapsed_matrix=t(matrix)
    collapsed_matrix=t(sapply(by(collapsed_matrix,rownames(collapsed_matrix),colSums),identity))
    assign(paste('filtered',taxa,study,sep='_'), collapsed_matrix)
    print(paste0(taxa,"/", study,": ",dim(get(paste('filtered',taxa,study,sep='_')))))
  } 
}

### Step 4: Taxa filtering based on prevalence (Important for Differential Abundance Analyses) ###

for (taxa in c("phylum","class","order","family","genus", "species")){
  for (study in c('all','MI','NY')) {
    matrix = get(paste('filtered',taxa,study,sep='_'))
    # Calculate the prevalence of each trait
    spe.nonZero  <- rowSums(matrix>0)/ncol(matrix) # Prevalence
    # Get the list of traits that are present in more than 50% samples
    spe.keep <- rownames(matrix)[spe.nonZero > 0.5]
    length(spe.keep) # Check how many traits are in the list.
    #Remove low abundance species
    matrix.filtered <- matrix[spe.keep,]
    assign(paste('filtered_prev',taxa,study,sep='_'), matrix.filtered)
    print(paste0(taxa,"/", study,": ",dim(get(paste('filtered_prev',taxa,study,sep='_')))))
  } 
}

ab_vs_prev_phylum=data.frame(trait=rownames(filtered_prev_phylum_all),
                             abundance=100*rowMeans(sweep(filtered_prev_phylum_all,2,colSums(filtered_prev_phylum_all),"/")),
                             prevalence=100*rowSums(filtered_prev_phylum_all>0)/ncol(filtered_prev_phylum_all))
ab_vs_prev_phylum <- ab_vs_prev_phylum %>% 
  mutate( Significance = case_when(trait %in% stats_phylum$phylum ~ "Significant traits", TRUE ~ "Non-significant traits"))
ab_vs_prev_family=data.frame(trait=rownames(filtered_prev_family_all),
                             abundance=100*rowMeans(sweep(filtered_prev_family_all,2,colSums(filtered_prev_family_all),"/")),
                             prevalence=100*rowSums(filtered_prev_family_all>0)/ncol(filtered_prev_family_all))
ab_vs_prev_family <- ab_vs_prev_family %>% 
  mutate( Significance = case_when(trait %in% stats_family$family ~ "Significant traits", TRUE ~ "Non-significant traits"))
ab_vs_prev_genus=data.frame(trait=rownames(filtered_prev_genus_all),
                             abundance=100*rowMeans(sweep(filtered_prev_genus_all,2,colSums(filtered_prev_genus_all),"/")),
                             prevalence=100*rowSums(filtered_prev_genus_all>0)/ncol(filtered_prev_genus_all))
ab_vs_prev_genus <- ab_vs_prev_genus %>% 
  mutate( Significance = case_when(trait %in% stats_genus$genus ~ "Significant traits", TRUE ~ "Non-significant traits"))

ab_vs_prev_species=data.frame(trait=rownames(filtered_prev_species_all),
                            abundance=100*rowMeans(sweep(filtered_prev_species_all,2,colSums(filtered_prev_species_all),"/")),
                            prevalence=100*rowSums(filtered_prev_species_all>0)/ncol(filtered_prev_species_all))

motch=match(all_shallow$trait,ab_vs_prev_species$trait)
ab_vs_prev_species=ab_vs_prev_species[motch,]

ab_vs_prev_species <- ab_vs_prev_species %>% 
  mutate( Significance = case_when(trait %in% all_shallow_sig$trait ~ "Significant traits", TRUE ~ "Non-significant traits"))

ab_vs_prev_species$genus<-sapply(strsplit(ab_vs_prev_species$trait, "_"), "[",3)



##NY
ab_vs_prev_species_NY=data.frame(trait=rownames(filtered_prev_species_NY),
                              abundance=100*rowMeans(sweep(filtered_prev_species_NY,2,colSums(filtered_prev_species_NY),"/")),
                              prevalence=100*rowSums(filtered_prev_species_NY>0)/ncol(filtered_prev_species_NY))
motch=match(all_shallow_NY$trait,ab_vs_prev_species_NY$trait)
ab_vs_prev_species_NY=ab_vs_prev_species_NY[motch,]
ab_vs_prev_species_NY <- ab_vs_prev_species_NY %>% 
  mutate( Significance = case_when(trait %in% all_shallow_sig_NY$trait ~ "Significant traits", TRUE ~ "Non-significant traits"))
ab_vs_prev_species_NY <- ab_vs_prev_species_NY %>% 
  mutate( Common = case_when(trait %in% stats_species$trait ~ "Common traits", TRUE ~ "Non-common traits"))
#ab_vs_prev_species_NY$genus<-sapply(strsplit(ab_vs_prev_species_NY$trait, "_"), "[",3)
##MI
ab_vs_prev_species_MI=data.frame(trait=rownames(filtered_prev_species_MI),
                                 abundance=100*rowMeans(sweep(filtered_prev_species_MI,2,colSums(filtered_prev_species_MI),"/")),
                                 prevalence=100*rowSums(filtered_prev_species_MI>0)/ncol(filtered_prev_species_MI))
motch=match(all_shallow_MI$trait,ab_vs_prev_species_MI$trait)
ab_vs_prev_species_MI=ab_vs_prev_species_MI[motch,]
ab_vs_prev_species_MI <- ab_vs_prev_species_MI %>% 
  mutate( Significance = case_when(trait %in% all_shallow_sig_MI$trait ~ "Significant traits", TRUE ~ "Non-significant traits"))
ab_vs_prev_species_MI <- ab_vs_prev_species_MI %>% 
  mutate( Common = case_when(trait %in% stats_species$trait ~ "Common traits", TRUE ~ "Non-common traits"))
#ab_vs_prev_species_MI$genus<-sapply(strsplit(ab_vs_prev_species_MI$trait, "_"), "[",3)

####Volcano plots
#Species
max_ab=max(ab_vs_prev_species$abundance)
max_prev=max(ab_vs_prev_species$prevalence)
plot(x=1, type = "n", ylab="Average Abundance (%)", xlab="Average Prevalence (%)", ylim=c(0,max_ab), xlim=c(50,max_prev))
points(x = ab_vs_prev_species$prevalence[ab_vs_prev_species$Significance == "Non-significant traits"],
       y = ab_vs_prev_species$abundance[ab_vs_prev_species$Significance == "Non-significant traits"],
       pch = 19,
       col = "grey")
points(x = ab_vs_prev_species$prevalence[ab_vs_prev_species$Significance == "Significant traits"],
       y = ab_vs_prev_species$abundance[ab_vs_prev_species$Significance == "Significant traits"],
       pch = 19,
       col = "red")

#NY
max_ab=max(ab_vs_prev_species_NY$abundance)
max_prev=max(ab_vs_prev_species_NY$prevalence)
plot(x=1, type = "n", ylab="Average Abundance (%)", xlab="Average Prevalence (%)", ylim=c(0,max_ab), xlim=c(50,max_prev))
points(x = ab_vs_prev_species_NY$prevalence[ab_vs_prev_species_NY$Significance == "Non-significant traits"],
       y = ab_vs_prev_species_NY$abundance[ab_vs_prev_species_NY$Significance == "Non-significant traits"],
       pch = 19,
       col = "grey")
points(x = ab_vs_prev_species_NY$prevalence[ab_vs_prev_species_NY$Significance == "Significant traits"],
       y = ab_vs_prev_species_NY$abundance[ab_vs_prev_species_NY$Significance == "Significant traits"],
       pch = 19,
       col = "red")
points(x = ab_vs_prev_species_NY$prevalence[ab_vs_prev_species_NY$Common == "Common traits" & ab_vs_prev_species_NY$Significance == "Significant traits"],
       y = ab_vs_prev_species_NY$abundance[ab_vs_prev_species_NY$Common == "Common traits" & ab_vs_prev_species_NY$Significance == "Significant traits"],
       pch = 19,
       col = "blue")

#MI
max_ab=max(ab_vs_prev_species_MI$abundance)
max_prev=max(ab_vs_prev_species_MI$prevalence)
plot(x=1, type = "n", ylab="Average Abundance (%)", xlab="Average Prevalence (%)", ylim=c(0,max_ab), xlim=c(50,max_prev))
points(x = ab_vs_prev_species_MI$prevalence[ab_vs_prev_species_MI$Significance == "Non-significant traits"],
       y = ab_vs_prev_species_MI$abundance[ab_vs_prev_species_MI$Significance == "Non-significant traits"],
       pch = 19,
       col = "grey")
points(x = ab_vs_prev_species_MI$prevalence[ab_vs_prev_species_MI$Significance == "Significant traits"],
       y = ab_vs_prev_species_MI$abundance[ab_vs_prev_species_MI$Significance == "Significant traits"],
       pch = 19,
       col = "red")
points(x = ab_vs_prev_species_MI$prevalence[ab_vs_prev_species_MI$Common == "Common traits" & ab_vs_prev_species_MI$Significance == "Significant traits"],
       y = ab_vs_prev_species_MI$abundance[ab_vs_prev_species_MI$Common == "Common traits" & ab_vs_prev_species_MI$Significance == "Significant traits"],
       pch = 19,
       col = "blue")





#Genus
max_ab=max(ab_vs_prev_genus$abundance)
max_prev=max(ab_vs_prev_genus$prevalence)
plot(x=1, type = "n", ylab="Abundance (%)", xlab="Prevalence (%)", ylim=c(0,max_ab), xlim=c(50,max_prev))
points(x = ab_vs_prev_genus$prevalence[ab_vs_prev_genus$Significance == "Non-significant traits"],
       y = ab_vs_prev_genus$abundance[ab_vs_prev_genus$Significance == "Non-significant traits"],
       pch = 19,
       col = "grey")

#Family
max_ab=max(ab_vs_prev_family$abundance)
max_prev=max(ab_vs_prev_family$prevalence)
plot(x=1, type = "n", ylab="Abundance (%)", xlab="Prevalence (%)", ylim=c(0,max_ab), xlim=c(50,max_prev))
points(x = ab_vs_prev_family$prevalence[ab_vs_prev_family$Significance == "Non-significant traits"],
       y = ab_vs_prev_family$abundance[ab_vs_prev_family$Significance == "Non-significant traits"],
       pch = 19,
       col = "grey")
points(x = ab_vs_prev_family$prevalence[ab_vs_prev_family$Significance == "Significant traits"],
       y = ab_vs_prev_family$abundance[ab_vs_prev_family$Significance == "Significant traits"],
       pch = 19,
       col = "red")
#Phylum
max_ab=max(ab_vs_prev_phylum$abundance)
max_prev=max(ab_vs_prev_phylum$prevalence)
plot(x=1, type = "n", ylab="Abundance (%)", xlab="Prevalence (%)", ylim=c(0,max_ab), xlim=c(50,max_prev))
points(x = ab_vs_prev_phylum$prevalence[ab_vs_prev_phylum$Significance == "Non-significant traits"],
       y = ab_vs_prev_phylum$abundance[ab_vs_prev_phylum$Significance == "Non-significant traits"],
       pch = 19,
       col = "grey")
points(x = ab_vs_prev_phylum$prevalence[ab_vs_prev_phylum$Significance == "Significant traits"],
       y = ab_vs_prev_phylum$abundance[ab_vs_prev_phylum$Significance == "Significant traits"],
       pch = 19,
       col = "red")

###All
all_shallow=all_VCs_full
all_shallow=all_shallow[grepl("all",all_shallow$study),]
all_shallow <- all_shallow %>% 
  mutate( Common = case_when(trait %in% y ~ "Common", TRUE ~ "Not common"))
all_shallow=all_shallow[grepl("Common",all_shallow$Common),]




#########GWAS########

suggestive_l=5.8

#Load shallow taxonomic
load("~/P50/microbiome/output/microbiome_profiles/taxonomic_profile/gtdb_profiles/207/step2/abundance_mean/not_norm_comp_16S_207_E1/no_beta/H5/subset_spe_gtdb.RData")
all_shallow=rbind(filtered_clr_counts_objs[[3]],
                  filtered_clr_counts_objs[[6]],
                  filtered_clr_counts_objs[[9]],
                  filtered_clr_counts_objs[[12]],
                  filtered_clr_counts_objs[[15]],
                  filtered_clr_counts_objs[[18]])
pca_shallow=prcomp(all_shallow,center=TRUE,scale=TRUE)$sdev
total_var = sum(pca_shallow^2)
for (k in 1:length(pca_shallow)){
  if (sum((pca_shallow[1:k])^2) >= 0.95*total_var) {
    cat("K is: ",k, "\n")
    break  # Exit the loop
  }
}
sig_l_shallow=-log10(10^(-suggestive_l)/k)
load("/users/abaud/fmorillo/P50/microbiome/output/microbiome_profiles/taxonomic_profile/gtdb_profiles/207/step2/abundance_mean/not_norm_comp_16S_207_E1/no_beta/GWAS/All_QTLs.RData")
#all_QTLs <- all_QTLs %>%
#  arrange(chr, pos)
#all_QTLs <- all_QTLs %>%
#  group_by(chr) %>%
#  mutate(cumulative_position = cumsum(pos))
all_QTLs$snp=paste0(all_QTLs$measure," (chr",all_QTLs$chr,":",all_QTLs$pos,")")
sig_snps=all_QTLs$snp[grepl("Significant",all_QTLs$Significance)]
all_shallow=all_QTLs
all_shallow=all_shallow[grepl("ALL",all_shallow$subset_id),]
all_shallow$method="Shallow Shotgun"
# Significance level
all_shallow <- all_shallow %>% 
  mutate( Significance = case_when(logP > sig_l_shallow ~ "Significant", TRUE ~ "Not significant"))

#Load 16S
load("~/P50/microbiome/output/microbiome_profiles/taxonomic_profile/16S/207/step2/mean/filtered_for_shallow_207_E1_comparison_30k_5k_exc_corrected/H5/subset_spe_gtdb.RData")
all_16S=rbind(filtered_clr_counts_objs[[3]],
                  filtered_clr_counts_objs[[6]],
                  filtered_clr_counts_objs[[9]],
                  filtered_clr_counts_objs[[12]],
                  filtered_clr_counts_objs[[15]],
                  filtered_clr_counts_objs[[18]])
pca_shallow=prcomp(all_16S,center=TRUE,scale=TRUE)$sdev
total_var = sum(pca_shallow^2)
for (k in 1:length(pca_shallow)){
  if (sum((pca_shallow[1:k])^2) >= 0.95*total_var) {
    cat("K is: ",k, "\n")
    break  # Exit the loop
  }
}
sig_l_16S=-log10(10^(-suggestive_l)/k)
load("/users/abaud/fmorillo/P50/microbiome/output/microbiome_profiles/taxonomic_profile/16S/207/step2/mean/filtered_for_shallow_207_E1_comparison_30k_5k_exc_corrected/GWAS/All_QTLs.RData")
#all_QTLs <- all_QTLs %>%
#  arrange(chr, pos)
#all_QTLs <- all_QTLs %>%
#  group_by(chr) %>%
#  mutate(cumulative_position = cumsum(pos))
all_QTLs$snp=paste0(all_QTLs$measure," (chr",all_QTLs$chr,":",all_QTLs$pos,")")
sig_snps=all_QTLs$snp[grepl("Significant",all_QTLs$Significance)]
all_16S=all_QTLs
all_16S=all_16S[grepl("ALL",all_16S$subset_id),]
all_16S$method="16S"
# Significance level
all_16S <- all_16S %>% 
  mutate( Significance = case_when(logP > sig_l_16S ~ "Significant", TRUE ~ "Not significant"))


#Load shallow functional
load("/users/abaud/fmorillo/P50/microbiome/output/microbiome_profiles/functional_profile/humann_profiles/eggnog/H5/filtered_prev.RData")
all_func=filtered_clr_counts_objs[[3]]
matrix_with_otu_id <- all_func %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU_id")
#library(biomformat)
#x=make_biom(data = all_func)
#write_biom(make_biom(data = all_func), "/users/abaud/fmorillo/P50/microbiome/output/microbiome_profiles/functional_profile/humann_profiles/eggnog/H5/eggnog_matrix.biom")
#write.table(matrix_with_otu_id, file = "/users/abaud/fmorillo/P50/microbiome/output/microbiome_profiles/functional_profile/humann_profiles/eggnog/H5/eggnog_matrix.txt", quote=F, sep = "\t", row.names = F, col.names = T)
#write.csv(all_func, file = "/users/abaud/fmorillo/P50/microbiome/output/microbiome_profiles/functional_profile/humann_profiles/eggnog/H5/eggnog_matrix.csv", quote=F)
pca_func=prcomp(all_func,center=TRUE,scale=TRUE)$sdev
total_var = sum(pca_func^2)
for (k in 1:length(pca_func)){
  if (sum((pca_func[1:k])^2) >= 0.95*total_var) {
    cat("K is: ",k, "\n")
    break  # Exit the loop
  }
}
sig_l_func=-log10(10^(-suggestive_l)/k)
load("/users/abaud/fmorillo/P50/microbiome/output/microbiome_profiles/functional_profile/humann_profiles/eggnog/GWAS/All_QTLs.RData")
#all_QTLs <- all_QTLs %>%
#  arrange(chr, pos)
#all_QTLs <- all_QTLs %>%
#  group_by(chr) %>%
#  mutate(cumulative_position = cumsum(pos))
all_QTLs$snp=paste0(all_QTLs$measure," (chr",all_QTLs$chr,":",all_QTLs$pos,")")
sig_snps=all_QTLs$snp[grepl("Significant",all_QTLs$Significance)]
all_func=all_QTLs
all_func=all_func[grepl("ALL",all_func$subset_id),]
all_func$method="Shallow Shotgun - Functional"
# Significance level
all_func <- all_func %>% 
  mutate( Significance = case_when(logP > sig_l_func ~ "Significant", TRUE ~ "Not significant"))

#measure_16S=unique(all_16S$measure)
#measure_shallow=unique(all_shallow$measure)
#m=intersect(measure_16S,measure_shallow)

#filt_shallow <- all_shallow[all_shallow$measure %in% m, ]
#filt_16S <- all_16S[all_16S$measure %in% m, ]

max_16S=max(all_16S$logP)
max_shallow=max(all_shallow$logP)
max_func=max(all_func$logP)
# Make Manhattan plots
max_logP=max(max_16S,max_shallow,max_func)
#max_logP=max(all_func$logP)
pdf('/users/abaud/fmorillo/P50/microbiome/output/microbiome_profiles/taxonomic_profile/gwas_manhattan.pdf',bg='white', width = 12, height = 7)
#par(mfrow = c(3,1), cex.lab=1.5, cex.axis=1.6)
#color_groups <- c("Not significant","Significant")
#color_labels <- c("black", "firebrick3")
#Shallow - Taxonomic
sig_snps=all_shallow$snp[grepl("Significant",all_shallow$Significance)]
manhattan(x=all_shallow,
          chr="chr",
          bp="pos",
          p="pvalue",
          snp="snp",
          suggestiveline =suggestive_l,
          genomewideline =sig_l_shallow,
          annotatePval = 10^(-(sig_l_shallow)),
          annotateTop = T,
          highlight=sig_snps,
          ylim=c(3.5,max(max_logP+2)),
          main = "Method: Shallow Shotgun Taxonomic (Study: ALL)",
          cex.main = 1.5,
          chrlabs=as.character(c(1:20)))
#myplot<-ggplot(all_shallow, aes(x=as.numeric(chr), y=as.numeric(logP))) + 
#  geom_point(aes(color = Significance), size = 1.5) +
#  xlab("Chromosomes (Shallow Shotgun - Taxonomic Profile) ") +
#  ylab(expression("-log"[10]*"P-Value")) +
#  theme_bw() +
#  theme(legend.position="bottom",
#        legend.key.size = unit(2, "lines"),
#        legend.text = element_text(size = 13),
#        legend.title = element_text(size = 13),
#        axis.line = element_line(),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        panel.border = element_blank(),
#        panel.background = element_blank(),
#        axis.text = element_text(size = 15),
#        axis.title = element_text(size = 15)) +
#  coord_cartesian(xlim = c(1, 20), ylim = c(4, max_logP+2)) +
#  scale_y_continuous(breaks = seq(4, max_logP+2, by = 1)) +
#  scale_x_continuous(breaks = seq(1, 20, by = 1)) +
#  scale_color_manual(values = color_labels,
#                     labels = color_groups,
#                     limits = color_groups) +
#  labs(color = "")
#print(myplot)
#16S
sig_snps=all_16S$snp[grepl("Significant",all_16S$Significance)]
manhattan(x=all_16S,
          chr="chr",
          bp="pos",
          p="pvalue",
          snp="snp",
          suggestiveline =suggestive_l,
          genomewideline =sig_l_16S,
          annotatePval = 10^(-(sig_l_16S)),
          annotateTop = T,
          highlight=sig_snps,
          ylim=c(3.5,max(max_logP+2)),
          main = "Method: 16S Taxonomic (Study: ALL)",
          cex.main = 1.5,
          chrlabs=as.character(c(1:20)))
#myplot<-ggplot(all_16S, aes(x=as.numeric(chr), y=as.numeric(logP))) + 
#  geom_point(aes(color = Significance), size = 1.5) +
#  xlab("Chromosomes (16S - Taxonomic Profile) ") +
#  ylab(expression("-log"[10]*"P-Value")) +
#  theme_bw() +
#  theme(legend.position="bottom",
#        legend.key.size = unit(2, "lines"),
#        legend.text = element_text(size = 13),
#        legend.title = element_text(size = 13),
#        axis.line = element_line(),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        panel.border = element_blank(),
#        panel.background = element_blank(),
#        axis.text = element_text(size = 15),
#        axis.title = element_text(size = 15)) +
#  coord_cartesian(xlim = c(1, 20), ylim = c(4, max_logP+2)) +
#  scale_y_continuous(breaks = seq(4, max_logP+2, by = 1)) +
#  scale_x_continuous(breaks = seq(1, 20, by = 1)) +
#  scale_color_manual(values = color_labels,
#                     labels = color_groups,
#                     limits = color_groups) +
#  labs(color = "")
#print(myplot)
#Shallow Functional
sig_snps=all_func$snp[grepl("Significant",all_func$Significance)]
manhattan(x=all_func,
          chr="chr",
          bp="pos",
          p="pvalue",
          snp="snp",
          suggestiveline =suggestive_l,
          genomewideline =sig_l_func,
          annotatePval = 10^(-(sig_l_func)),
          annotateTop = T,
          highlight=sig_snps,
          ylim=c(3.5,max(max_logP+2)),
          main = "Method: Shallow Shotgun Functional (Study: ALL)",
          cex.main = 1.5,
          chrlabs=as.character(c(1:20)))
#myplot<-ggplot(all_func, aes(x=as.numeric(chr), y=as.numeric(logP))) + 
#  geom_point(aes(color = Significance), size = 1.5) +
#  xlab("Chromosomes (Shallow Shotgun - Functional Profile) ") +
#  ylab(expression("-log"[10]*"P-Value")) +
#  theme_bw() +
#  theme(legend.position="bottom",
#        legend.key.size = unit(2, "lines"),
#        legend.text = element_text(size = 13),
#        legend.title = element_text(size = 13),
#        axis.line = element_line(),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        panel.border = element_blank(),
#        panel.background = element_blank(),
#        axis.text = element_text(size = 15),
#        axis.title = element_text(size = 15)) +
#  coord_cartesian(xlim = c(1, 20), ylim = c(4, max_logP+2)) +
#  scale_y_continuous(breaks = seq(4, max_logP+2, by = 1)) +
#  scale_x_continuous(breaks = seq(1, 20, by = 1)) +
#  scale_color_manual(values = color_labels,
#                     labels = color_groups,
#                     limits = color_groups) +
#  labs(color = "")
#print(myplot)
dev.off()


#### Pie chart

#"s__Paraprevotella_sp900548345" (Human gut) "s__Paraprevotella_sp900760455" (Human gut) "s__Paraprevotella_sp910584955" (Mouse gut)

#Load shallow taxonomic
load("~/P50/microbiome/output/microbiome_profiles/taxonomic_profile/gtdb_profiles/207/step2/abundance_mean/not_norm_comp_16S_207_E1/no_beta/H5/filtered_prev.RData")
shallow_species=filtered_prev_objs[[18]]
shallow_species <- sweep(shallow_species,2,colSums(shallow_species),"/")
shallow_species <-shallow_species[grepl("s__Paraprevotella",rownames(shallow_species)),]
shallow_species <-shallow_species*100
shallow_paraprev=data.frame(species=rownames(shallow_species))
shallow_paraprev$rel_abundance[1]=median(shallow_species[1,])
shallow_paraprev$rel_abundance[2]=median(shallow_species[2,])
shallow_paraprev$rel_abundance[3]=median(shallow_species[3,])
# Create a pie chart with custom colors
colors <- c("darkgreen", "lightgreen", "red")
pie(shallow_paraprev$rel_abundance, col = colors, labels = round(shallow_paraprev$rel_abundance,3))


######

#Load shallow taxonomic
load("~/P50/microbiome/output/microbiome_profiles/taxonomic_profile/gtdb_profiles/207/step2/abundance_mean/not_norm_comp_16S_207_E1/no_beta/H5/subset_spe_gtdb.RData")
shallow_species=filtered_clr_counts_objs[[18]]
shallow_paraprevotella=shallow_species[grepl("s__Paraprevotella",rownames(shallow_species)),]

#boxplot(shallow_paraprevotella[1,],
#        shallow_paraprevotella[2,],
#        shallow_paraprevotella[3,],col = "transparent",
#        names = c("s__Paraprevotella_sp900548345", "s__Paraprevotella_sp900760455", "s__Paraprevotella_sp910584955"))
#rownames(shallow_paraprevotella)

p1=data.frame(species=rownames(shallow_paraprevotella)[1],rel_abundance=shallow_paraprevotella[1,])
p2=data.frame(species=rownames(shallow_paraprevotella)[2],rel_abundance=shallow_paraprevotella[2,])
p3=data.frame(species=rownames(shallow_paraprevotella)[3],rel_abundance=shallow_paraprevotella[3,])
all_para=rbind(p1,p2,p3)

pdf('/users/abaud/fmorillo/P50/microbiome/output/microbiome_profiles/taxonomic_profile/paraprevotella_species.pdf',bg='white')
order=c("s__Paraprevotella_sp900548345", "s__Paraprevotella_sp900760455", "s__Paraprevotella_sp910584955")
myplot <- ggplot(all_para, aes(x=factor(species,level = order), y = rel_abundance, color=species)) + 
  geom_violin(width = 0.5) +  # Use geom_violin to make it narrow
  #geom_boxplot(outlier.shape = NA, fill = "transparent", color=c("red","black","black"), width = 0.2) +  # Add boxplot for reference
  xlab("") +
  ylab("Relative abundances (After CLR)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = NULL) +
  scale_color_manual(values=c("darkgrey", "darkgrey","red")) +
  ggtitle("Paraprevotella species")
print(myplot)

dev.off()


####### Correlations between 16S and Shallow (taxonomy)
#invrank= function(row) {qnorm((rank(row,na.last="keep",ties.method = 'random')-0.5)/sum(!is.na(row)))}
ranks=c("Phylum","Class","Order","Family","Genus","Species")

#### Method - CLR

#Load shallow taxonomic
load("~/P50/microbiome/output/microbiome_profiles/taxonomic_profile/gtdb_profiles/207/step2/abundance_mean/not_norm_comp_16S_207_E1/no_beta/H5/subset_spe_gtdb.RData")
all_shallow_Phylum=filtered_clr_counts_objs[[3]]
all_shallow_Class=filtered_clr_counts_objs[[6]]
all_shallow_Order=filtered_clr_counts_objs[[9]]
all_shallow_Family=filtered_clr_counts_objs[[12]]
all_shallow_Genus=filtered_clr_counts_objs[[15]]
all_shallow_Species=filtered_clr_counts_objs[[18]]

#Load 16S taxonomic
load("~/P50/microbiome/output/microbiome_profiles/taxonomic_profile/16S/207/step2/mean/filtered_for_shallow_207_E1_comparison_30k_5k_exc_corrected/H5/subset_spe_gtdb.RData")
all_16S_Phylum=filtered_clr_counts_objs[[3]]
all_16S_Class=filtered_clr_counts_objs[[6]]
all_16S_Order=filtered_clr_counts_objs[[9]]
all_16S_Family=filtered_clr_counts_objs[[12]]
all_16S_Genus=filtered_clr_counts_objs[[15]]
all_16S_Species=filtered_clr_counts_objs[[18]]

correlation_clr=data.frame(rank=ranks,r2_pearson=0,pvalue_pearson=0,r2_spearman=0,pvalue_spearman=0,method="CLR")
for (rank in ranks){
  tab_shallow=get(paste("all_shallow",rank,sep="_"))
  tab_16S=get(paste("all_16S",rank,sep="_"))
  m=intersect(rownames(tab_16S),rownames(tab_shallow))
  filt_shallow=tab_shallow[m,]
  filt_16S=tab_16S[m,]
  #Transform data to make it parametric
  #filt_shallow = t(apply(filt_shallow, FUN = invrank, MAR = 1))
  #filt_16S = t(apply(filt_16S, FUN = invrank, MAR = 1))
  plot(filt_16S,filt_shallow,xlab=paste0("CLR abundances 16S - ",rank),ylab=paste0("CLR Abundances Shallow Shotgun - ",rank))
  correlation_clr$r2_pearson[which(ranks==rank)]=cor.test(filt_16S,filt_shallow,method="pearson")$estimate*100
  correlation_clr$pvalue_pearson[which(ranks==rank)]=cor.test(filt_16S,filt_shallow,method="pearson")$p.value
  correlation_clr$r2_spearman[which(ranks==rank)]=cor.test(filt_16S,filt_shallow,method="spearman")$estimate*100
  correlation_clr$pvalue_spearman[which(ranks==rank)]=cor.test(filt_16S,filt_shallow,method="spearman")$p.value
  print(rank)
}


#### Method - Relative abundance

#Load shallow taxonomic
load("~/P50/microbiome/output/microbiome_profiles/taxonomic_profile/gtdb_profiles/207/step2/abundance_mean/not_norm_comp_16S_207_E1/no_beta/H5/filtered_prev.RData")
all_shallow_Phylum=filtered_prev_objs[[3]]
all_shallow_Phylum <- sweep(all_shallow_Phylum,2,colSums(all_shallow_Phylum),"/")
all_shallow_Class=filtered_prev_objs[[6]]
all_shallow_Class <- sweep(all_shallow_Class,2,colSums(all_shallow_Class),"/")
all_shallow_Order=filtered_prev_objs[[9]]
all_shallow_Order <- sweep(all_shallow_Order,2,colSums(all_shallow_Order),"/")
all_shallow_Family=filtered_prev_objs[[12]]
all_shallow_Family <- sweep(all_shallow_Family,2,colSums(all_shallow_Family),"/")
all_shallow_Genus=filtered_prev_objs[[15]]
all_shallow_Genus <- sweep(all_shallow_Genus,2,colSums(all_shallow_Genus),"/")
all_shallow_Species=filtered_prev_objs[[18]]
all_shallow_Species <- sweep(all_shallow_Species,2,colSums(all_shallow_Species),"/")

#Load 16S taxonomic
load("~/P50/microbiome/output/microbiome_profiles/taxonomic_profile/16S/207/step2/mean/filtered_for_shallow_207_E1_comparison_30k_5k_exc_corrected/H5/filtered_prev.RData")
all_16S_Phylum=filtered_prev_objs[[3]]
all_16S_Phylum <- sweep(all_16S_Phylum,2,colSums(all_16S_Phylum),"/")
all_16S_Class=filtered_prev_objs[[6]]
all_16S_Class <- sweep(all_16S_Class,2,colSums(all_16S_Class),"/")
all_16S_Order=filtered_prev_objs[[9]]
all_16S_Order <- sweep(all_16S_Order,2,colSums(all_16S_Order),"/")
all_16S_Family=filtered_prev_objs[[12]]
all_16S_Family <- sweep(all_16S_Family,2,colSums(all_16S_Family),"/")
all_16S_Genus=filtered_prev_objs[[15]]
all_16S_Genus <- sweep(all_16S_Genus,2,colSums(all_16S_Genus),"/")
all_16S_Species=filtered_prev_objs[[18]]
all_16S_Species <- sweep(all_16S_Species,2,colSums(all_16S_Species),"/")

correlation_relabandance=data.frame(rank=ranks,r2_pearson=0,pvalue_pearson=0,r2_spearman=0,pvalue_spearman=0,method="Relative Abundance")
for (rank in ranks){
  tab_shallow=get(paste("all_shallow",rank,sep="_"))
  tab_16S=get(paste("all_16S",rank,sep="_"))
  m=intersect(rownames(tab_16S),rownames(tab_shallow))
  filt_shallow=tab_shallow[m,]
  filt_16S=tab_16S[m,]
  #Transform data to make it parametric
  #filt_shallow = t(apply(filt_shallow, FUN = invrank, MAR = 1))
  #filt_16S = t(apply(filt_16S, FUN = invrank, MAR = 1))
  plot(filt_16S,filt_shallow,xlab=paste0("Relative abundances 16S - ",rank),ylab=paste0("Relative Abundances Shallow Shotgun - ",rank))
  correlation_relabandance$r2_pearson[which(ranks==rank)]=cor.test(filt_16S,filt_shallow,method="pearson")$estimate*100
  correlation_relabandance$pvalue_pearson[which(ranks==rank)]=cor.test(filt_16S,filt_shallow,method="pearson")$p.value
  correlation_relabandance$r2_spearman[which(ranks==rank)]=cor.test(filt_16S,filt_shallow,method="spearman")$estimate*100
  correlation_relabandance$pvalue_spearman[which(ranks==rank)]=cor.test(filt_16S,filt_shallow,method="spearman")$p.value
  print(rank)
}

####
library(ape)
immeDB_annot=read.gff("/users/abaud/data/secondary/indexes/immeDB/Data2_Annotation_of_genes_on_MGEs.gff3", na.strings = c(".", "?"), GFF3 = TRUE)

#################


# Compare proportions

# 100% 16S
load("~/NextFlow/Microbiome_Profiling/step2/Output_16S_Full_Harm/Heritability/heritability.RData")
herit_16S_100=all_VCs_full
herit_16S_100=herit_16S_100[grepl("ALL",herit_16S_100$subset_id),]
herit_16S_100$method="16S_100%"

# 75% 16S
load("~/NextFlow/Microbiome_Profiling/step2/Output_16S_75/Heritability/heritability.RData")
herit_16S_75=all_VCs_full
herit_16S_75=herit_16S_75[grepl("ALL",herit_16S_75$subset_id),]
herit_16S_75$method="16S_75%"

# 50% 16S
load("~/NextFlow/Microbiome_Profiling/step2/Output_16S_50/Heritability/heritability.RData")
herit_16S_50=all_VCs_full
herit_16S_50=herit_16S_50[grepl("ALL",herit_16S_50$subset_id),]
herit_16S_50$method="16S_50%"

# 25% 16S
load("~/NextFlow/Microbiome_Profiling/step2/Output_16S_25/Heritability/heritability.RData")
herit_16S_25=all_VCs_full
herit_16S_25=herit_16S_25[grepl("ALL",herit_16S_25$subset_id),]
herit_16S_25$method="16S_25%"

# 0% 16S
load("~/NextFlow/Microbiome_Profiling/step2/Output_Shallow_Full_Harm/Heritability/heritability.RData")
herit_16S_0=all_VCs_full
herit_16S_0=herit_16S_0[grepl("ALL",herit_16S_0$subset_id),]
herit_16S_0$method="16S_0%"

all_methods_common=rbind(herit_16S_100,
                         herit_16S_75,
                         herit_16S_50,
                         herit_16S_25,
                         herit_16S_0)

all_methods_common=all_methods_common[!grepl("g__Limivicinus",all_methods_common$trait),]

order=c("16S_100%","16S_75%","16S_50%","16S_25%","16S_0%")
myplot<-ggplot(all_methods_common, aes(x=factor(method,level = order), y=var_Ad*100)) + 
  geom_boxplot() +
  xlab("Method") +
  ylab("Heritability (%)") +
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
  ylim(0,30) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  geom_signif(comparisons = list(c("16S_100%","16S_75%")),
                                 #c("16S_75%","16S_50%"),
                                 #c("16S_50%","16S_25%"),
                                 #c("16S_25%","16S_0%"),
                                 #c("16S_100%","16S_0%")
              map_signif_level=TRUE)
print(myplot)



###################

load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Shallow_Sample_Harm/Heritability/heritability.RData")

subset_id="ALL" ### Remove
sig_f=0.1
sig_b=0.05
all_VCs_full=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),] ### Remove

all_VCs_full$bonf_pvalue_DGE = all_VCs_full$pvalue_DGE*dim(all_VCs_full)[1]
all_VCs_full$qvalue_DGE = qvalue(all_VCs_full$pvalue_DGE)$qvalue
all_VCs_full$qp_diff = all_VCs_full$qvalue_DGE - all_VCs_full$pvalue_DGE
all_VCs_full <- all_VCs_full %>% 
  mutate( Significance = case_when(
    qvalue_DGE < sig_f & bonf_pvalue_DGE < sig_b ~ "Heritable traits (Bonferroni < 5%)",
    qvalue_DGE < sig_f & bonf_pvalue_DGE >= sig_b ~ "Heritable traits (FDR < 10%)",
    TRUE ~ "Non-heritable traits"))
all_VCs_full = all_VCs_full[order(all_VCs_full$logP, decreasing = T),]
all_VCs_full$qvalue_DGE = round(all_VCs_full$qvalue_DGE, digits = 2)

all_shallow=all_VCs_full
all_shallow=all_shallow[grepl("ALL",all_shallow$subset_id),]
all_shallow$method="Shallow Shotgun"


load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_16S_Sample_Harm/Heritability/heritability.RData")

subset_id="ALL" ### Remove
sig_f=0.1
sig_b=0.05
all_VCs_full=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),] ### Remove

all_VCs_full$bonf_pvalue_DGE = all_VCs_full$pvalue_DGE*dim(all_VCs_full)[1]
all_VCs_full$qvalue_DGE = qvalue(all_VCs_full$pvalue_DGE)$qvalue
all_VCs_full$qp_diff = all_VCs_full$qvalue_DGE - all_VCs_full$pvalue_DGE
all_VCs_full <- all_VCs_full %>% 
  mutate( Significance = case_when(
    qvalue_DGE < sig_f & bonf_pvalue_DGE < sig_b ~ "Heritable traits (Bonferroni < 5%)",
    qvalue_DGE < sig_f & bonf_pvalue_DGE >= sig_b ~ "Heritable traits (FDR < 10%)",
    TRUE ~ "Non-heritable traits"))
all_VCs_full = all_VCs_full[order(all_VCs_full$logP, decreasing = T),]
all_VCs_full$qvalue_DGE = round(all_VCs_full$qvalue_DGE, digits = 2)

all_16S=all_VCs_full
all_16S=all_16S[grepl("ALL",all_16S$subset_id),]
all_16S$method="16S"

all_methods=rbind(all_shallow,all_16S)
all_methods$heritability=100*all_methods$var_Ad

# Specify the desired order for ranks and Significance
desired_rank_order <- c("Phylum","Class","Order","Family","Genus","Species","ASV")
desired_significance_order <- c('Non-heritable traits', 'Heritable traits (FDR < 10%)', 'Heritable traits (Bonferroni < 5%)')

# Convert the 'rank' column to a factor with the specified order
all_methods$rank <- factor(all_methods$rank, levels = desired_rank_order)

# Convert the 'Significance' column to a factor with the specified order
all_methods$Significance <- factor(all_methods$Significance, levels = desired_significance_order)

# Summarize the data: count occurrences of each Significance category per method and rank
summary_df <- all_methods %>%
  group_by(rank, method, Significance) %>%
  summarize(count = n()) %>%
  ungroup()

summary_df=na.omit(summary_df)

# Compute bar top positions by grouping and summing
top_labels <- summary_df %>%
  group_by(method, rank) %>%
  summarise(total_count = sum(count), .groups = "drop") %>%
  inner_join(
    summary_df %>% filter(Significance == "Heritable traits (FDR < 10%)"),
    by = c("method", "rank")
  )

pdf("~/paper_figures/16s_vs_SS_sig_counts_norm_maternal.pdf")
myplot<-ggplot(summary_df, aes(x = factor(rank), y = count, fill = Significance)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  geom_text(data = top_labels,
            aes(x = factor(rank), y = total_count, label = count),
            color="darkblue",
            size = 5,
            vjust = -0.2,
            inherit.aes = FALSE) +
  facet_wrap(~method, strip.position = "bottom") +
  scale_fill_manual(values = c("Heritable traits (Bonferroni < 5%)" = "blue",
                               "Heritable traits (FDR < 10%)" = "darkblue", 
                               "Non-heritable traits" = "lightblue")) + 
  labs(x = "", y = "Number of traits", fill = "", title = "Number of heritable traits by rank and method") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    axis.line = element_line(color = "grey"),   
    panel.grid.major = element_blank(),         
    panel.grid.minor = element_blank(),         
    axis.text.x = element_text(size = 12, angle = 60, hjust = 1), 
    axis.text.y = element_text(size = 12), 
    axis.title.y = element_text(size = 14),
    strip.text = element_text(face = "bold", size = 14),  
    strip.placement = "outside",                
    panel.spacing = unit(0.1, "lines"),
    plot.title = element_text(hjust = 0.5)
    )
print(myplot)
# Zoom in on genus

plot_data <- summary_df %>%
  filter(grepl("Genus", rank))

# Compute bar top positions by grouping and summing
top_labels <- plot_data %>%
  group_by(method, rank) %>%
  summarise(total_count = sum(count), .groups = "drop") %>%
  inner_join(
    plot_data %>% filter(Significance == "Heritable traits (FDR < 10%)"),
    by = c("method", "rank")
  )

myplot<-ggplot(summary_df[grepl("Genus",summary_df$rank),], aes(x = factor(rank), y = count, fill = Significance)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  geom_text(data = top_labels,
            aes(x = factor(rank), y = total_count, label = count),
            color="darkblue",
            size = 5,
            vjust = -0.2,
            inherit.aes = FALSE) +
  facet_wrap(~method, strip.position = "bottom") +
  scale_fill_manual(values = c("Heritable traits (Bonferroni < 5%)" = "blue",
                               "Heritable traits (FDR < 10%)" = "darkblue", 
                               "Non-heritable traits" = "lightblue")) + 
  labs(x = "", y = "Number of traits", fill = "", title = "Number of heritable genera by method") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    axis.line = element_line(color = "grey"),   
    panel.grid.major = element_blank(),         
    panel.grid.minor = element_blank(),         
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, angle = 60, hjust = 1), 
    axis.title.y = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 14),
    strip.placement = "outside",                
    panel.spacing = unit(0.1, "lines"),
    plot.title = element_text(hjust = 0.5)
  )
print(myplot)
dev.off()


#########

load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Shallow_Full_Harm/Matrix_Processing/residuals_qned_counts.RData")
residuals_SS=residuals_qned_counts_objs
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2//Output_16S_Full_Harm/Matrix_Processing/residuals_qned_counts.RData")
residuals_16S=residuals_qned_counts_objs

ranks=c("Phylum","Class","Order","Family","Genus","Species")
names(ranks)=c(3,6,9,12,15,18)
pdf('~/paper_figures/clr_16s_shallow.pdf')
for (rank in ranks){
  number=names(ranks)[which(ranks==rank)]
  resids_ss=residuals_SS[[as.numeric(number)]]
  resids_16S=residuals_16S[[as.numeric(number)]]
  flat_A <- as.vector(resids_ss)
  flat_B <- as.vector(resids_16S)
  plot(flat_A,flat_B,xlab="Relative abundances (CLR residuals): Shallow Shotgun",ylab="Relative abundances (CLR residuals): 16S")
  cor_result=cor.test(resids_ss,resids_16S,method="pearson")
  r_value <- cor_result$estimate
  p_value <- cor_result$p.value
  title(main = paste0("Rank: ",rank), line = 2.5, cex.main = 1.5)
  legend("bottomright", legend = c(paste("rho =", round(r_value, 3)), "P < 0.01"), col = c("transparent", "transparent"), pch = 1, bg="white")
}
dev.off()

#mtext(paste("rho =", round(r_value, 4)), side = 3, line = 3, cex = 0.9, col = "black")
#mtext(paste("p-value =", format(p_value, digits = 4, scientific = TRUE)), side = 3, line = 4, cex = 0.9, col = "black")


cor.test(resids_ss,resids_16S,method="pearson")
cor.test(resids_ss,resids_16S,method="spearman")
