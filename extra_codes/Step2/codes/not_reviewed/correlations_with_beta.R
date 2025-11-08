load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Cluster_Cage/Heritability/heritability.RData")
filtered_VCs=all_VCs_full[grepl("ALL",all_VCs_full$subset_id),]

load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Cluster_Cage/Cluster_Analyses/residuals_qned_counts.RData")
all_traits_ALL=residuals_qned_counts_objs[[18]]

##### Beta

load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Cluster_Cage/Matrix_Processing/beta_diversity.RData")
beta=gp.dist_objs[[6]]
genome.pcoa <- ape::pcoa(beta)
genome.pc <- data.frame(genome.pcoa$vectors)
#genome.pc[1:5,1:5]
#all_traits_ALL[1:5,1:5]

# Calculate correlation between CLR counts and beta diversity
correlations=data.frame()
for(i in 1:length(rownames(all_traits_ALL))){
  test_spearman_alpha_pd = cor.test(all_traits_ALL[i,],genome.pc$Axis.1, method = "spearman", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_ALL)[i],cor_beta_pd=test_spearman_alpha_pd$estimate,pval=-log(test_spearman_alpha_pd$p.value))
  correlations=rbind(correlations,new_row)
}

###### Alpha

load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Cluster_Cage/Matrix_Processing/alpha_diversity.RData")
alpha=alpha_table_objs[[3]]

# Calculate correlation between CLR counts and alpha diversity
correlations=data.frame()
for(i in 1:length(rownames(all_traits_ALL))){
  test_spearman_alpha_pd = cor.test(all_traits_ALL[i,],alpha$qD, method = "spearman", exact=FALSE)
  new_row=data.frame(trait=rownames(all_traits_ALL)[i],cor_alpha_pd=test_spearman_alpha_pd$estimate,pval=-log(test_spearman_alpha_pd$p.value))
  correlations=rbind(correlations,new_row)
}

#######

correlations$heritability=100*filtered_VCs$var_Ad[match(correlations$trait,filtered_VCs$trait)]
correlations$Phylum=taxonomy$Phylum[match(correlations$trait,taxonomy$Species)]
correlations$Family=taxonomy$Family[match(correlations$trait,taxonomy$Species)]
correlations$Genus=taxonomy$Genus[match(correlations$trait,taxonomy$Species)]


correlations_plot=correlations
correlations_plot=correlations[grepl("Firmicutes_A",correlations$Phylum) | grepl("Bacteroidota",correlations$Phylum),]
correlations_plot=correlations[grepl("Lachnospiraceae",correlations$Family) | grepl("Bacteroidaceae",correlations$Family),]
correlations_plot=correlations[grepl("CAG-95",correlations$Genus) | grepl("Prevotella",correlations$Genus),]

myplot<-ggplot(correlations_plot, aes(x=reorder(Genus,cor_beta_pd,FUN = median), y=cor_beta_pd)) +
  geom_boxplot() +
  xlab("") +
  ylab("Correlation with beta diversity (PCoA1)") +
  #ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.text=element_text(size=10),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 45))
  #geom_signif(comparisons = list(c("Non-significant","Significant")),
  #            map_signif_level = TRUE)
print(myplot)

myplot<-ggplot(correlations_plot, aes(x=reorder(Genus,heritability,FUN = median), y=heritability)) +
  geom_boxplot() +
  xlab("") +
  ylab("Heritability (%)") +
  #ggtitle(paste0("Subset: ",subset_id)) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.text=element_text(size=10),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(guide = guide_axis(angle = 45))
#geom_signif(comparisons = list(c("Non-significant","Significant")),
#            map_signif_level = TRUE)
print(myplot)



plot(correlations$cor_beta_pd,correlations$heritability,xlab="Correlation with beta diversity (R2)",ylab="Heritability")
cor.test(correlations$cor_beta_pd,correlations$heritability,method="pearson")

plot(correlations$pval,correlations$heritability,xlab="Correlation with beta diversity (-logP)",ylab="Heritability")
cor.test(correlations$pval,correlations$heritability,method="pearson")


##########

cor.test(all_traits_ALL[which(rownames(all_traits_ALL)=="s__CAG-95_sp014803415"),],all_traits_ALL[which(rownames(all_traits_ALL)=="s__Prevotella_sp002251365"),], method = "pearson", exact=FALSE)
plot(all_traits_ALL[which(rownames(all_traits_ALL)=="s__CAG-95_sp014803415"),],all_traits_ALL[which(rownames(all_traits_ALL)=="s__Prevotella_sp002251365"),],xlab="CAG-95 sp014803415",ylab="Prevotella sp002251365")                                                  
abline(lm(all_traits_ALL[which(rownames(all_traits_ALL)=="s__Prevotella_sp002251365"),] ~ all_traits_ALL[which(rownames(all_traits_ALL)=="s__CAG-95_sp014803415"),]), col = "red", lwd = 2)

cor.test(all_traits_ALL[which(rownames(all_traits_ALL)=="s__Acetatifactor_sp910578995"),],all_traits_ALL[which(rownames(all_traits_ALL)=="s__Prevotella_sp900769055"),], method = "spearman", exact=FALSE)
cor.test(all_traits_ALL[which(rownames(all_traits_ALL)=="s__Acetatifactor_sp910578995"),],all_traits_ALL[which(rownames(all_traits_ALL)=="s__Prevotella_sp900545525"),], method = "spearman", exact=FALSE)

correlations_acet=data.frame()
all_traits_prev=all_traits_ALL[rownames(all_traits_ALL) %in% x$trait,]
for(i in 1:length(rownames(all_traits_prev))){
  test_spearman_alpha_pd = cor.test(all_traits_prev[i,],all_traits_ALL[which(rownames(all_traits_ALL)=="s__Acetatifactor_sp910578995"),], method = "spearman", exact=FALSE)
  new_row=data.frame(Species=rownames(all_traits_prev)[i],cor_beta_pd=test_spearman_alpha_pd$estimate,pval=-log(test_spearman_alpha_pd$p.value))
  correlations_acet=rbind(correlations_acet,new_row)
}
