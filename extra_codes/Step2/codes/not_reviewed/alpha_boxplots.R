covariate="Study"
subset_id="ALL"
order=c("MI","NY")

pdf("/users/abaud/fmorillo/paper_figures/alpha_boxplots.pdf")

# Simpson
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Cluster_Alpha/Diversity/alpha_diversity.RData")
alpha_table=alpha_table_objs[[8]]
myplot<-ggplot(alpha_table, aes(x=factor(!!as.name(covariate),level=order), y=qD)) +
  geom_boxplot() +
  xlab(paste0("Subset: ",subset_id," / Covariate: ",covariate)) +
  ylab("Asymptotic Alpha Diversity (Shannon Index)") +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12))+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  geom_signif(comparisons = list(unique(alpha_table[[covariate]])),
              map_signif_level = TRUE)
print(myplot)

# Rao's
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Cluster_Final/Diversity/alpha_diversity.RData")
alpha_table=alpha_table_objs[[3]]
myplot<-ggplot(alpha_table, aes(x=factor(!!as.name(covariate),level=order), y=qD)) +
  geom_boxplot() +
  xlab(paste0("Subset: ",subset_id," / Covariate: ",covariate)) +
  ylab("Asymptotic Alpha Diversity (Rao's Quadratic Index)") +
  theme_bw() +
  theme(legend.position="bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12))+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  geom_signif(comparisons = list(unique(alpha_table[[covariate]])),
              map_signif_level = TRUE)
print(myplot)

dev.off()