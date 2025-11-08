library(Matrix)
library(vegan)
library(ggplot2)
library(ggsignif)

subset_id="ALL"
number=1
#if (subset_id=="MI"){
#  number=number-2
#} else {
#  if (subset_id=="NY"){
#    number=number-1
#  }
#}

# Load Alpha
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Tax_Mat_Net/Diversity/residuals_qned_counts_alpha.RData")
alpha=residuals_qned_counts_alpha_objs[[number]]
rownames(alpha)
# [1] "alpha__TD_q0" "alpha__TD_q1" "alpha__TD_q2" "alpha__PD_q0" "alpha__PD_q1" "alpha__PD_q2"
alpha=alpha[6,]

# Load beta
load("/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Tax_Mat_Net/Diversity/residuals_qned_counts_beta.RData")
beta=residuals_qned_counts_beta_objs[[1]]
beta=as.matrix(beta)
rownames(beta)
beta_df=data.frame(Axis.1=beta[5,],Axis.2=beta[6,])
beta=beta[5,]
#beta <- beta - min(beta, na.rm = TRUE)
#rownames(beta)
# [1] "beta__TD_dispersion" "beta__TD_PC1" "beta__TD_PC2"  
#beta <- as.data.frame(beta[2,,drop=F])
gp.dist=gp.dist_objs[[2]]

# Load enterotypes
load("/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Tax_Mat_Net/Enterotypes/enterotypes.RData")
ent_samples=enterotypes_samples_objs[[1]]
ent_samples$alpha=alpha[match(ent_samples$sample,names(alpha))]
ent_samples$beta=beta[match(ent_samples$sample,names(beta))]
#ent_samples$beta=as.numeric(ent_samples$beta)
#any(is.na(ent_samples))

beta_df=beta_df[rownames(beta_df) %in% ent_samples$sample,]
beta_df=beta_df[ent_samples$sample,]

#hist(ent_samples$alpha)
#hist(ent_samples$beta)

#t_test_alpha <- t.test(alpha ~ Enterotype, data = ent_samples)
#t_test_beta <- t.test(beta ~ Enterotype, data = ent_samples)
#
#t_test_alpha
#t_test_beta

formula <- as.formula("gp.dist ~ Enterotype")
permanova.beta <- adonis2(formula, data = ent_samples, permutations = 999, seed = 123)
formula <- as.formula("alpha ~ Enterotype")
anova.alpha <- summary(aov(formula, data = ent_samples))[[1]]
anova.alpha$R2=anova.alpha$`Sum Sq`/sum(anova.alpha$`Sum Sq`)

# Boxplot Alpha
order=c("enterotype__1", "enterotype__2")
myplot<-ggplot(ent_samples, aes(x=factor(Enterotype, levels = order), y=alpha)) +
  geom_boxplot() +
  xlab("") +
  ylab("Rao's quadratic index (residuals)") +
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
  geom_signif(comparisons = list(c("enterotype__1", "enterotype__2")),
              map_signif_level = TRUE)
print(myplot)

# Boxplot Beta
order=c("enterotype__1", "enterotype__2")
myplot<-ggplot(ent_samples, aes(x=factor(Enterotype, levels = order), y=beta)) +
  geom_boxplot() +
  xlab("") +
  ylab("Beta diversity - PCoA1 (residuals)") +
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
  geom_signif(comparisons = list(c("enterotype__1", "enterotype__2")),
              map_signif_level = TRUE)
print(myplot)

# Make plot: Enterotype Beta
myplot<-ggplot(beta_df, aes(x=Axis.1, y=Axis.2, colour = as.factor(ent_samples$Enterotype))) +
  geom_point() + 
  #ggtitle(paste0("Subset: ",subset_id," / Cov: ",covariate, " / Dist. Type: ",beta)) +
  xlab("PCo1") +
  ylab("PCo2") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 17),
        axis.line = element_line(),
        plot.margin = margin(5.5,5.5, 5.5, 6.5, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        legend.text=element_text(size=13),
        legend.title = element_blank(),
        legend.position="bottom")
print(myplot)

# Make plot: Study
myplot<-ggplot(beta_df, aes(x=Axis.1, y=Axis.2, colour = as.factor(ent_samples$study))) +
  geom_point() + 
  #ggtitle(paste0("Subset: ",subset_id," / Cov: ",covariate, " / Dist. Type: ",beta)) +
  xlab("PCo1") +
  ylab("PCo2") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 17),
        axis.line = element_line(),
        plot.margin = margin(5.5,5.5, 5.5, 6.5, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        legend.text=element_text(size=13),
        legend.title = element_blank(),
        legend.position="bottom")
print(myplot)

