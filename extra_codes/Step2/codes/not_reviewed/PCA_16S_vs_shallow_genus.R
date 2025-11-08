#load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/old_Outputs/old_comparisons_16S_SS/Output_16S_Sample_Harm/Harmonization/harmonized_samples_and_taxa.RData")
#genus_16S=filtered_prev_objs[[15]]
#genus_shallow=filtered_prev_objs[[36]]

load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Shallow_Full_Harm/Matrix_Processing/residuals_qned_counts.RData")
residuals_SS=residuals_qned_counts_objs
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_16S_Full_Harm/Matrix_Processing/residuals_qned_counts.RData")
residuals_16S=residuals_qned_counts_objs

ranks=c("Phylum","Class","Order","Family","Genus","Species")
names(ranks)=c(3,6,9,12,15,18)
pdf('~/paper_figures/pca_16s_shallow.pdf')
for (rank in ranks) {
  number=names(ranks)[which(ranks==rank)]
  resids_ss=residuals_SS[[as.numeric(number)]]
  resids_16S=residuals_16S[[as.numeric(number)]]
  # Proportion of variance explained
  pca1=prcomp(t(resids_ss), scale=F, center=F)
  pca2=prcomp(t(resids_16S), scale=F, center=F)
  prop_var1 <- pca1$sdev^2 / sum(pca1$sdev^2)
  prop_var2 <- pca2$sdev^2 / sum(pca2$sdev^2)
  # Plotting
  #max=max(max(prop_var1),max(prop_var2))
  plot(prop_var1, type = "b", col = "blue", xlab = "Principal Component", ylab = "Proportion of Variance", ylim=c(0, 0.40), main = paste0("Variance Explained by PCs (Rank: ",rank,")"))
  points(prop_var2, type = "b", col = "red")
  legend("topright", legend = c("Shallow Shotgun", "16S"), col = c("blue", "red"), pch = 1)
}
dev.off()






#######

s_centered_pca <- prcomp(t(genus_shallow), scale=F, center=F)
g16s_centered_pca <- prcomp(t(genus_16S), scale=F, center=F)

summary_shallow=as.data.frame(t(summary(s_centered_pca)$importance))
summary_shallow$method="Shallow Shotgun"
summary_shallow$PC=rownames(summary_shallow)
summary_shallow <- summary_shallow %>%
  mutate(PC = factor(PC, levels = paste0("PC", sort(as.numeric(gsub("PC", "", PC))))))
summary_16S=as.data.frame(t(summary(g16s_centered_pca)$importance))
summary_16S$method="16S"
summary_16S$PC=rownames(summary_16S)
summary_16S <- summary_16S %>%
  mutate(PC = factor(PC, levels = paste0("PC", sort(as.numeric(gsub("PC", "", PC))))))
plot=rbind(summary_shallow,summary_16S)
plot=plot[plot$PC=="PC1" | plot$PC=="PC2" | plot$PC=="PC3" | plot$PC=="PC4" | plot$PC=="PC5",]

ggplot(plot, aes(x = PC, y = `Cumulative Proportion`, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = NULL,
       x = "",
       y = "Cummulative proportion of variance explained",
       fill = "Method") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.line = element_line(),legend.title = element_text(size = 10),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.0)) +
  scale_y_continuous(breaks = seq(0, max(plot$`Cumulative Proportion`, na.rm = TRUE), by = 0.1))

########

rcond(t(genus_shallow))
rcond(t(genus_16S))

#########

# Determine the number of principal components to retain
# Use a scree plot or cumulative variance to decide
# For example, retain components that explain up to 90% of the variance
num_components_1 <- which(cumsum(s_centered_pca$sdev^2 / sum(s_centered_pca$sdev^2)) > 0.90)[1]
num_components_2 <- which(cumsum(g16s_centered_pca$sdev^2 / sum(g16s_centered_pca$sdev^2)) > 0.90)[1]

# Reconstruct the data using the selected principal components
reconstructed_method1 <- s_centered_pca$x[, 1:num_components_1] %*% t(s_centered_pca$rotation[, 1:num_components_1])
reconstructed_method2 <- g16s_centered_pca$x[, 1:num_components_2] %*% t(g16s_centered_pca$rotation[, 1:num_components_2])

# Calculate the residues (differences between original data and reconstructed data)
residues_method1 <- t(genus_shallow) - reconstructed_method1
residues_method2 <- t(genus_16S) - reconstructed_method2

# Statistical analysis of residues

# 1. Calculate the variance of the residues for each method
#variance_residues_method1 <- apply(residues_method1, 2, var)
#variance_residues_method2 <- apply(residues_method2, 2, var)

# 2. Compare the means of residues
mean_residues_method1 <- mean(residues_method1)
mean_residues_method2 <- mean(residues_method2)

# Print the means
cat("Mean of residues for Method 1:", mean_residues_method1, "\n")
cat("Mean of residues for Method 2:", mean_residues_method2, "\n")

# 3. Perform a statistical test to compare the residues
# Using paired t-test to compare the means of variances
t_test_result <- t.test(residues_method1, residues_method2, paired = TRUE)

# Print t-test results
print(t_test_result)

# 4. Visualize the residues' variance distribution
residues_df <- data.frame(
  Method = rep(c("Method1", "Method2"), each = ncol(data)),
  Variance = c(variance_residues_method1, variance_residues_method2)
)

ggplot(residues_df, aes(x = Method, y = Variance)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  theme_minimal() +
  labs(title = "Variance of Residues for Each Method", y = "Variance")


####

# Proportion of variance explained

pca1=s_centered_pca
pca2=g16s_centered_pca

prop_var1 <- pca1$sdev^2 / sum(pca1$sdev^2)
prop_var2 <- pca2$sdev^2 / sum(pca2$sdev^2)

# Plotting
pdf('~/paper_figures/pca_16s_shallow.pdf')
plot(prop_var1, type = "b", col = "blue", xlab = "Principal Component", ylab = "Proportion of Variance", main = "Variance Explained by PCs")
points(prop_var2, type = "b", col = "red")
legend("topright", legend = c("Shallow Shotgun", "16S"), col = c("blue", "red"), pch = 1)
dev.off()

# Bartlett's test
bartlett_test <- bartlett.test(list(pca1$sdev, pca2$sdev))
print(bartlett_test)

# MANOVA
pcs1 <- as.data.frame(pca1$x)
pcs2 <- as.data.frame(pca2$x)

# Combine datasets for MANOVA
combined_pcs <- rbind(pcs1, pcs2)
group <- factor(c(rep(1, nrow(pcs1)), rep(2, nrow(pcs2))))

manova_test <- manova(cbind(PC1, PC2, PC3) ~ group, data = combined_pcs)
summary(manova_test, test = "Pillai")
