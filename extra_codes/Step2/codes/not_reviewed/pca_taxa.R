load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Shallow/Matrix_Processing/residuals_qned_counts.RData")
matrix=t(residuals_qned_counts_objs[[3]])

pca_matrix <- prcomp(matrix, scale=F, center=F)
summary(pca_matrix)
pca_scores <-pca_matrix$x

#pca_scores[,1:2]

# Plot the first two principal components
plot(pca_scores[, 1], pca_scores[, 2],
     xlab = "Principal Component 1",
     ylab = "Principal Component 2",
     #main = "PCA of UBA4372 KEGG KO",
     pch = 19, col = "lightgrey")
#text(pca_scores[, 1], pca_scores[, 2], labels = rownames(pca_scores), pos = 4, cex = 0.7, col = "black")

# Extract loadings (eigenvectors)
pca_loadings <- pca_matrix$rotation
loadings_pc <- pca_loadings[, 1]
sorted_loadings_pc <- sort(abs(loadings_pc), decreasing = TRUE)
top_contributors <- names(sorted_loadings_pc)[1:5]  # Adjust the number to see more/less top contributors
top_contributions <- loadings_pc[top_contributors]
print(top_contributions)
