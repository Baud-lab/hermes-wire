Rscript /users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/bin/matrix_processing_network.R \
-met "Shallow" \
-phe "microbiome_taxonomic" \
-mod_comp "NO" \
-residues "/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Norm_nocluster_30k_full/Cluster_Analyses/residuals_qned_counts_new.RData" \
-ent "/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Norm_nocluster_30k_full/Enterotypes/enterotypes.RData" \ 
-covariates "/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/covariates_Shallow.txt" \
-cov_types "/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/covariates_classification.txt" \
-tax "/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/taxonomy_GTDB207.txt" \
-tax_ranks "Phylum,Class,Order,Family,Genus,Species" \
-func_ranks "EC1,EC2,EC3,EC4" \
-measure "bicor" \
-cluster "cluster_walktrap" \
-hub "degree,closeness,betweenness,eigenvector" \
-umap "/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Norm_nocluster_30k_full/Cluster_Analyses/umap_dir" \


met="Shallow"
phe="microbiome_taxonomic"
mod_comp="NO"
residues="/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Norm_nocluster_30k_full/Cluster_Analyses/residuals_qned_counts_new.RData" 
ent="/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Norm_nocluster_30k_full/Enterotypes/enterotypes.RData"
covariates="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/covariates_Shallow.txt"
cov_types="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/covariates_classification.txt"
tax="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/taxonomy_GTDB207.txt"
tax_ranks="Phylum,Class,Order,Family,Genus,Species"
func_ranks="EC1,EC2,EC3,EC4"
measure="bicor"
cluster="cluster_walktrap"
hub="degree,closeness,betweenness,eigenvector"
umap="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Norm_nocluster_30k_full/Cluster_Analyses/umap_dir"
umap_res="/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Norm_nocluster_30k_full/Matrix_Processing/matrices_residues/"
Rscript /users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/bin/matrix_processing_network.R -met $met -phe $phe -mod_comp $mod_comp -residues $residues -ent $ent -covariates $covariates -cov_types $cov_types -umap $umap -umap_res $umap_res -tax $tax -tax_ranks $tax_ranks -func_ranks $func_ranks -measure $measure -cluster $cluster -hub $hub



##########


# 1.1. Jaccard indices for centralities
jacc_deg  <- cmp$jaccDeg   # Named vector: value "jacc", "p.greater", "p.less"
jacc_betw <- cmp$jaccBetw
jacc_close<- cmp$jaccClose
jacc_eigen<- cmp$jaccEigen
jacc_hub  <- cmp$jaccHub   # overlap of hub sets by chosen hub quantile

# 1.2. Adjusted Rand index for clustering similarity
rand_idx      <- cmp$randInd       # Adjusted Rand index for full networks
rand_idx_lcc  <- cmp$randIndLCC    # Adjusted Rand index for largest connected component

# 1.3. Print results in a tidy format
jacc_df <- tibble(
  centrality = c("degree","betweenness","closeness","eigenvector","hubs"),
  jacc       = c(jacc_deg["jacc"], jacc_betw["jacc"], jacc_close["jacc"], jacc_eigen["jacc"], jacc_hub["jacc"]),
  p_greater  = c(jacc_deg["p.greater"], jacc_betw["p.greater"], jacc_close["p.greater"], jacc_eigen["p.greater"], jacc_hub["p.greater"]),
  p_less     = c(jacc_deg["p.less"], jacc_betw["p.less"], jacc_close["p.less"], jacc_eigen["p.less"], jacc_hub["p.less"])
)
print("Jaccard indices for centralities and hubs:")
print(jacc_df)

rand_df <- tibble(
  metric = c("Rand index (full)", "Rand index (LCC)"),
  value  = c(rand_idx["value"], rand_idx_lcc["value"]),
  pval   = c(rand_idx["pval"], rand_idx_lcc["pval"])
)
print("Adjusted Rand indices for cluster similarity:")
print(rand_df)

# 2.1. Whole-network GCD
gcd_obj <- cmp$gcd         # This is an object (list) with elements $gcd, $ocount1, $ocount2, $gcm1, $gcm2
gcd_value <- gcd_obj$gcd
print(paste("Graphlet Correlation Distance (full network):", round(gcd_value, 3)))

# 2.2. LCC GCD
gcd_lcc_obj <- cmp$gcdLCC
gcd_lcc_value <- gcd_lcc_obj$gcd
print(paste("Graphlet Correlation Distance (LCC):", round(gcd_lcc_value, 3)))


# 3.1. Absolute global properties for each network
props_list <- cmp$properties
# Example names in props_list: deg1, deg2, betw1, betw2, avDiss1, avDiss2, avPath1, avPath2, density1, density2, clustCoef1, clustCoef2, modularity1, modularity2, etc.
# Extract a selection of key metrics:
global_props_df <- tibble(
  metric = c("average dissimilarity", "average path length", "density", "clustering coefficient", "modularity"),
  network1 = c(props_list$avDiss1, props_list$avPath1, props_list$density1, props_list$clustCoef1, props_list$modularity1),
  network2 = c(props_list$avDiss2, props_list$avPath2, props_list$density2, props_list$clustCoef2, props_list$modularity2)
)
print("Global network properties for each network:")
print(global_props_df)

# 3.2. Differences
diffG <- cmp$diffGlobal
diff_global_df <- tibble(
  metric = c("nComp","avgDiss","avgPath","density","vertConnect","edgeConnect","natConnect","PEP","clustCoef","modularity"),
  diff   = c(diffG$diffnComp, diffG$diffavDiss, diffG$diffavPath, diffG$diffDensity, diffG$diffVertConnect, diffG$diffEdgeConnect, diffG$diffNatConnect, diffG$diffPEP, diffG$diffClustCoef, diffG$diffModul)
)
print("Differences (Network2 minus Network1) in global properties:")
print(diff_global_df)

# 3.3. Largest connected component (LCC) properties
propsLCC <- cmp$propertiesLCC
lcc_props_df <- tibble(
  metric = c("LCC size (abs)", "LCC size (rel)", "avgDiss","avgPath","density","vertConnect","edgeConnect","natConnect","PEP","clustCoef","modularity"),
  network1 = c(propsLCC$lccSize1, propsLCC$lccSizeRel1, propsLCC$avDiss1, propsLCC$avPath1, propsLCC$density1, propsLCC$vertConnect1, propsLCC$edgeConnect1, propsLCC$natConnect1, propsLCC$pep1, propsLCC$clustCoef1, propsLCC$modularity1),
  network2 = c(propsLCC$lccSize2, propsLCC$lccSizeRel2, propsLCC$avDiss2, propsLCC$avPath2, propsLCC$density2, propsLCC$vertConnect2, propsLCC$edgeConnect2, propsLCC$natConnect2, propsLCC$pep2, propsLCC$clustCoef2, propsLCC$modularity2)
)
print("LCC properties for each network:")
print(lcc_props_df)

diffGLCC <- cmp$diffGlobalLCC
diff_lcc_df <- tibble(
  metric = c("lccSize","lccSizeRel","avgDiss","avgPath","density","vertConnect","edgeConnect","natConnect","PEP","clustCoef","modularity"),
  diff   = c(diffGLCC$difflccSize, diffGLCC$difflccSizeRel, diffGLCC$diffavDiss, diffGLCC$diffavPath, diffGLCC$diffDensity, diffGLCC$diffVertConnect, diffGLCC$diffEdgeConnect, diffGLCC$diffNatConnect, diffGLCC$diffPEP, diffGLCC$diffClustCoef, diffGLCC$diffModul)
)
print("Differences (Network2 minus Network1) in LCC properties:")
print(diff_lcc_df)


# 4.1. Convert diffCentr lists to data frames
diff_deg  <- data.frame(taxon = names(cmp$diffCentr$diffDeg),  diff = cmp$diffCentr$diffDeg,  stringsAsFactors = FALSE)
diff_betw <- data.frame(taxon = names(cmp$diffCentr$diffBetw), diff = cmp$diffCentr$diffBetw, stringsAsFactors = FALSE)
diff_close<- data.frame(taxon = names(cmp$diffCentr$diffClose),diff = cmp$diffCentr$diffClose,stringsAsFactors = FALSE)
diff_eig  <- data.frame(taxon = names(cmp$diffCentr$diffEigen),diff = cmp$diffCentr$diffEigen,stringsAsFactors = FALSE)

# 4.2. Rank taxa by absolute change in centrality
top_n <- 10  # number of top taxa to display
top_deg  <- diff_deg  %>% mutate(absdiff = abs(diff)) %>% arrange(desc(absdiff)) %>% slice_head(n = top_n)
top_betw <- diff_betw %>% mutate(absdiff = abs(diff)) %>% arrange(desc(absdiff)) %>% slice_head(n = top_n)
top_close<- diff_close%>% mutate(absdiff = abs(diff)) %>% arrange(desc(absdiff)) %>% slice_head(n = top_n)
top_eig  <- diff_eig  %>% mutate(absdiff = abs(diff)) %>% arrange(desc(absdiff)) %>% slice_head(n = top_n)

print("Top taxa by absolute change in degree centrality (network2 minus network1):")
print(top_deg)
print("Top taxa by absolute change in betweenness centrality:")
print(top_betw)
print("Top taxa by absolute change in closeness centrality:")
print(top_close)
print("Top taxa by absolute change in eigenvector centrality:")
print(top_eig)



# 5.1. Check if diffEdges exists
if (!is.null(cmp$diffEdges)) {
  diff_edges <- cmp$diffEdges
  # Preview top differential edges by absolute difference in association strength
  # Assuming diffEdges has columns: taxon1, taxon2, asso1, asso2, diff, p.value (check actual column names)
  # Adapt column names based on actual structure:
  # E.g., names(diff_edges) might include "taxa1", "taxa2", "as1", "as2", "diff", "p.adjust"
  print("Head of diffEdges:")
  print(head(diff_edges))
  
  # If columns known, rank by absolute difference:
  if (all(c("as1","as2") %in% colnames(diff_edges))) {
    diff_edges <- diff_edges %>%
      mutate(absdiff = abs(as2 - as1)) %>%
      arrange(desc(absdiff))
    print(paste("Top", top_n, "differential edges by absolute association change:"))
    print(head(diff_edges, n = top_n))
  }
} else {
  message("No diffEdges element found in cmp; consider computing differential network via diffnet().")
}



# Example for degree centrality:
deg1 <- cmp$properties$deg1
deg2 <- cmp$properties$deg2
df_deg <- tibble(
  taxon = names(deg1),
  degree = deg1,
  network = "Network1"
) %>%
  bind_rows(
    tibble(taxon = names(deg2), degree = deg2, network = "Network2")
  )
ggplot(df_deg, aes(x = degree, fill = network)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  theme_minimal() +
  labs(title = "Degree centrality distributions", x = "Degree", y = "Count")


# Example for betweeness centrality:
betw1 <- cmp$properties$betw1
betw2 <- cmp$properties$betw2
df_deg <- tibble(
  taxon = names(betw1),
  betweeness = betw1,
  network = "Network1"
) %>%
  bind_rows(
    tibble(taxon = names(betw2), betweeness = betw2, network = "Network2")
  )
ggplot(df_deg, aes(x = betweeness, fill = network)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  theme_minimal() +
  labs(title = "Betweeness centrality distributions", x = "Betweeness", y = "Count")


# Example for closeness centrality:
close1 <- cmp$properties$close1
close2 <- cmp$properties$close2
df_close <- tibble(
  taxon = names(close1),
  close = close1,
  network = "Network1"
) %>%
  bind_rows(
    tibble(taxon = names(close2), close = close2, network = "Network2")
  )
ggplot(df_close, aes(x = close, fill = network)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  theme_minimal() +
  labs(title = "Closeness centrality distributions", x = "Closeness", y = "Count")

# Example for closeness centrality:
eigen1 <- cmp$properties$eigen1
eigen2 <- cmp$properties$eigen2
df_eigen <- tibble(
  taxon = names(eigen1),
  eigen = eigen1,
  network = "Network1"
) %>%
  bind_rows(
    tibble(taxon = names(eigen2), eigen = eigen2, network = "Network2")
  )
ggplot(df_eigen, aes(x = eigen, fill = network)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  theme_minimal() +
  labs(title = "Eigenvalue centrality distributions", x = "Eigenvalue", y = "Count")



cl1 <- cmp$properties$clust1
cl2 <- cmp$properties$clust2
# Build a contingency table
ct <- table(clust1 = cl1, clust2 = cl2)
print("Contingency table of module memberships between networks:")
print(ct)
# Compute adjusted Rand index separately if desired:
library(mclust)  # for adjustedRandIndex
ari_val <- adjustedRandIndex(cl1, cl2)
print(paste("Adjusted Rand Index (mclust) between clusterings:", round(ari_val, 3)))
