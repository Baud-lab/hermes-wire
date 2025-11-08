library(pheatmap, lib.loc = "/users/abaud/fmorillo/R/x86_64-pc-linux-gnu-library/4.2/")
#library(stats, lib.loc = "/users/abaud/fmorillo/R/x86_64-pc-linux-gnu-library/4.2/")
#library(compositions)
library(tidyverse)
library(WGCNA, lib.loc = "/users/abaud/fmorillo/R/x86_64-pc-linux-gnu-library/4.2/") 
library(dplyr, lib.loc = "/users/abaud/fmorillo/R/x86_64-pc-linux-gnu-library/4.2/")
library(ggplot2, lib.loc = "/users/abaud/fmorillo/R/x86_64-pc-linux-gnu-library/4.2/")
library(ggsignif, lib.loc = "/users/abaud/fmorillo/R/x86_64-pc-linux-gnu-library/4.2/")
library(patchwork, lib.loc = "/users/abaud/fmorillo/R/x86_64-pc-linux-gnu-library/4.2/")

# Read Enterotypes
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Tax_Mat_Net/Enterotypes/enterotypes.RData")
number=1
ent_samples=enterotypes_samples_objs[[number]]
stats_ent_samples=count(ent_samples,Enterotype,study,sex)


# Load distances
#print("Load taxonomic distances")
#load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Tax_Mat_Net/Diversity/residuals_qned_counts_beta.RData")
#taxon_dist=gp.dist_objs[[1]]

print("Load taxonomic residual values")

load("/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Tax_Mat_Net/Cluster_Analyses/residuals_qned_counts_clusters.RData")
clust_matrix=residuals_qned_counts_clusters_objs[[1]]
#clust_matrix_two=clust_matrix[rownames(clust_matrix)=="cluster__3" |
#                                rownames(clust_matrix)=="cluster__7", ]
#clust_matrix_four=clust_matrix[rownames(clust_matrix)=="cluster__3" |
#                                 rownames(clust_matrix)=="cluster__1" |
#                                 rownames(clust_matrix)=="cluster__7" |
#                                 rownames(clust_matrix)=="cluster__8",
#]

# Read Heritability
#heritability="/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Norm_nocluster_30k_full/Heritability/heritability.RData"
heritability="/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Tax_Mat_Net/Heritability/heritability.RData"
load(heritability)

# Load functional annotations
load("/users/abaud/data/secondary/indexes/gtdb/genes/genomes_comparison/eggnog_annotations/all_annotations/func_eggNOG_OGs.RData")
#load("/users/abaud/data/secondary/indexes/gtdb/genes/genomes_comparison/eggnog_annotations/all_annotations/func_CAZy.RData")
#load("/users/abaud/data/secondary/indexes/gtdb/genes/genomes_comparison/eggnog_annotations/all_annotations/func_PFAMs.RData")
load("/users/abaud/data/secondary/indexes/gtdb/genes/genomes_comparison/eggnog_annotations/all_annotations/func_matrices_v4.RData")
load("/users/abaud/data/secondary/indexes/gtdb/genes/genomes_comparison/eggnog_annotations/all_annotations/func_COG_pathway_category.RData")

load("/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Tax_Mat_Net/Cluster_Analyses/residuals_qned_counts_new.RData")
genus_matrix=residuals_qned_counts_objs[[5]]
genus_matrix_two=genus_matrix[rownames(genus_matrix)=="g__Bacteroides" | rownames(genus_matrix)=="g__Prevotella",]
rownames(genus_matrix_two)=gsub("g__","",rownames(genus_matrix_two))
#genus_matrix_seven=genus_matrix[rownames(genus_matrix)=="g__Bacteroides" | rownames(genus_matrix)=="g__Prevotella" |
#                                  rownames(genus_matrix)=="g__CAG-485" |
#                                  rownames(genus_matrix)=="g__CAG-873" |
#                                  rownames(genus_matrix)=="g__Acetatifactor" |
#                                  rownames(genus_matrix)=="g__CAG-95" |
#                                  rownames(genus_matrix)=="g__Helicobacter_D",
#]

print("Load functional residual values")
#load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Func_Mat_Net/Diversity/residuals_qned_counts_beta.RData")
#func_dist=gp.dist_objs[[4]]
load("/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Func_Mat_Net/Matrix_Processing/residuals_qned_counts.RData")
func_matrix=as.matrix(residuals_qned_counts_objs[[4]])

tax_matrices=c(
  "clust_matrix",
  #"clust_matrix_four"
  #"clust_matrix_two",
  "genus_matrix_two"
  #"genus_matrix_seven"
)

EC4_functions <- read_excel("paper_figures/EC4_functions.xlsx")
EC4_Guilds <- read_excel("paper_figures/EC4_Guilds.xlsx")
annotation=rbind(EC4_functions,EC4_Guilds)


total_bioprocesses=unique(c(annotation$Bioprocess))
colours <- metafolio::gg_color_hue(
  n = length(total_bioprocesses), 
  hue_min = 0,    # start at red
  hue_max = 330,  # nearly full circle
  c = 100,        # full chroma for vivid colors
  l = 60          # moderate luminance
)
shuffled_indices <- sample(length(colours))
colours <- colours[shuffled_indices]
names(colours) <- total_bioprocesses

# Define a named list as required by pheatmap
ann_colors <- list(
  Bioprocess = colours
  #Origin = setNames(c("coral", "lightgreen", "purple"), c("Correlation", "Enrichment - High heritability", "Enrichment - Low cohousing effects"))
)

## Spearman is common for ecological distances; 9 999 perms is standard
#mant <- mantel(taxon_dist, func_dist,
#               method      = "spearman",
#               permutations = 9999, 
#               na.rm = TRUE)      # handles any accidental NAs
#print(mant)


#shared <- intersect(colnames(clust_matrix), colnames(func_matrix))
## ── bind and correlate rows ──────────────────────────────────
##   cor() expects variables in columns, so we transpose first
#corr_full <- cor(t(rbind(clust_matrix, func_matrix)), method = "pearson")  # :contentReference[oaicite:3]{index=3}
#
## ── extract only the 2 × ~1000 block we need ─────────────────
#corr_sub  <- corr_full[ rownames(clust_matrix), rownames(func_matrix) ]

# Heuristics: ~0.30–0.40 in per column; ~0.30 in per row + padding

#cP <- corAndPvalue(t(clust_matrix), t(clust_matrix), method = "pearson")
#corr_plot <- cP$cor
#rownames(corr_plot)=gsub("cluster__","Guild ",rownames(corr_plot))
#colnames(corr_plot)=gsub("cluster__","Guild ",colnames(corr_plot))
## custom diverging palette (blue-white-red)
#col_fun <- colorRampPalette(c("#313695", "white", "#A50026"))
#
#myplot=pheatmap(
#  corr_plot,
#  color               = col_fun(100),
#  cluster_rows        = TRUE,        # only 2 rows; no need to cluster
#  cluster_cols        = TRUE,         # dendrogram on ECs helps grouping
#  clustering_distance_cols = "euclidean",
#  clustering_method   = "complete",
#  show_colnames       = TRUE,        # hide long EC labels; add later if needed
#  main                = "Pearson r correlations"
#)
#print(myplot)

for (mat in tax_matrices){
  
  print(mat)
  if(mat=="genus_matrix_two"){
    EC4_functions <- read_excel("paper_figures/EC4_functions.xlsx")
    w <- 9.5         # width in inches
    h <- 3.6
    num_heat=10
  } else {
    EC4_functions <- read_excel("paper_figures/EC4_Guilds.xlsx")
    w <- max(14, 0.35 * 33)          # width in inches
    h <- max(7,  4 + 0.30 * 12)
    num_heat=30
  }

  pdf(paste0("~/paper_figures/heatmap_taxa_and_function_",mat,".pdf"),
      #width = w, height = h,
      paper = "special", onefile = TRUE,
      useDingbats = FALSE, pointsize = 10)
  
  #EC4_functions=EC4_functions[EC4_functions$EC4 %in% ec_codes,]
  annotation=data.frame(
    Bioprocess = as.factor(EC4_functions$Bioprocess),
    row.names=EC4_functions$EC4)
  
  # ── calculate correlations & raw p-values in one call ───────
  cP <- corAndPvalue(t(get(mat)), t(func_matrix), method = "pearson")
  # ── Benjamini-Hochberg adjustment over *all* tests ───────────
  p_adj <- matrix(p.adjust(cP$p, method = "BH"), nrow = nrow(get(mat)),
                  dimnames = dimnames(cP$p))
  # ── filter EC columns significant for at least one cluster ───
  alpha <- 0.05
  keep  <- apply(p_adj, 2, function(x) any(x <= alpha))
  corr_sig <- cP$cor[, keep]        # 2 × n_sig_EC matrix
  
  ### 1.  For every function, find the largest |R²| across clusters
  #max_abs_r2 <- apply(corr_sig, 2, function(x) max(abs(x), na.rm = TRUE))
  ## names(max_abs_r2) are the same as row-names(mat)
  ### 2.  Order functions by that statistic (descending) …
  #top_idx <- order(max_abs_r2, decreasing = TRUE)[1:50]
  ### 3.  …and slice the matrix
  #corr_plot <- corr_sig[,top_idx, drop = FALSE]
  
  # ── filter EC columns based on biological relevance threshold ───
  thr <- 0.20
  keep_cols <- rep(TRUE, ncol(corr_sig))
  while (sum(keep_cols) > 50) {          # tighten until ≤ 80 columns remain
    keep_cols <- apply(corr_sig, 2, function(x) any(abs(x) >= thr))
    corr_sig  <- corr_sig[, keep_cols, drop = FALSE]
    thr <- thr + 0.05
  }
  corr_plot <- corr_sig
  #corr_plot <- cP$cor
  print(paste0("Final threshold: ",thr))
  print(paste0("Final functions: ",ncol(corr_plot)))
  keep_cols=colnames(corr_plot)
  #if (mat==tax_matrices[1] || mat==tax_matrices[2] || mat==tax_matrices[3]){
    rownames(corr_plot)=gsub("cluster__","Guild ",rownames(corr_plot)) 
  #}  
  
  if (nrow(corr_plot)<3){
    clust=F
  } else {clust=T}
  
  # custom diverging palette (blue-white-red)
  col_fun <- colorRampPalette(c("#313695", "white", "#A50026"))
  
  myplot=pheatmap(
    corr_plot,
    #treeheight_row = num_heat, treeheight_col = num_heat,
    color               = col_fun(100),
    cluster_rows        = clust,        # only 2 rows; no need to cluster
    cluster_cols        = TRUE,         # dendrogram on ECs helps grouping
    #clustering_distance_cols = "euclidean",
    clustering_method   = "ward.D2",
    show_colnames       = TRUE,        # hide long EC labels; add later if needed
    #annotation_col = annotation,
    #annotation_colors = ann_colors,
    #annotation_legend = TRUE,
    legend = TRUE,
    fontsize_row = 10,       # font size for row names
    fontsize_col = 10,       # font size for column names
    main                = "Pearson correlations (rho)"
  )
  print(myplot)
  
  write.table(corr_plot,file=paste0("/users/abaud/fmorillo/paper_figures/tables/ec_correlations_",mat,".csv"),sep=",",quote=F,col.names=T,row.names=T)
  
  # Define enterotypes
  enterotypes <- sort(as.character(unique(ent_samples$Enterotype)))
  
  # Split matrices by enterotype
  for (ent in enterotypes) {
    motch <- ent_samples$sample[ent_samples$Enterotype == ent]
    this_matrix <- get(mat)
    this_matrix <- this_matrix[, colnames(this_matrix) %in% motch]
    this_func_matrix <- func_matrix[, colnames(func_matrix) %in% motch]
    assign(paste(mat, "matrix", ent, sep = "_"), this_matrix)
    assign(paste("func_matrix", ent, sep = "_"), this_func_matrix)
  }
  
  # Calculate correlations for each enterotype
  mat_list <- list()
  for (ent in enterotypes) {
    ent_name <- if (ent == "enterotype__1") "Enterotype 1" else "Enterotype 2"
    message(ent_name)
    
    this_matrix <- get(paste(mat, "matrix", ent, sep = "_"))
    this_func_matrix <- get(paste("func_matrix", ent, sep = "_"))
    
    # Calculate correlations
    cP <- corAndPvalue(t(this_matrix), t(this_func_matrix), method = "pearson")
    corr_plot <- cP$cor
    corr_plot <- corr_plot[, colnames(corr_plot) %in% keep_cols]
    
    if (mat %in% tax_matrices[1:3]) {
      rownames(corr_plot) <- gsub("Guild ", "cluster__", rownames(corr_plot))
    }
    
    assign(paste("corr_plot", ent, sep = "_"), corr_plot)
    mat_list[[ent]] <- corr_plot
  }
  
  # Create a common color gradient across both enterotypes
  all_vals <- unlist(lapply(mat_list, function(m) as.numeric(m)))
  all_vals <- all_vals[is.finite(all_vals)]
  
  # Symmetric limits around zero for correlations
  if (mat!="genus_matrix_two"){
    vlim <- max(abs(all_vals), na.rm = TRUE)
    cols <- colorRampPalette(c("#313695", "white", "#A50026"))(100)
    brks <- seq(-vlim, vlim, length.out = length(cols) + 1)
    # Optional: nice legend ticks
    leg_brks <- seq(-vlim, vlim, length.out = 5)
    leg_labs <- round(leg_brks, 1)
  } else {
    vmin <- 0.2
    vmax <- 0.8
    center <- 0.5
    cols <- colorRampPalette(c("#313695", "white", "#A50026"))(100)
    
    # Create breaks with 0.6 as the "white" center
    n_colors <- 100
    n_lower <- round(n_colors * (center - vmin) / (vmax - vmin))
    n_upper <- n_colors - n_lower
    
    lower <- seq(vmin, center, length.out = n_lower)
    upper <- seq(center, vmax, length.out = n_upper + 1)
    brks <- c(lower, upper[-1])  # merge without duplicating center
  }
  
  mat_list=list()
  for (ent in enterotypes){
    
    if(ent=="enterotype__1"){
      ent_name="Enterotype 1"
    } else {ent_name="Enterotype 2"} #### Fix
    print(ent_name)
    this_matrix=get(paste(mat,"matrix",ent,sep="_"))
    this_func_matrix=get(paste("func_matrix",ent,sep="_"))
    
    # ── calculate correlations & raw p-values in one call ───────
    cP <- corAndPvalue(t(this_matrix), t(this_func_matrix), method = "pearson")
    corr_plot <- cP$cor
    corr_plot <- corr_plot[, colnames(corr_plot) %in% keep_cols]
    if (mat==tax_matrices[1] || mat==tax_matrices[2] || mat==tax_matrices[3]){
      rownames(corr_plot)=gsub("Guild ","cluster__",rownames(corr_plot)) 
    }
    assign(paste("corr_plot",ent,sep="_"),corr_plot)
    mat_list[[ent]] <- corr_plot
  }
  
  diff_plot=abs(corr_plot_enterotype__2)-abs(corr_plot_enterotype__1)
  diff_plot=t(diff_plot)
  diff_plot=as.data.frame(diff_plot)
  diff_plot$EC <- rownames(diff_plot)
  #diff_plot$Genus=rownames(diff_plot)
  
  # Reshape from wide to long
  diff_long <- diff_plot %>%
    pivot_longer(
      cols = c(colnames(diff_plot[-c(length(colnames(diff_plot)))])),
      names_to = "Cluster",
      values_to = "DeltaRho"
    )
  
  myplot<-ggplot(diff_long, aes(x=reorder(Cluster,DeltaRho,FUN = median), y=DeltaRho)) +
    geom_boxplot() +
    xlab("") +
    ylab("Ent2 – Ent1") +
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
    geom_signif(comparisons = list(c("Bacteroides","Prevotella")),
                map_signif_level = TRUE)
  print(myplot)
  
  ## ── calculate correlations & raw p-values in one call ───────
  #cP <- corAndPvalue(t(this_matrix), t(this_func_matrix), method = "pearson")
  ## ── Benjamini-Hochberg adjustment over *all* tests ───────────
  #p_adj <- matrix(p.adjust(cP$p, method = "BH"), nrow = nrow(this_matrix),
  #                dimnames = dimnames(cP$p))
  ## ── filter EC columns significant for at least one cluster ───
  #alpha <- 0.05
  #keep  <- apply(p_adj, 2, function(x) any(x <= alpha))
  #corr_sig <- cP$cor[, keep]        # 2 × n_sig_EC matrix
  #
  ## ── filter EC columns based on biological relevance threshold ───
  #thr <- 0.20
  #keep_cols <- rep(TRUE, ncol(corr_sig))
  #while (sum(keep_cols) > 80) {          # tighten until ≤ 80 columns remain
  #  keep_cols <- apply(corr_sig, 2, function(x) any(abs(x) >= thr))
  #  corr_sig  <- corr_sig[, keep_cols, drop = FALSE]
  #  thr <- thr + 0.05
  #}
  #corr_plot <- corr_sig 
  #print(paste0("Final threshold: ",thr))
  #print(paste0("Final functions: ",ncol(corr_plot)))
  
  for (ent in enterotypes){
    corr_plot <- get(paste("corr_plot", ent, sep = "_"))
    clust <- nrow(corr_plot) >= 3
    ent_name <- if (ent == "enterotype__1") "Enterotype 1" else "Enterotype 2"
    message(paste0("Final functions: ", ncol(corr_plot)))
    
    myplot=pheatmap(
      corr_plot,
      #treeheight_row = num_heat, treeheight_col = num_heat,
      color = cols,
      breaks = brks,   # << ensures identical gradient for both plots
      cluster_rows = clust,
      cluster_cols = TRUE,
      clustering_distance_cols = "euclidean",
      clustering_method = "complete",
      show_colnames = TRUE,
      main = paste0("Pearson correlations (rho): ", ent_name),
      legend_breaks = leg_brks,
      legend_labels = leg_labs,
      #annotation_col = annotation,
      #annotation_colors = ann_colors,
      #annotation_legend = TRUE,
      legend = TRUE,
      fontsize_row = 10,       # font size for row names
      fontsize_col = 10,       # font size for column names
    )
    print(myplot)
    
    write.table(corr_plot,file=paste0("/users/abaud/fmorillo/paper_figures/tables/ec_correlations_",mat,"_",ent,".csv"),sep=",",quote=F,col.names=T,row.names=T)
  }
  #
  #dE1 <- as.dist(1 - mat_list[[1]])
  #dE2 <- as.dist(1 - mat_list[[2]])
  #mant <- mantel(dE1, dE2, permutations = 9999)
  #print("Mantel test")
  #print(mant)          # reports Mantel r and permutation p-value
  #
  #perm  <- permustats(mant)                        # grab permutation stats  :contentReference[oaicite:1]{index=1}
  #obs.r <- mant$statistic                     # observed Mantel r
  #pval  <- mant$signif                       # permutation p-value
  #
  ### ── 2. Permutation density / histogram ──────────────────────────────────
  #g1 <- ggplot(data.frame(r = perm$perm), aes(x = r)) +
  #  geom_histogram(aes(y = after_stat(density)),
  #                 bins = 40, fill = "grey70", colour = "white") +
  #  geom_vline(xintercept = obs.r, colour = "red", linewidth = 1) +
  #  labs(x = "Mantel statistic (r)", y = "Density",
  #       title = "Permutation distribution of Mantel r") +
  #  theme_bw()
  #  #annotate("text", x = obs.r, y = Inf, vjust = 1.1,
  #  #         label = sprintf("Observed r = %.2f\np = %.4f", obs.r, pval),
  #  #         colour = "red", hjust = 0)
  #print(g1)
  #
  ### ── 3. Distance–distance scatter (upper-triangle only) ──────────────────
  ###     Convert both dist objects into the same long table
  #mat1  <- as.matrix(dE1)
  #mat2  <- as.matrix(dE2)
  #upper <- upper.tri(mat1)
  #df    <- data.frame(d1 = mat1[upper], d2 = mat2[upper])
  #
  #g2 <- ggplot(df, aes(d1, d2)) +
  #  geom_hex(bins = 40) +
  #  geom_smooth(method = "lm", colour = "black", se = FALSE) +
  #  labs(x = "Distance matrix 1", y = "Distance matrix 2",
  #       title = "Pairwise distances coloured by density") +
  #  theme_bw()
  #print(g2)
  #
  ### ── 4. Combine and display ──────────────────────────────────────────────
  #(g1 | g2) + plot_annotation(
  #  caption = "Mantel test: correlation between guild–guild distance matrices\n"  ,
  #  theme = theme(plot.caption = element_text(hjust = 0))
  #)
  dev.off()
}



#########################################################################

#mat1 <- read.csv("~/paper_figures/tables/ec_correlations_genus_matrix_seven_enterotype__1_full.csv", row.names = 1, check.names = FALSE)
#mat2 <- read.csv("~/paper_figures/tables/ec_correlations_genus_matrix_seven_enterotype__2_full.csv", row.names = 1, check.names = FALSE)
#
## Distance among guilds (rows).  Use (1 – Pearson r) so 0 = identical, 2 = opposite
#d_g1 <- as.dist(1 - cor(t(mat1), use = "pairwise.complete.obs"))
#d_g2 <- as.dist(1 - cor(t(mat2), use = "pairwise.complete.obs"))
#
## Distance among EC4 functions (columns)
#d_e1 <- as.dist(1 - cor(mat1,  use = "pairwise.complete.obs"))
#d_e2 <- as.dist(1 - cor(mat2,  use = "pairwise.complete.obs"))
#
## Global similarity: Mantel test
#
#library(vegan)
#mantel_guild  <- mantel(d_g1, d_g2, permutations = 999, method = "pearson")
#mantel_ec4    <- mantel(d_e1, d_e2, permutations = 999, method = "pearson")
#print(mantel_guild); print(mantel_ec4)
#
## Tree-shape similarity: dendrogram correlation
#
#library(dendextend)
#tree_g1 <- as.dendrogram(hclust(d_g1, method = "complete"))
#tree_g2 <- as.dendrogram(hclust(d_g2, method = "complete"))
#tree_e1 <- as.dendrogram(hclust(d_e1, method = "complete"))
#tree_e2 <- as.dendrogram(hclust(d_e2, method = "complete"))
#
## Baker's γ correlation and permutation test (10 000 perms)
#baker_guild <- cor_bakers_gamma(dendlist(tree_g1, tree_g2),
#                            #method = "baker",
#                            #type   = "permutation",
#                            nperm  = 10000)
#baker_ec4   <- cor.dendlist(dendlist(tree_e1, tree_e2),
#                            method = "baker",
#                            type   = "permutation",
#                            nperm  = 10000)
#print(baker_guild); print(baker_ec4)
#
## Partition Agreement
#library(mclust)
#k <- 4                            # pick k informed by your dendrogram
#cl_g1 <- cutree(tree_g1, k)
#cl_g2 <- cutree(tree_g2, k)
#ari_g <- adjustedRandIndex(cl_g1, cl_g2)
#
## Permutation test for ARI
#ari_null <- replicate(10000, adjustedRandIndex(sample(cl_g1), cl_g2))
#p_ari <- mean(ari_null >= ari_g)

## classical MDS of each distance matrix
#ord_g1 <- cmdscale(d_g1, k = 2)
#ord_g2 <- cmdscale(d_g2, k = 2)
#proc_g <- protest(ord_g1, ord_g2, permutations = 9999)
#print(proc_g)

## 1. Load packages
#suppressPackageStartupMessages({
#  library(dplyr)      # data wrangling
#  library(readr)      # fast read_csv
#  library(pheatmap)   # optional heatmaps
#})
#
## 2. Read correlation matrices
#mat_all <- read.csv("~/paper_figures/tables/ec_correlations_clust_matrix.csv", row.names = 1, check.names = FALSE)
#mat_all= as.matrix(mat_all)
#mat_e1  <- read.csv("~/paper_figures/tables/ec_correlations_clust_matrix_enterotype__1.csv", row.names = 1, check.names = FALSE)
#mat_e1= as.matrix(mat_e1)
#mat_e2  <- read.csv("~/paper_figures/tables/ec_correlations_clust_matrix_enterotype__2.csv", row.names = 1, check.names = FALSE)
#mat_e2= as.matrix(mat_e2)
#
## 3. Δρ matrix (Ent1 minus Ent2)
#delta   <- mat_e1 - mat_e2
#
## 4. Helper: identify top ↑ and ↓ EC4 per guild
#top_change <- function(v) {
#  c(up_EC4   = names(sort(v, decreasing = TRUE))[1],
#    up_val   = sort(v, decreasing = TRUE)[1],
#    down_EC4 = names(sort(v))[1],
#    down_val = sort(v)[1])
#}
#
#top_tbl <- t(apply(delta, 1, top_change)) %>%
#  as.data.frame() %>%
#  mutate(across(ends_with("val"), as.numeric),
#         guild = rownames(.)) %>%
#  dplyr::select(guild, everything())
#
## 5. Define three guild clusters
#clusters <- list(
#  mucus        = c("Guild 0","Guild 1","Guild 3","Guild 5","Guild 6"),
#  lumen        = c("Guild 7","Guild 8","Guild 10","Guild 11"),
#  intermediate = c("Guild 2","Guild 4","Guild 9")
#)
#
## 6. Compute average absolute change per cluster
#cluster_stats <- lapply(clusters, function(g) {
#  sub <- delta[g, , drop = FALSE]
#  data.frame(mean_abs = mean(abs(sub)),
#             mean_pos = mean(sub[sub >  0]),
#             mean_neg = mean(sub[sub <  0]))
#}) %>% bind_rows(.id = "cluster")
#
# #7. (Optional) Visualise Δρ heat-map
# pheatmap(delta,
#          clustering_distance_rows = "euclidean",
#          clustering_distance_cols = "euclidean",
#          color = colorRampPalette(c("navy","white","firebrick"))(100),
#          main = expression(paste(Delta, rho," (Ent1 - Ent2)")))
#
## 8. Save results
#write.csv(top_tbl, "guild_top_changes_Ent1_vs_Ent2.csv", row.names = FALSE)
#write.csv(cluster_stats, "cluster_summary_Ent1_vs_Ent2.csv", row.names = FALSE)
