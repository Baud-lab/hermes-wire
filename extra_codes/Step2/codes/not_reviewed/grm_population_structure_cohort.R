## ---- packages ----
suppressPackageStartupMessages({
  library(RMTstat)
  library(ggplot2)
  library(vegan)
})

################

# Arguments

subset_id="ALL"
# "ALL" "MI" "NY"
method="Shallow"

set.seed(123)

###############

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

check_grm_centered <- function(GRM) {
  n  <- nrow(grm)
  J  <- diag(n) - matrix(1/n, n, n)        # centring matrix
  tol <- 1e-10                              # numerical tolerance
  
  rowZero  <- max(abs(rowMeans(grm)))  < tol
  colZero  <- max(abs(colMeans(grm)))  < tol
  nullOne  <- abs(drop(rep(1, n) %*% grm %*% rep(1, n))) < tol
  idem     <- max(abs(grm - J %*% grm %*% J)) < tol
  offDiag  <- abs(mean(grm[row(grm) != col(grm)])) < tol
  
  centered=rowZero && colZero && nullOne && idem && offDiag
  return(centered)
}

tw_from_grm <- function(grm, centered, alpha = 0.05, drop = 3,
                        adjust = c("bonferroni", "fdr")) {
  
  ## --- 1. eigen-decomposition of centred GRM ------------------------------
  
  if (centered){
    grm_cent=grm
  } else{
    grm_cent <- scale(grm, center = TRUE, scale = FALSE)
  }
  
  eig      <- sort(eigen(grm_cent, symmetric = TRUE,
                         only.values = TRUE)$values, decreasing = TRUE)
  N        <- length(eig)
  
  ## --- 2. estimate M by fitting MP bulk -----------------------------------
  lambda_plus_hat <- max(eig[-seq_len(drop)])        # skip top 'drop' PCs
  c_hat           <- (sqrt(lambda_plus_hat) - 1)^2
  M_hat           <- round(N / c_hat)
  
  ## --- 3. Tracy-Widom test with RMTstat -----------------------------------
  library(RMTstat)
  pars   <- WishartMaxPar(ndf = M_hat, pdim = N, var = 1, beta = 1)
  z      <- (eig - pars$centering) / pars$scaling
  p_raw  <- ptw(z, beta = 1, lower.tail = FALSE)
  
  adjust <- match.arg(adjust)
  p_adj  <- p.adjust(p_raw, method = ifelse(adjust == "bonferroni",
                                            "bonferroni", "BH"))
  k      <- sum(p_adj < alpha)
  
  list(k         = k,
       p.values  = p_adj,
       z         = z,
       M.hat     = M_hat,
       lambda.plus  = lambda_plus_hat)
}


########################

# Load

metadata <- read.delim("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/metadata_Shallow.txt")

residuals="~/NextFlow/Microbiome_Profiling/step2/Output_Norm_nocluster_30k_full/Matrix_Processing/residuals_qned_counts.RData"
load(residuals)
residuals_ALL=residuals_qned_counts_objs[[3]]
colnames(residuals_ALL)=sapply(strsplit(colnames(residuals_ALL), "_"),"[",1)
samples=unique(c(colnames(residuals_ALL)))

# GRM
print("Load GRM")
load("/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/grm.RData")
grm=grm[intersect(rownames(grm), samples),intersect(rownames(grm), samples)]

this_metadata=metadata[metadata$RFID %in% rownames(grm),]
rownames(this_metadata)=this_metadata$RFID

if (!identical(rownames(grm), rownames(this_metadata))) {
  this_metadata <- this_metadata[rownames(grm), , drop = FALSE]
}

################

# Test population structure

centered=check_grm_centered(grm)
cat(paste0("Centered: ",centered,"\n\n"))

out <- tw_from_grm(grm, centered, alpha = 0.05, drop = 3)
k=out$k
cat("Significant PCs =", out$k, "\nEffective M  =", out$M.hat, "\n")

# Eigen decomposition of GRM
if (centered){
  grm_cent=grm
} else{
  grm_cent <- scale(grm, center = TRUE, scale = FALSE)
}

eig <- eigen(grm_cent, symmetric = TRUE)
eigenvalues <- eig$values
eigenvectors <- eig$vectors

# Prepare MANOVA data
df_manova <- data.frame(
  eigenvectors[, 1:out$k, drop = FALSE], 
  Study = this_metadata$Study
)
colnames(df_manova)[1:k] <- paste0("PC", 1:k)

# MANOVA with Pillai's test
manova_res <- summary.manova(
  manova(as.matrix(df_manova[, 1:k]) ~ Study, data = df_manova),
  test = "Pillai"
)

# Results
variance_explained <- 100 * sum(eigenvalues[1:k]) / sum(eigenvalues)
p_value <- manova_res$stats["Study", "Pr(>F)"]

cat(sprintf(
  "Top %d PCs explain %.1f%% variance\nMANOVA p-value = %.3e",
  k, variance_explained, p_value
))

# 4. Permutation MANOVA (999 permutations)

adonis_result <- adonis2(df_manova[, 1:k] ~ Study, 
                         data = df_manova, 
                         method = "euclidean", 
                         permutations = 999)

# 5. Extract and report results
p_value <- adonis_result[1, "Pr(>F)"]
p_lab <- bquote(italic(p) == .(formatC(p_value, digits = 3)))
r_squared <- adonis_result[1, "R2"]

cat("PERMANOVA Results:\n",
    "F =", round(adonis_result[1, "F"], 2), 
    "R² =", round(r_squared, 3),
    "p =", p_value, "\n")

save(manova_res,
     adonis_result,
     file="/users/abaud/fmorillo/paper_figures/cohort_population_structure.RData")

pdf("/users/abaud/fmorillo/paper_figures/cohort_population_structure.pdf")
myplot<-ggplot(df_manova, aes(x=PC1, y=PC2, colour = as.factor(Study))) +
  geom_point() + 
  ggtitle("Principal components of genetic relatedness") +
  scale_colour_manual(values = c("MI" = "#1f77b4", "NY" = "#ff7f0e")) +
  xlab("PC1") +
  ylab("PC2") +
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
        legend.position="bottom") +
  annotate("label", 
           x = max(df_manova$PC1), 
           y = max(df_manova$PC2), 
           label = as.expression(p_lab),   # supply the expression itself
           parse = TRUE,                 # tell ggplot to interpret it
           hjust = 1, vjust = 1,
           size = 5,
           label.size = 0.5,
           label.r = unit(0.2, "lines"),
           label.padding = unit(0.25, "lines"),
           fill = "white",
           color = "black"
           )
print(myplot)
dev.off()
