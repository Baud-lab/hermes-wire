# Load required libraries
library(rhdf5, lib.loc = "/users/abaud/fmorillo/R/x86_64-pc-linux-gnu-library/4.2/")
library(sommer)  # For mixed models with GRM
library(tidyverse, lib.loc = "/users/abaud/fmorillo/R/x86_64-pc-linux-gnu-library/4.2/")  # For data manipulation

# Load relative abundances (CLR-transformed residuals)
print("Load relative abundances (CLR-transformed residuals)")
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Norm_nocluster_30k_full/Cluster_Analyses/residuals_qned_counts_new.RData")
resids_NY=residuals_qned_counts_objs[[14]]
resids_MI=residuals_qned_counts_objs[[13]]
rownames(resids_MI)<- gsub("[-+.]", "_", rownames(resids_MI))
rownames(resids_NY)<- gsub("[-+.]", "_", rownames(resids_NY))
colnames(resids_MI)=sapply(strsplit(colnames(resids_MI), "_"), "[",1)
colnames(resids_NY)=sapply(strsplit(colnames(resids_NY), "_"), "[",1)
samples=sapply(strsplit(unique(c(colnames(resids_MI),colnames(resids_NY))), "_"), "[",1)

# Match taxa between cohorts (assuming rows are taxa, columns are samples)
shared_taxa <- intersect(rownames(resids_MI), rownames(resids_NY))
if(length(shared_taxa) == 0) {
  stop("No common taxa found between cohorts")
}

# Create metadata for modeling
print("Load metadata")
metadata <- read.delim("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/metadata_Shallow.txt")
metadata <- metadata %>% filter(RFID %in% samples)
meta_MI <- metadata[metadata$Study == "MI", ]
meta_NY <- metadata[metadata$Study == "NY", ]
meta_MI$Sample <- paste0("MI_", meta_MI$RFID)
meta_NY$Sample <- paste0("NY_", meta_NY$RFID)
metadata_combined <- rbind(meta_MI, meta_NY)

# Define maternal and cohousing effects
print("Define maternal and cohousing effects")
df=data.frame(RFID=metadata_combined$RFID,
              maternal=as.factor(metadata_combined$dam),
              cohousing=as.factor(metadata_combined$cage))
df=na.omit(df)
samples=df$RFID
metadata_combined <- metadata_combined %>% filter(RFID %in% samples)

# Ensure GRM row/colnames match RFID in metadata
print("Load GRM")
h5="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/P50_rats_Rn7_cp2.h5"
h5f = H5Fopen(h5)
grm=h5read(h5f, paste0("/GRM/P50_Rn7_pruned/","/matrix"))
ids=h5read(h5f, paste0("/GRM/P50_Rn7_pruned/row_header/","/sample_ID"))
H5Fclose(h5f)
rownames(grm)=ids
colnames(grm)=ids
grm=grm[intersect(rownames(grm), samples),intersect(rownames(grm), samples)]
metadata_combined <- metadata_combined %>% filter(RFID %in% rownames(grm))
rownames(grm) <- colnames(grm)  <- metadata_combined$RFID
new_ids <- metadata_combined$Sample
grm_expanded <- matrix(0, nrow = length(new_ids), ncol = length(new_ids))
rownames(grm_expanded) <- colnames(grm_expanded) <- new_ids
print("Expand GRM according to the cohorts")
for (i in seq_along(new_ids)) {
  for (j in seq_along(new_ids)) {
    id1 <- sub("^(MI|NY)_", "", new_ids[i])
    id2 <- sub("^(MI|NY)_", "", new_ids[j])
    grm_expanded[i, j] <- grm[id1, id2]
  }
}

# Ensure residuals are aligned with metadata
print("Ensure residuals are aligned with metadata")
resids_MI <- resids_MI[,colnames(resids_MI) %in% metadata_combined$RFID[metadata_combined$Study == "MI"]]  # Rows = taxa, columns = samples
resids_NY <- resids_NY[,colnames(resids_NY) %in% metadata_combined$RFID[metadata_combined$Study == "NY"]]
rownames(resids_MI)<- gsub("[-+.]", "_", rownames(resids_MI))
rownames(resids_NY)<- gsub("[-+.]", "_", rownames(resids_NY))
resids_MI=t(resids_MI)
resids_NY=t(resids_NY)

# Prepare results storage
genetic_cor_results <- data.frame(
  Taxon = character(),
  Genetic_Correlation = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each shared taxon
print("Loop through each shared taxon")
for (taxon in shared_taxa) {
  print(taxon)
  
  # Extract vectors
  print("Extract vectors")
  y_MI <- resids_MI[, taxon, drop = FALSE]
  y_NY <- resids_NY[, taxon, drop = FALSE]
  
  # Add cohort identifiers to sample names
  print("Add cohort identifiers to sample names")
  rownames(y_MI) <- paste0("MI_", rownames(resids_MI))
  rownames(y_NY) <- paste0("NY_", rownames(resids_NY))
  
  # Combine responses
  print("Combine responses")
  y_combined <- rbind(y_MI, y_NY)
  colnames(y_combined) <- "y"
  
  # Add trait info (MI or NY version of the taxon)
  print("Add trait info")
  metadata_combined$Trait <- ifelse(metadata_combined$Study == "MI",
                                    paste0(taxon, "_MI"), paste0(taxon, "_NY"))
  
  # Prepare final dataset
  print("Prepare final dataset")
  dat <- data.frame(
    y = y_combined,
    Sample = rownames(y_combined),
    Trait = metadata_combined$Trait
  )
  
  # Reshape to multivariate format
  dat$trait_group <- factor(dat$Trait)
  dat$Sample <- factor(dat$Sample)
  dat$maternal = as.factor(metadata_combined$dam[match(dat$Sample,metadata_combined$Sample)])
  dat$cohousing = as.factor(metadata_combined$cage[match(dat$Sample,metadata_combined$Sample)])
  
  # Fit bivariate model
  print("Fit the full bivariate model")
  full_mod <- mmer(
    y      ~ Trait,
    random = ~ vsr(usr(trait_group), Sample, Gu = grm_expanded) + vsr(dsr(trait_group), maternal) + vsr(dsr(trait_group), cohousing),
    rcov   = ~ units,
    data   = dat
  )
    
  # Extract genetic correlation
  print("Extract genetic correlation")
  cov_full <- matrix(c(
    full_mod$sigma[[1]],   # Var(MI)
    full_mod$sigma[[2]],   # Cov(NY,MI)
    full_mod$sigma[[2]],   # Cov(MI,NY)
    full_mod$sigma[[3]]    # Var(NY)
  ), nrow = 2, byrow = TRUE)
  gcor_est <- cov2cor(cov_full)[1, 2]
  
  # Extract final log-likelihoods
  logL_full <- full_mod$monitor[1, ncol(full_mod$monitor)]
  
  # — Fit reduced model (covariance = 0) —
  red_mod <- mmer(
    y      ~ Trait,
    random = ~ vsr(dsr(trait_group), Sample, Gu = grm_expanded) + vsr(dsr(trait_group), maternal) + vsr(dsr(trait_group), cohousing),
    rcov   = ~ units,
    data   = dat
  )
  logL_red  <-  red_mod$monitor[1, ncol( red_mod$monitor)]
  
  # LRT p-value (1 df)
  LRT_stat <- 2 * (logL_full - logL_red)
  p_val    <- pchisq(LRT_stat, df = 1, lower.tail = FALSE)
  
  cor_results=data.frame(
    Taxon               = taxon,
    Genetic_Correlation = gcor_est,
    p_value             = p_val,
    stringsAsFactors    = FALSE
  )
  print(cor_results)
  
  # Store results
  print("Store results")
  genetic_cor_results <- rbind(
    genetic_cor_results,
    cor_results
  )
}

print("Multiple-testing correction")
genetic_cor_results$bonf_pvalue = genetic_cor_results$p_value*dim(genetic_cor_results)[1]
genetic_cor_results$qvalue =p.adjust(genetic_cor_results$p_value, method = "fdr")
genetic_cor_results <- genetic_cor_results %>% 
  mutate( Significance = case_when(
    qvalue < 0.1 & bonf_pvalue < 0.05 ~ "Significant (Bonferroni)",
    qvalue < 0.1 & bonf_pvalue >= 0.05 ~ "Significant (FDR)",
    TRUE ~ "Non-significant"))
genetic_cor_results = genetic_cor_results[order(genetic_cor_results$p_value, decreasing = T),]

print("Saving")
save(genetic_cor_results,file="/users/abaud/fmorillo/paper_figures/genetic_correlations_cohorts.RData")

load("/users/abaud/fmorillo/paper_figures/genetic_correlations_cohorts.RData")
