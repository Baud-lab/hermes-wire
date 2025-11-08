#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table); library(GENESIS); library(SeqArray); library(SeqVarTools)
})

# --- args ---
args <- commandArgs(trailingOnly=TRUE)
kv <- function(k) sub(paste0(".*",k,"="),"",grep(paste0("^",k,"="),args,val=TRUE))
gds <- kv("--gds")
scan_lu <- kv("--scan_lookup")
scan_ids <- kv("--scan_ids_test"); label_ids <- kv("--label_ids_test")
resid_rda <- kv("--resid_rda"); beta_rda <- kv("--beta_rda"); clusters_rda <- kv("--clusters_rda")
phen_rda <- kv("--phen_rda"); meta <- kv("--metadata")
grm <- kv("--grm_rdata"); dam<-kv("--srm_dam"); cage<-kv("--srm_cage")
loco <- kv("--loco_rds"); bins<-kv("--loco_bins"); Ns<-kv("--loco_Ns"); ids<-kv("--loco_ids")
genus_prefix <- kv("--genus_prefix")
exp_genus <- kv("--exposure_genus"); out_genus <- kv("--outcome_genus")
extra_out <- kv("--extra_outcomes")
snps_tsv <- kv("--snps_tsv"); outfile <- kv("--out")

# --- load residual traits ---
load(resid_rda)      # residuals_qned_counts_objs (matrix-like)
resid_last <- residuals_qned_counts_objs[[length(residuals_qned_counts_objs)]]
load(beta_rda)       # should contain matrix with row named params.med_beta_trait_name
load(clusters_rda)   # same idea
load(phen_rda)       # expected a data.frame with a column = params.med_glucose_column

# labels/test IDs
lab <- readLines(label_ids); scn <- readLines(scan_ids)

# trait getter from residual RDatas
get_row <- function(M, nm) {
  v <- as.numeric(M[nm, ])
  names(v) <- colnames(M); v[lab]
}

# extra outcomes list
outs <- if (nzchar(extra_out)) strsplit(extra_out,",")[[1]] else character()
outs <- trimws(unique(tolower(outs)))

# exposure + outcomes names
exp_trait <- paste0(genus_prefix, exp_genus)

out_traits <- c()
if ("genus_other" %in% outs) out_traits <- c(out_traits, paste0(genus_prefix, out_genus))
if ("beta"        %in% outs) out_traits <- c(out_traits, "BETA_SENTINEL")     # placeholder
if ("guild3"      %in% outs) out_traits <- c(out_traits, "GUILD3_SENTINEL")   # placeholder
if ("glucose"     %in% outs) out_traits <- c(out_traits, "GLUCOSE_SENTINEL")  # placeholder
out_traits <- unique(out_traits)

# compose list of trait vectors
traits_list <- list()
traits_src  <- list()

traits_list[[exp_trait]] <- get_row(resid_last, exp_trait); traits_src[[exp_trait]] <- "resid"
for (ot in out_traits) {
  if (ot == "BETA_SENTINEL") {
    # 'beta' row is in beta_rda (assume object is last matrix)
    obj <- mget(ls()[grepl("residuals", ls())], ifnotfound=list(), inherits=TRUE)
    # safer: read from beta_rda directly again
    load(beta_rda)
    beta_obj <- get(ls()[1])
    traits_list[[ot]] <- get_row(beta_obj, Sys.getenv("MED_BETA_TRAIT_NAME", unset="beta__TD_PCoA1"))
    traits_src[[ot]]  <- "beta_rda"
  } else if (ot == "GUILD3_SENTINEL") {
    load(clusters_rda)
    clu_obj <- get(ls()[1])
    traits_list[[ot]] <- get_row(clu_obj, Sys.getenv("MED_GUILD3_NAME", unset="cluster__3"))
    traits_src[[ot]]  <- "clusters_rda"
  } else if (ot == "GLUCOSE_SENTINEL") {
    # phen_rda is a data.frame with row 'sample' or similar; accept vector by label order
    load(phen_rda)
    ph <- get(ls()[ls()!="residuals_qned_counts_objs"][1])
    col <- Sys.getenv("MED_GLUCOSE_COL", unset="glucose_at_dissection")
    stopifnot(col %in% names(ph))
    x <- ph[[col]]; names(x) <- ph$sample
    traits_list[[ot]] <- x[lab]; traits_src[[ot]] <- "phen_rda"
  } else {
    # genus_other case (already built above via resid_last)
    traits_list[[ot]] <- get_row(resid_last, ot)
    traits_src[[ot]]  <- "resid"
  }
}

# --- SNPs to test
S <- fread(snps_tsv)  # gene_name, SNP, CHR, BP, A1, A2

# --- Build null models and run tests (mirror your GENESIS setup)
# NOTE: Reuse your own helper to construct ScanAnnotationDataFrame + null model with LOCO GRMs
# For brevity, assuming a helper function get_null_model(trait_vector, ...) exists.
get_null_model <- function(y) {
  stop("Implement null model setup as in your run_genesis_genus.R (LOCO + SRMs)")
}

gdsf <- seqOpen(gds); on.exit(seqClose(gdsf))
res_list <- list()

# exposure first (to carry EAF etc. for harmonisation)
exp_nm <- exp_trait
null_exp <- get_null_model(traits_list[[exp_nm]])
for (gn in unique(S$gene_name)) {
  ids <- S[gene_name==gn, SNP]
  if (!length(ids)) next
  seqSetFilter(gdsf, variant.id=NULL, sample.id=scn, verbose=FALSE)
  seqSetFilterAnnot(gdsf, include=list(id=ids))
  AE <- assocTestSingle(null_exp, test="Score", verbose=FALSE)
  E <- as.data.table(AE)
  setnames(E, tolower(names(E)))
  E <- E[, .(SNP=id, A1=alt, A2=ref, maf=maf, nsamp=n.obs,
             beta.exposure=est, se.exposure=se, gene_name=gn)]
  # outcomes
  for (ot in out_traits) {
    null_out <- get_null_model(traits_list[[ot]])
    AO <- assocTestSingle(null_out, test="Score", verbose=FALSE)
    O <- as.data.table(AO)
    setnames(O, tolower(names(O)))
    O <- O[, .(SNP=id, beta.outcome=est, se.outcome=se, outcome=ot)]
    M <- merge(E, O, by="SNP", all.x=TRUE)
    M[, eaf := maf]
    res_list[[length(res_list)+1]] <- M
  }
}
RES <- rbindlist(res_list, fill=TRUE)
fwrite(RES, sub("\\.gz$","", outfile), sep="\t"); system(paste("gzip -f", shQuote(sub("\\.gz$","", outfile))))
