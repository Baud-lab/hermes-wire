#!/usr/bin/env Rscript
## meta_beta_cf_across_folds.R
## Combine per-fold beta_cf tables into a cross-fit meta table.
## Build a per-gene CONSENSUS sentinel from per-fold sentinel TSVs.
## Usage (preferred):
##   meta_beta_cf_across_folds.R \
##     --inputs-list inputs.txt \
##     --sentinels-list sentinels.txt \
##     --out crossfit.tsv \
##     --consensus consensus.tsv
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
})

LOG  <- function(...) cat("[meta_cf]", paste(..., collapse=" "), "\n")
STOP <- function(...) { cat("FATAL:", paste(..., collapse=" "), "\n"); quit(save="no", status=1) }

opt <- parse_args(OptionParser(option_list = list(
  make_option("--inputs",        type="character", default="", help="CSV list of per-fold beta files (legacy)"),
  make_option("--sentinels",     type="character", default="", help="CSV list of per-fold sentinel files (legacy)"),
  make_option("--inputs-list",   type="character", default="", help="Text file; one beta file per line"),
  make_option("--sentinels-list",type="character", default="", help="Text file; one sentinel file per line"),
  make_option("--out",           type="character", default="meta_beta_cf.tsv", help="Cross-fit meta output"),
  make_option("--consensus",     type="character", default="consensus_sentinels.tsv", help="Consensus sentinel output")
)))

## ---------- Read file lists ----------
read_list <- function(csv, lst) {
  ff <- character()
  if (nzchar(lst)) {
    if (!file.exists(lst)) STOP("--inputs-list/--sentinels-list file not found: ", lst)
    ff <- trimws(readLines(lst))
  }
  if (nzchar(csv)) {
    ff2 <- strsplit(csv, ",", fixed=TRUE)[[1]]
    ff2 <- trimws(ff2)
    ff <- c(ff, ff2)
  }
  ff <- unique(ff[nzchar(ff)])
  ff
}

in_files   <- read_list(opt$inputs,    opt$`inputs-list`)
sent_files <- read_list(opt$sentinels, opt$`sentinels-list`)
if (!length(in_files))   STOP("No input beta files provided")
if (!length(sent_files)) STOP("No sentinel files provided")

pick_col <- function(nms, wanted) {
  low <- tolower(nms); wanted <- tolower(wanted)
  hit <- which(low %in% wanted); if (!length(hit)) return(NA_character_)
  nms[hit[1]]
}

## ---------- 1) Read & stack per-fold beta tables ----------
all_rows <- list(); idx <- 0L
for (f in in_files) {
  if (!file.exists(f)) { LOG("WARN missing beta:", f); next }
  dt <- tryCatch(fread(f), error=function(e) NULL)
  if (is.null(dt) || !nrow(dt)) { LOG("WARN empty/unreadable beta:", f); next }
  
  nms <- names(dt)
  col_trait <- pick_col(nms, c("trait","species","taxon","taxa"))
  col_chr   <- pick_col(nms, c("chr","chrom","chromosome"))
  col_pos   <- pick_col(nms, c("pos","position","bp"))
  col_beta  <- pick_col(nms, c("beta_meta","Est","beta","effect","beta_mean"))
  col_se    <- pick_col(nms, c("se_meta","Est.SE","se","stderr","StdErr"))
  col_p     <- pick_col(nms, c("p_meta","Score.pval","pval","p"))
  
  if (any(is.na(c(col_trait,col_chr,col_pos,col_beta)))) {
    LOG("WARN missing essential columns in", basename(f)); next
  }
  
  out <- data.table(
    trait = as.character(dt[[col_trait]]),
    chr   = as.character(dt[[col_chr]]),
    pos   = suppressWarnings(as.integer(dt[[col_pos]])),
    beta  = suppressWarnings(as.numeric(dt[[col_beta]]))
  )
  out[, se := if (!is.na(col_se)) suppressWarnings(as.numeric(dt[[col_se]])) else NA_real_]
  out[, p  := if (!is.na(col_p))  suppressWarnings(as.numeric(dt[[col_p]]))  else NA_real_]
  
  idx <- idx + 1L
  b <- basename(f)
  foldid <- sub(".*fold\\s*([0-9]+).*", "\\1", b, perl=TRUE)
  if (identical(foldid, b) || !nzchar(foldid)) foldid <- as.character(idx)
  out[, fold := foldid]
  
  setkey(out, trait, chr, pos)
  out <- out[, .SD[order(ifelse(is.na(p), Inf, p))][1], by=key(out)]  # de-dup within fold
  all_rows[[length(all_rows)+1L]] <- out
}
if (!length(all_rows)) STOP("No valid beta inputs after reading")
DT <- rbindlist(all_rows, use.names=TRUE, fill=TRUE)
LOG("stacked rows:", nrow(DT),
    "| unique trait×site:", DT[, uniqueN(paste(trait, chr, pos))],
    "| folds seen:", paste(unique(DT$fold), collapse=","))

## ---------- 2) Read per-fold sentinels & build CONSENSUS (1 SNP per gene) ----------
S_list <- list()
for (sf in sent_files) {
  if (!file.exists(sf)) { LOG("WARN missing sentinel:", sf); next }
  s <- tryCatch(fread(sf), error=function(e) NULL)
  if (is.null(s) || !nrow(s)) { LOG("WARN empty/unreadable sentinel:", sf); next }
  setDT(s); setnames(s, tolower(names(s)))
  nm_gene <- intersect(c("gene","gene_id"), names(s))[1]
  nm_chr  <- intersect(c("chr","chrom","chromosome"), names(s))[1]
  nm_pos  <- intersect(c("bp","pos","position"), names(s))[1]
  nm_p    <- intersect(c("p","pval","p.value","score.pval"), names(s))[1]
  if (is.na(nm_gene) || is.na(nm_chr) || is.na(nm_pos)) {
    LOG("WARN sentinel missing gene/chr/pos:", sf, "-> skipping file"); next
  }
  tmp <- unique(data.table(
    gene = as.character(s[[nm_gene]]),
    chr  = as.character(s[[nm_chr]]),
    pos  = suppressWarnings(as.integer(s[[nm_pos]])),
    p    = if (!is.na(nm_p)) suppressWarnings(as.numeric(s[[nm_p]])) else NA_real_
  ))
  S_list[[length(S_list)+1L]] <- tmp
}
if (!length(S_list)) STOP("No valid sentinel rows across files")
S <- rbindlist(S_list, use.names=TRUE, fill=TRUE)
S <- S[!is.na(gene) & !is.na(chr) & !is.na(pos)]
S_freq <- S[, .(n=.N), by=.(gene, chr, pos)]
S_medp <- S[, .(medp = suppressWarnings(median(p, na.rm=TRUE))), by=.(gene, chr, pos)]
CONS   <- merge(S_freq, S_medp, by=c("gene","chr","pos"), all.x=TRUE)
setorder(CONS, gene, -n, medp, chr, pos)
CONS_1 <- CONS[, .SD[1], by=gene]  # 1 SNP per gene (mode; tie-break by median p; then chr:pos)
LOG("consensus size (genes):", nrow(CONS_1))
fwrite(CONS_1[, .(gene, chr, pos, n, medp)], opt$consensus, sep="\t")
LOG("WROTE consensus:", opt$consensus)

## ---------- 3) Filter betas to CONSENSUS & add 'gene' ----------
before <- nrow(DT)
DT <- merge(DT, CONS_1[, .(gene, chr, pos)], by=c("chr","pos"))
LOG("After consensus filter: rows:", nrow(DT), "from", before,
    "| unique genes:", DT[, uniqueN(gene)],
    "| unique sites:", DT[, uniqueN(paste(chr,pos))])
if (nrow(DT) == 0L) STOP("No overlap between betas and consensus sentinels")

## ---------- 4) Meta across folds ----------
meta <- DT[, {
  n_folds <- .N
  if (all(is.na(se))) {
    beta_meta <- mean(beta, na.rm=TRUE); se_meta <- NA_real_; z_meta <- NA_real_; p_meta <- NA_real_
  } else {
    w <- 1/(se^2); w[!is.finite(w)] <- NA_real_
    if (all(is.na(w))) {
      beta_meta <- mean(beta, na.rm=TRUE); se_meta <- NA_real_; z_meta <- NA_real_; p_meta <- NA_real_
    } else {
      w[is.na(w)] <- 0
      beta_meta <- sum(w*beta)/sum(w)
      se_meta   <- sqrt(1/sum(w))
      z_meta    <- beta_meta / se_meta
      p_meta    <- 2*pnorm(abs(z_meta), lower.tail=FALSE)
    }
  }
  list(n_folds=n_folds, beta_meta=beta_meta, se_meta=se_meta, z_meta=z_meta, p_meta=p_meta,
       beta_mean=mean(beta, na.rm=TRUE), beta_sd=sd(beta, na.rm=TRUE),
       p_min_fold=suppressWarnings(min(p, na.rm=TRUE)))
}, by=.(gene, trait, chr, pos)]

setorder(meta, p_meta, na.last=TRUE)
top_n <- min(5L, nrow(meta))
if (top_n > 0L) {
  msg <- paste(sprintf("%s | %s chr%s:%s p=%.3g beta=%.3g (n=%d)",
                       meta$gene[1:top_n], meta$trait[1:top_n],
                       meta$chr[1:top_n], meta$pos[1:top_n],
                       meta$p_meta[1:top_n], meta$beta_meta[1:top_n],
                       meta$n_folds[1:top_n]), collapse=" | ")
  LOG("Top meta hits:", msg)
}

meta$bonf_pvalue = meta$p_meta*dim(meta)[1]
meta$qvalue =p.adjust(meta$p_meta, method = "fdr")
meta <- meta %>% 
  mutate( Significance = case_when(
    qvalue < 0.1 & bonf_pvalue < 0.05 ~ "Significant (Bonferroni)",
    qvalue < 0.1 & bonf_pvalue >= 0.05 ~ "Significant (FDR)",
    TRUE ~ "Non-significant"))
meta$Direction <- ifelse(meta$estimate > 0, "Positive", "Negative")

fwrite(meta, opt$out, sep="\t")
LOG("WROTE meta:", opt$out, "rows:", nrow(meta))
LOG("DONE")
