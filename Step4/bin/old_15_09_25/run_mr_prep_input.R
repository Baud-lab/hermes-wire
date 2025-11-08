#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(TwoSampleMR) })

args <- commandArgs(trailingOnly=TRUE)
kv <- function(k) sub(paste0(".*",k,"="),"",grep(paste0("^",k,"="),args,val=TRUE))
infile <- kv("--assoc_gz")
do_steiger <- as.logical(kv("--do_steiger"))
f_min <- as.numeric(kv("--f_min"))
ld_prune <- as.logical(kv("--ld_prune"))
ld_r2 <- as.numeric(kv("--ld_r2")); ld_kb <- as.numeric(kv("--ld_kb"))
prefix <- kv("--out_prefix")

DT <- fread(infile)
# Build exp/out tables
expd <- unique(DT[, .(SNP, effect_allele=A1, other_allele=A2, eaf, samplesize=nsamp,
                      beta=beta.exposure, se=se.exposure,
                      phenotype=paste0("EXP__", gene_name),
                      gene=gene_name)])
outd <- DT[, .(SNP, effect_allele=A1, other_allele=A2, eaf, samplesize=nsamp,
               beta=beta.outcome, se=se.outcome,
               phenotype=paste0("OUT__", outcome,"__",gene_name),
               gene=gene_name, outcome)]

exp <- format_data(as.data.frame(expd), type="exposure",
                   snp_col="SNP", beta_col="beta", se_col="se",
                   effect_allele_col="effect_allele", other_allele_col="other_allele",
                   eaf_col="eaf", samplesize_col="samplesize", phenotype_col="phenotype")
out <- format_data(as.data.frame(outd), type="outcome",
                   snp_col="SNP", beta_col="beta", se_col="se",
                   effect_allele_col="effect_allele", other_allele_col="other_allele",
                   eaf_col="eaf", samplesize_col="samplesize", phenotype_col="phenotype")

harm <- harmonise_data(exp, out, action=2)  # allele/strand harmonisation
# carry explicit gene/outcome
harm$gene    <- sub("^EXP__","",harm$id.exposure)
harm$outcome <- sub("^OUT__","",harm$id.outcome); harm$outcome <- sub("__.*$","",harm$outcome)

# Steiger
if (do_steiger) harm <- steiger_filtering(harm)

# F-stat filter
harm$F_stat <- (harm$beta.exposure / harm$se.exposure)^2
harm <- harm[harm$F_stat >= f_min]

# (Optional) light LD prune – with no LD matrix, at least drop duplicates
if (ld_prune) harm <- harm[!duplicated(harm$SNP)]

outfile <- paste0("mr_input__", prefix, ".tsv.gz")
fwrite(harm, outfile, sep="\t")
