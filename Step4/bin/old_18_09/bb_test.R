#!/usr/bin/env Rscript
# bb_across_genera.R
suppressPackageStartupMessages({ library(data.table);
  library(optparse);
  library(pheatmap);
  library(rstatix);
  library(dplyr)})

opt_list <- list(
  make_option("--corr_prev", type="character", default="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step4/results_prev/phase7_pgls/correlations_by_gene__Prevotella.tsv"),
  make_option("--corr_bac",  type="character", default="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step4/results_bac/phase7_pgls/correlations_by_gene__Bacteroides.tsv"),
  make_option("--geno_prev", type="character", default="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step4/results_prev/phase8_synthesis_simple/Prevotella/selected_genes.RData"),
  make_option("--geno_bac",  type="character", default="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step4/results_bac/phase8_synthesis_simple/Bacteroides/selected_genes.RData"),
  make_option("--select_by", type="character", default="qvalue<0.10", help="selection rule for families (e.g., 'qvalue<0.10' or 'p_used<0.05')"),
  make_option("--out",  type="character", default="bb_combined.tsv")
)
opt <- parse_args(OptionParser(option_list=opt_list))

rbind_safe <- function(...) rbindlist(list(...), use.names=TRUE, fill=TRUE)
x1 <- fread(opt$corr_prev, sep="\t")
x2 <- fread(opt$corr_bac,  sep="\t")
dt <- rbind_safe(x1, x2)
setnames(dt, tolower(names(dt)))
stopifnot(all(c("genus","gene","p_used") %in% names(dt)))
if (!("qvalue" %in% names(dt))) dt[, qvalue := p.adjust(p_used, method="fdr"), by=genus]

# --- select families (genera) ---
sel <- try(dt[, any(eval(parse(text=opt$select_by))), by=genus], silent=TRUE)
if (inherits(sel,"try-error")) stop("Bad --select_by expression: ", opt$select_by)
sel_genera <- sel[V1==TRUE, genus]
m <- uniqueN(dt$genus); R <- length(sel_genera)
dt[, selected_family := genus %in% sel_genera]

# --- BB reweight inside selected families ---
# build the per-genus BB results for selected genera
bb <- dt[selected_family == TRUE, {
  w    <- m / max(1L, R)
  p_bb <- pmin(1, p_used * w)
  .(bb_p = p_bb, bb_q = p.adjust(p_bb, method = "BH"))
}, by = genus]

# update dt by joining bb INTO dt (bb is the i= table)
dt[bb, on = .(genus), `:=`(bb_p = i.bb_p, bb_q = i.bb_q)]

dt <- dt %>% 
  mutate( BB_Sig = case_when(
    bb_q < 0.1 ~ "Significant (FDR)",
    TRUE ~ "Non-significant"))

#fwrite(dt, opt$out, sep="\t")
#cat(sprintf("[BB] Families m=%d, selected R=%d, wrote %s\n", m, R, opt$out))

keep_genes = keep_genes=dt$gene[dt$BB_Sig=="Significant (FDR)"]

#######################################

load(opt$geno_prev)
geno_Prevotella=genotypes
correlations_prev=correlations
load(opt$geno_bac)
geno_Bacteroides=genotypes
correlations_bac=correlations
correlations=rbind(correlations_prev,correlations_bac)
correlations=correlations[correlations$gene %in% keep_genes,]

#genes=unique(correlations$Gene_name)
#dup_genes=correlations$gene[which(duplicated(correlations$gene))]

#########################

x1$Bioprocess=correlations_prev$Bioprocess[match(x1$gene,correlations_prev$gene)]
x1$Category=correlations_prev$Category[match(x1$gene,correlations_prev$gene)]

x2$Bioprocess=correlations_bac$Bioprocess[match(x2$gene,correlations_bac$gene)]
x2$Category=correlations_bac$Category[match(x2$gene,correlations_bac$gene)]

count(x1,Bioprocess,Direction)
count(x2,Bioprocess,Direction)

write.table(x1,file="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step4/results_prev/results_prev.tsv",quote=F,sep="\t",col.names=T,row.names=F)
write.table(x2,file="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step4/results_bac/results_bac.tsv",quote=F,sep="\t",col.names=T,row.names=F)

#########################

residuals_ALL=data.frame()

# Read count tables after CLR transformation
print("Read residuals")
load("~/NextFlow/Microbiome_Profiling/step2/Output_Tax_Mat_Net/Matrix_Processing/residuals_qned_counts.RData")
residuals_genus=residuals_qned_counts_objs[[5]]

# Load beta
print("Read beta diversity")
load("/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Tax_Mat_Net/Diversity/residuals_qned_counts_beta.RData")
beta=residuals_qned_counts_beta_objs[[1]]
residuals_beta=beta[5,,drop=F]
rownames(residuals_beta)[1]="Beta_PD_PC1"
residuals_ALL=rbind(residuals_ALL,residuals_beta)

## Load cluster residuals
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Tax_Mat_Net/Cluster_Analyses/residuals_qned_counts_clusters.RData")
all_traits_cluster=residuals_qned_counts_clusters_objs[[1]]
all_traits_cluster=all_traits_cluster[grepl("cluster__3",rownames(all_traits_cluster)),,drop=F]
rownames(all_traits_cluster)="Guild_3"
residuals_ALL=rbind(residuals_ALL,all_traits_cluster)

genera=c("Prevotella","Bacteroides")
for (genus in genera){
  print(paste0("Genus: ",genus))
  this_residuals_genus=residuals_genus[rownames(residuals_genus)==paste0("g__",genus),,drop=F]
  rownames(this_residuals_genus)=gsub("g__","",rownames(this_residuals_genus))
  this_residuals_ALL=rbind(residuals_ALL,this_residuals_genus)
  y_total=get(paste("geno",genus,sep="_"))
  print(paste0("Columns before: ",(ncol(y_total)-1)))
  # We'll create a vector to store the names of the columns to remove.
  cols_to_remove <- c()
  # Loop through all columns starting from the second one (index 2).
  # The '2:ncol(y_total)' part ensures we skip the first column.
  for (i in 2:ncol(y_total)) {
    # Get the number of unique values in the current column.
    n_unique <- length(unique(y_total[, i]))
    # Check if the number of unique values is less than 3.
    if (n_unique < 3) {
      # If the condition is met, add the column's name to our removal list.
      cols_to_remove <- c(cols_to_remove, colnames(y_total)[i])
    }
  }
  # Remove the identified columns from the data frame.
  y_total <- y_total[, !(colnames(y_total) %in% cols_to_remove)]
  # Print the result to see which columns were kept.
  print(paste0("Columns after: ",(ncol(y_total)-1)))
  genes=unique(colnames(y_total)[-1])
  genes=sort(genes)
  colnames(this_residuals_ALL)=sapply(strsplit(colnames(this_residuals_ALL), "_"),"[",1)
  this_residuals_ALL=this_residuals_ALL[,colnames(this_residuals_ALL) %in% y_total$sample]
  traits=rownames(this_residuals_ALL)
  for (trait in traits){
    print(paste0("Trait: ",trait))
    residuals=as.numeric(this_residuals_ALL[rownames(this_residuals_ALL)==trait,])
    y_total[[trait]]=residuals
  }
  y_total=na.omit(y_total)
  genotype_total=data.frame()
  samples=y_total$sample
  pdf(paste0("/users/abaud/fmorillo/paper_figures/putative_genes_genotypes_",genus,".pdf"))
  for (gene in genes){
    for (trait in traits){
      print(paste0("Trait: ",trait," / Gene: ",gene))
      comb_filtered=data.frame(sample=samples,
                               gene=gene,
                               trait=trait,
                               genotypes=y_total[,colnames(y_total)==gene,drop=F],
                               residuals=y_total[,colnames(y_total)==trait,drop=F])
      colnames(comb_filtered)[4]="genotype"
      colnames(comb_filtered)[5]="residuals"
      if (trait == "Beta_PD_PC1"){
        y_text="Beta diversity (PCo1)"
        trait_text="Beta diversity"
      } else if (trait == "alpha__PD_q2"){
        y_text="Alpha diversity (Rao's quadratic index)"
        trait_text="Alpha diversity"
      } else if (trait == "Guild_3"){
        trait_text="Guild 3"
        y_text="Pseudo-relative abundance"
      } else {
        y_text="Pseudo-relative abundance"
        trait_text <- trait
      }
      p_value1 <- wilcox.test(comb_filtered$residuals[comb_filtered$genotype == "2"], 
                              comb_filtered$residuals[comb_filtered$genotype == "1"])$p.value
      p_value2 <- wilcox.test(comb_filtered$residuals[comb_filtered$genotype == "1"], 
                              comb_filtered$residuals[comb_filtered$genotype == "0"])$p.value
      comb_filtered_v2=comb_filtered
      comb_filtered_v2$genotype=as.character(comb_filtered_v2$genotype)
      comb_filtered_v2$genotype[comb_filtered_v2$genotype == "0"]="BB"
      comb_filtered_v2$genotype[comb_filtered_v2$genotype == "1"]="AB"
      comb_filtered_v2$genotype[comb_filtered_v2$genotype == "2"]="AA"
      maxv=max(comb_filtered_v2$residuals)
      order=c("AA", "AB", "BB")
      myplot<-ggplot(comb_filtered_v2, aes(x=factor(genotype, levels = order), y=residuals, fill=genotype)) +
        geom_boxplot() +
        xlab("Genotypes") +
        ylab(y_text) +
        ggtitle(paste0("Trait: ",trait_text," / Gene: ",gene)) +
        theme_bw() +
        scale_fill_manual(
          values = c(
            "AA" = "#0077be",
            "AB" = "lightblue",
            "BB" = "cyan"
          ), guide="none")+
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
        #scale_x_discrete(guide = guide_axis(angle = 60)) +
        geom_signif(comparisons = list(c("AA","AB"),c("AB","BB")),
                    y_position = maxv + 0.2,
                    annotations = sprintf("p = %.2e", c(p_value1, p_value2)),
                    map_signif_level = TRUE)
      print(myplot)
      median_residuals_df <- comb_filtered %>%
        group_by(genotype) %>%
        summarise(median_residuals = median(residuals, na.rm = TRUE))
      median_residuals_df = median_residuals_df[order(median_residuals_df$median_residuals, decreasing = T),]
      median_residuals_df$gene=gene
      median_residuals_df$trait=trait
      if (trait=="Beta_PD_PC1"){
        if (length(median_residuals_df$genotype)<3){
          median_residuals_df=median_residuals_df[c(1,2),c(1,3,4),drop=F]
        } else {
          median_residuals_df=median_residuals_df[c(1,3),c(1,3,4),drop=F]
        }
        median_residuals_df$trait[1]=paste(trait,"Ent2",sep="_")
        median_residuals_df$trait[2]=paste(trait,"Ent1",sep="_")
      } else {
        median_residuals_df=median_residuals_df[1,c(1,3,4),drop=F]
      }
      genotype_total=rbind(genotype_total,median_residuals_df)
    }
  }
  dev.off()
  genotypes_total <- genotype_total %>%
    mutate(genotype = as.character(genotype))
  genotypes_wide <- genotype_total %>%
    pivot_wider(
      names_from = gene,
      values_from = genotype
    )
  genotypes_wide=t(genotypes_wide)
  colnames(genotypes_wide)=genotypes_wide[1,]
  genotypes_wide=genotypes_wide[-1,]
  genotypes_wide=as.data.frame(genotypes_wide)
  
  # Calculate Cramer's V for all pairs
  cramers_v_results <- data.frame(
    Pair = character(),
    Cramers_V = numeric(),
    Chi_sq_pvalue = numeric(),
    stringsAsFactors = FALSE
  )
  # Guild_3 vs. Beta_PD_PC1_Ent1
  v_guild3_ent1 <- cramer_v(genotypes_wide$Guild_3, genotypes_wide$Beta_PD_PC1_Ent1)
  chi_sq_guild3_ent1 <- chisq.test(genotypes_wide$Guild_3, genotypes_wide$Beta_PD_PC1_Ent1)
  cramers_v_results <- rbind(cramers_v_results, data.frame(Pair = "Guild_3_vs_Ent1", Cramers_V = v_guild3_ent1, Chi_sq_pvalue=chi_sq_guild3_ent1$p.value))
  # Guild_3 vs. Beta_PD_PC1_Ent2
  v_guild3_ent2 <- cramer_v(genotypes_wide$Guild_3, genotypes_wide$Beta_PD_PC1_Ent2)
  chi_sq_guild3_ent2 <- chisq.test(genotypes_wide$Guild_3, genotypes_wide$Beta_PD_PC1_Ent2)
  cramers_v_results <- rbind(cramers_v_results, data.frame(Pair = "Guild_3_vs_Ent2", Cramers_V = v_guild3_ent2, Chi_sq_pvalue=chi_sq_guild3_ent2$p.value))
  # Genus vs. Beta_PD_PC1_Ent1
  v_genus_ent1 <- cramer_v(genotypes_wide[,4], genotypes_wide$Beta_PD_PC1_Ent1)
  chi_sq_genus_ent1 <- chisq.test(genotypes_wide[,4], genotypes_wide$Beta_PD_PC1_Ent1)
  cramers_v_results <- rbind(cramers_v_results, data.frame(Pair = paste0(genus,"_vs_Ent1"), Cramers_V = v_genus_ent1, Chi_sq_pvalue=chi_sq_genus_ent1$p.value))
  # Genus vs. Beta_PD_PC1_Ent2
  v_genus_ent2 <- cramer_v(genotypes_wide[,4], genotypes_wide$Beta_PD_PC1_Ent2)
  chi_sq_genus_ent2 <- chisq.test(genotypes_wide[,4], genotypes_wide$Beta_PD_PC1_Ent2)
  cramers_v_results <- rbind(cramers_v_results, data.frame(Pair = paste0(genus,"_vs_Ent2"), Cramers_V = v_genus_ent2, Chi_sq_pvalue=chi_sq_genus_ent2$p.value))
  # bac vs. Guild3
  v_genus_guild3 <- cramer_v(genotypes_common[,4], genotypes_common$Guild_3)
  chi_sq_genus_guild3 <- chisq.test(genotypes_common[,4], genotypes_common$Guild_3)
  cramers_v_results <- rbind(cramers_v_results, data.frame(Pair = paste0(genus,"_vs_Guild_3"), Cramers_V = v_genus_guild3, Chi_sq_pvalue=chi_sq_genus_guild3$p.value))
  # Print the results
  print(cramers_v_results)
  
  save(genotypes_wide,
       cramers_v_results,
       y_total,
       file=paste0("/users/abaud/fmorillo/paper_figures/genotypes_wide_",genus,".RData"))
}


####################################################
load("~/paper_figures/genotypes_widePrevotella.RData")
genotypes_wide_prev=genotypes_wide
load("~/paper_figures/genotypes_wideBacteroides.RData")
genotypes_wide_bac=genotypes_wide

load(opt$geno_prev)
correlations_prev=correlations
load(opt$geno_bac)
correlations_bac=correlations
correlations=rbind(correlations_prev,correlations_bac)
dup_genes=correlations$gene[which(duplicated(correlations$gene))]
dup_genes=correlations$Gene_name[match(dup_genes,correlations$gene)]

genotypes_wide_bac_filt=genotypes_wide_bac[rownames(genotypes_wide_bac) %in% dup_genes,]
genotypes_wide_prev_filt=genotypes_wide_prev[rownames(genotypes_wide_prev) %in% dup_genes,]
#head(genotypes_wide_bac_filt)
#head(genotypes_wide_prev_filt)
genotypes_common=genotypes_wide_bac_filt
genotypes_common$Prevotella=genotypes_wide_prev_filt$Prevotella

# Calculate Cramer's V for all pairs
cramers_v_results <- data.frame(
  Pair = character(),
  Cramers_V = numeric(),
  Chi_sq_pvalue = numeric(),
  stringsAsFactors = FALSE
)
# Guild_3 vs. Beta_PD_PC1_Ent1
v_guild3_ent1 <- cramer_v(genotypes_common$Guild_3, genotypes_common$Beta_PD_PC1_Ent1)
chi_sq_guild3_ent1 <- chisq.test(genotypes_common$Guild_3, genotypes_common$Beta_PD_PC1_Ent1)
cramers_v_results <- rbind(cramers_v_results, data.frame(Pair = "Guild_3_vs_Ent1", Cramers_V = v_guild3_ent1, Chi_sq_pvalue=chi_sq_guild3_ent1$p.value))
# Guild_3 vs. Beta_PD_PC1_Ent2
v_guild3_ent2 <- cramer_v(genotypes_common$Guild_3, genotypes_common$Beta_PD_PC1_Ent2)
chi_sq_guild3_ent2 <- chisq.test(genotypes_common$Guild_3, genotypes_common$Beta_PD_PC1_Ent2)
cramers_v_results <- rbind(cramers_v_results, data.frame(Pair = "Guild_3_vs_Ent2", Cramers_V = v_guild3_ent2, Chi_sq_pvalue=chi_sq_guild3_ent2$p.value))
# bac vs. Beta_PD_PC1_Ent1
v_bac_ent1 <- cramer_v(genotypes_common$Bacteroides, genotypes_common$Beta_PD_PC1_Ent1)
chi_sq_bac_ent1 <- chisq.test(genotypes_common$Bacteroides, genotypes_common$Beta_PD_PC1_Ent1)
cramers_v_results <- rbind(cramers_v_results, data.frame(Pair = paste0("Bacteroides_vs_Ent1"), Cramers_V = v_bac_ent1, Chi_sq_pvalue=chi_sq_bac_ent1$p.value))
# bac vs. Beta_PD_PC1_Ent2
v_bac_ent2 <- cramer_v(genotypes_common$Bacteroides, genotypes_common$Beta_PD_PC1_Ent2)
chi_sq_bac_ent2 <- chisq.test(genotypes_common$Bacteroides, genotypes_common$Beta_PD_PC1_Ent2)
cramers_v_results <- rbind(cramers_v_results, data.frame(Pair = paste0("Bacteroides_vs_Ent2"), Cramers_V = v_bac_ent2, Chi_sq_pvalue=chi_sq_bac_ent2$p.value))

# prev vs. Beta_PD_PC1_Ent1
v_prev_ent1 <- cramer_v(genotypes_common$Prevotella, genotypes_common$Beta_PD_PC1_Ent1)
chi_sq_prev_ent1 <- chisq.test(genotypes_common$Prevotella, genotypes_common$Beta_PD_PC1_Ent1)
cramers_v_results <- rbind(cramers_v_results, data.frame(Pair = paste0("Prevotella_vs_Ent1"), Cramers_V = v_prev_ent1, Chi_sq_pvalue=chi_sq_prev_ent1$p.value))
# prev vs. Beta_PD_PC1_Ent2
v_prev_ent2 <- cramer_v(genotypes_common$Prevotella, genotypes_common$Beta_PD_PC1_Ent2)
chi_sq_prev_ent2 <- chisq.test(genotypes_common$Prevotella, genotypes_common$Beta_PD_PC1_Ent2)
cramers_v_results <- rbind(cramers_v_results, data.frame(Pair = paste0("Prevotella_vs_Ent2"), Cramers_V = v_prev_ent2, Chi_sq_pvalue=chi_sq_prev_ent2$p.value))
# prev vs. Guild3
v_prev_guild3 <- cramer_v(genotypes_common$Prevotella, genotypes_common$Guild_3)
chi_sq_prev_guild3 <- chisq.test(genotypes_common$Prevotella, genotypes_common$Guild_3)
cramers_v_results <- rbind(cramers_v_results, data.frame(Pair = paste0("Prevotella_vs_Guild_3"), Cramers_V = v_prev_guild3, Chi_sq_pvalue=chi_sq_prev_guild3$p.value))
# bac vs. Guild3
v_bac_guild3 <- cramer_v(genotypes_common$Bacteroides, genotypes_common$Guild_3)
chi_sq_bac_guild3 <- chisq.test(genotypes_common$Bacteroides, genotypes_common$Guild_3)
cramers_v_results <- rbind(cramers_v_results, data.frame(Pair = paste0("Bacteroides_vs_Guild_3"), Cramers_V = v_bac_guild3, Chi_sq_pvalue=chi_sq_bac_guild3$p.value))
# prev vs. prev
v_prev_bac <- cramer_v(genotypes_common$Prevotella, genotypes_common$Bacteroides)
chi_sq_prev_bac <- chisq.test(genotypes_common$Prevotella, genotypes_common$Bacteroides)
cramers_v_results <- rbind(cramers_v_results, data.frame(Pair = paste0("Prevotella_vs_Bacteroides"), Cramers_V = v_prev_bac, Chi_sq_pvalue=chi_sq_prev_bac$p.value))

# Print the results
print(cramers_v_results)

##################################
