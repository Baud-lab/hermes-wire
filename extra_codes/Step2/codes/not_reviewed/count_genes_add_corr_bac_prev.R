putative_genes <- read_xlsx("/users/abaud/fmorillo/paper_figures/putative_genes_final_full.xlsx")
putative_genes=putative_genes[,c(1:3)]
putative_genes$Gene=gsub("\\*","",putative_genes$Gene)
colnames(putative_genes)[3]="gene"


#correlation_prev <- read.delim("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step4/results_prev/phase7_pgls/correlations_by_gene__Prevotella.tsv")
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step4/results_prev/phase8_synthesis_simple/Prevotella/selected_genes.RData")
correlation_prev <- correlations
#correlation_bac <- read.delim("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step4/results_bac/phase7_pgls/correlations_by_gene__Bacteroides.tsv")
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step4/results_bac/phase8_synthesis_simple/Bacteroides/selected_genes.RData")
correlation_bac <- correlations
correlations_total=rbind(correlation_bac,correlation_prev)
correlations_total$gene=correlations_total$Gene_name
correlations_total=correlations_total[,1:4]
correlations_total <- correlations_total %>% 
  mutate( Direction = case_when(
    beta_used > 0 ~ "Positive",
    beta_used < 0 ~ "Negative"))
#correlations_total$genus=gsub("s__","",correlations_total$genus)


correlations_total$Bioprocess=putative_genes$Bioprocess[match(correlations_total$gene,putative_genes$gene)]
correlations_total$Category=putative_genes$Category[match(correlations_total$gene,putative_genes$gene)]

write.table(correlations_total,file="~/paper_figures/correlations_total_beta_h2.tsv",sep="\t",quote=F,col.names=T,row.names = F)

#correlations_total=correlations_total[correlations_total$gene!="Pax7",]


# packages
library(dplyr)
library(tidyr)
library(stringr)

## --- 0) Optional: ensure expected column names exist -------------------------
# putative_genes: Bioprocess, Category, gene
# correlations_total: genus, gene, All_Sig, Direction

## --- 1) Build per-gene flags from correlations_total -------------------------
gene_flags <- correlations_total %>%
  filter(genus %in% c("Bacteroides", "Prevotella")) %>%
  mutate(
    #Additive_Bacteroides    = genus == "Bacteroides",
    #Additive_Prevotella     = genus == "Prevotella",
    Correlation_Bacteroides = genus == "Bacteroides",# & Significance != "Significant",
    Correlation_Prevotella  = genus == "Prevotella"#  & ignificance != "Significant"
  ) %>%
  group_by(gene) %>%
  summarise(
    #Additive_Bacteroides    = any(Additive_Bacteroides),
    #Additive_Prevotella     = any(Additive_Prevotella),
    Correlation_Bacteroides = any(Correlation_Bacteroides),
    Correlation_Prevotella  = any(Correlation_Prevotella),
    .groups = "drop"
  )

## --- 2) Attach flags to putative_genes (fill missing with FALSE) ------------
pg <- putative_genes %>%
  left_join(gene_flags, by = "gene") %>%
  mutate(across(c(#Additive_Bacteroides, Additive_Prevotella,
                  Correlation_Bacteroides, Correlation_Prevotella),
                ~ replace_na(.x, FALSE)))

## --- 3) Category-level counts (inside each Bioprocess, sorted decreasing) ---
tab_category <- pg %>%
  group_by(Bioprocess, Category) %>%
  summarise(
    Count_genes             = n_distinct(gene),
    #Additive_Bacteroides    = sum(Additive_Bacteroides),
    #Additive_Prevotella     = sum(Additive_Prevotella),
    Correlation_Bacteroides = sum(Correlation_Bacteroides),
    Correlation_Prevotella  = sum(Correlation_Prevotella),
    .groups = "drop"
  )

## --- 4) Bioprocess totals (sorted decreasing) -------------------------------
tab_bioprocess <- pg %>%
  group_by(Bioprocess) %>%
  summarise(
    Count_genes             = n_distinct(gene),
    #Additive_Bacteroides    = sum(Additive_Bacteroides),
    #Additive_Prevotella     = sum(Additive_Prevotella),
    Correlation_Bacteroides = sum(Correlation_Bacteroides),
    Correlation_Prevotella  = sum(Correlation_Prevotella),
    .groups = "drop"
  ) %>%
  arrange(desc(Count_genes)) %>%
  mutate(Category = "(Total)")

## --- 5) Order categories within each bioprocess by their counts -------------
bp_rank <- tab_bioprocess %>%
  transmute(Bioprocess, BP_rank = row_number())  # order of bioprocess blocks

tab_category_ordered <- tab_category %>%
  left_join(bp_rank, by = "Bioprocess") %>%
  arrange(BP_rank, desc(Count_genes)) %>%
  select(-BP_rank)

tab_bioprocess_ordered <- tab_bioprocess %>%
  left_join(bp_rank, by = "Bioprocess")

## --- 6) Combine into a single “pivot-like” table with an indented label -----
pivot_like <- bind_rows(
  tab_bioprocess_ordered %>% mutate(Level = "Bioprocess"),
  tab_category_ordered %>% mutate(Level = "Category")
) %>%
  arrange(BP_rank, match(Level, c("Bioprocess", "Category")), desc(Count_genes)) %>%
  mutate(
    Row_Label = ifelse(Level == "Bioprocess",
                       Bioprocess,
                       paste0("    ", Category))
  ) %>%
  select(
    Row_Label,
    Count_genes,
    #Additive_Bacteroides, Additive_Prevotella,
    Correlation_Bacteroides, Correlation_Prevotella,
    Bioprocess, Category, Level  # keep for filtering if you want; drop if not needed
  )

## --- 7) (Optional) Add a Grand Total row ------------------------------------
grand_total <- tibble(
  Row_Label               = "Grand Total",
  Count_genes             = n_distinct(pg$gene),
  #Additive_Bacteroides    = sum(pg$Additive_Bacteroides),
  #Additive_Prevotella     = sum(pg$Additive_Prevotella),
  Correlation_Bacteroides = sum(pg$Correlation_Bacteroides),
  Correlation_Prevotella  = sum(pg$Correlation_Prevotella),
  Bioprocess = NA_character_,
  Category  = NA_character_,
  Level     = "Total"
)

pivot_like <- bind_rows(pivot_like, grand_total)


##############

library(dplyr)
library(forcats)
library(ggplot2)
library(scales)

# 1) Round to 3 decimals and rebuild the joined df
correlations_total <- correlations_total %>%
  mutate(r_value_round = round(beta_used, 3))

putative_genes_bac_prev <- correlations_total %>%
  filter(genus %in% c("Bacteroides", "Prevotella")) %>%
  left_join(putative_genes %>% select(gene, Bioprocess, Category), by = "gene") %>%
  mutate(Full_bioprocess = paste0(Bioprocess, " (", Category, ")"))

# 2) Keep only Significant rows (and de-dup if desired)
sig <- putative_genes_bac_prev# %>%
  #filter(All_Sig == "Significant", !is.na(Full_bioprocess))

# If you truly want up to 3 entries per gene (e.g., different Full_bioprocess),
# keep n = 3; otherwise set n = 1
sig <- sig %>%
  group_by(genus, gene) %>%
  slice_max(order_by = abs(r_value_round), n = 3, with_ties = FALSE) %>%
  ungroup()

# 3) Global symmetric limit so Bacteroides and Prevotella share the same scale
lim_global <- max(abs(sig$r_value_round), na.rm = TRUE)
lim_global <- ceiling(lim_global * 1000) / 1000  # round up to 0.001

make_barplot <- function(df, genus_name, lim = lim_fixed){
  df_g <- df %>% filter(genus == genus_name) %>%
    mutate(gene_ord = fct_reorder(gene, r_value_round, .desc = TRUE))
  if(genus_name=="Bacteroides"){
    lim_fixed=c(0, lim)
  } else {
    lim_fixed=c(-lim, lim)
  }
  ggplot(df_g, aes(x = gene_ord, y = r_value_round, fill = Bioprocess)) +
    geom_col(width = 0.6) +
    geom_hline(yintercept = 0, linewidth = 0.35) +
    coord_flip(clip = "off") +
    scale_y_continuous(
      limits = lim_fixed,
      breaks = seq(-lim, lim, by = 0.1),
      labels = number_format(accuracy = 0.01)
    ) +
    labs(x = NULL, y = "rho")+
    ggtitle(paste0("Spearman correlation (rho) between species additive effects and heritability: ",genus_name)) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "bottom",        # move below plot for more horizontal space
      legend.text = element_text(size = 10),
      legend.key.size = unit(0.5, "lines"),
      axis.text.y = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 12, face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5),
      plot.title = element_text(hjust = 0.5, size = 13, face = "bold")
    )
}

lim_fixed <- 10  # target y-limit
p_prev <- make_barplot(sig, "Prevotella",  lim_fixed)
lim_fixed <- 10  # target y-limit
p_bac  <- make_barplot(sig, "Bacteroides", lim_fixed)

pdf("~/paper_figures/correlation_herit_bac_prev.pdf")
print(p_prev)
print(p_bac)
dev.off()
