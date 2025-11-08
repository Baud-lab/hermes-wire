#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggsignif)
  library(forcats)
})

## ============================ CLI ============================ ##
option_list <- list(
  make_option(c("-i","--input_tsv"), type="character", default=NULL,
              help="Input TSV with species-level metrics (columns documented below)."),
  make_option(c("-o","--outpdf"), type="character", default="heritable_characteristics.pdf",
              help="Output PDF path [default: %default]"),
  make_option(c("-s","--subset_id"), type="character", default="All samples",
              help="Label to appear in plot titles [default: %default]"),
  make_option(c("--traits_phylum_tsv"), type="character", default=NULL,
              help="(Optional) Phylum abundance matrix (rows=phyla, cols=samples) for sample-level Bacteroidota vs beta."),
  make_option(c("--diversity_tsv"), type="character", default=NULL,
              help="(Optional) Sample-level diversity file with columns: sample, beta, alpha."),
  make_option(c("--seed"), type="integer", default=123,
              help="Random seed [default: %default]")
)
opt <- parse_args(OptionParser(option_list=option_list))
set.seed(opt$seed)

## ======================= Data requirements =================== ##
required_cols <- c(
  "All_Sig",                  # "Heritable" / "Non-heritable"
  "Abundance",                # mean relative abundance (0-1 or %; this script assumes fraction 0-1)
  "Heritability",             # %
  "Genome_size",              # AA count (or nt -> label accordingly)
  "Number_of_genes",
  "Mean_gene_size",           # AA
  "Degree","Closeness","Eigenvector","Betweenness",
  "correlations_beta","correlations_alpha","correlations_glucose",
  "Genus","Family","Phylum","Enterotype","Cluster",
  "Maternal","Cohousing","Variance","Hubs"
)

if (is.null(opt$input_tsv) || !file.exists(opt$input_tsv)) {
  stop("Please provide --input_tsv with a valid TSV file.")
}

df <- suppressMessages(readr::read_tsv(opt$input_tsv, guess_max = 1e5))

missing <- setdiff(required_cols, colnames(df))
if (length(missing)) {
  stop(paste0(
    "Missing required columns in input TSV: ",
    paste(missing, collapse=", "),
    "\nRequired set:\n", paste(required_cols, collapse=", ")
  ))
}

## ===================== Pre-processing ======================== ##
# coerce numeric columns safely
num_cols <- c("Abundance","Heritability","Genome_size","Number_of_genes","Mean_gene_size",
              "Degree","Closeness","Eigenvector","Betweenness",
              "correlations_beta","correlations_alpha","correlations_glucose",
              "Maternal","Cohousing","Variance")
df <- df %>%
  mutate(across(all_of(num_cols), ~ suppressWarnings(as.numeric(.x))))

# normalise All_Sig spelling + ordering
df <- df %>%
  mutate(All_Sig = case_when(
    str_to_lower(All_Sig) %in% c("heritable","significant","sig","yes") ~ "Heritable",
    TRUE ~ "Non-heritable"
  )) %>%
  mutate(All_Sig = factor(All_Sig, levels=c("Non-heritable","Heritable")))

# Clean labels but keep originals if needed
df <- df %>%
  mutate(
    Genus_clean  = str_remove(Genus, "^g__"),
    Family_clean = str_remove(Family,"^f__")
  )

subset_id <- opt$subset_id

## ===================== Plot helpers ========================== ##
theme_base <- theme_bw() +
  theme(
    legend.position="bottom",
    axis.line = element_line(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 15),
    plot.title = element_text(hjust = 0.5)
  )

lab_xy <- function(x="", y="") list(xlab(x), ylab(y))

# auto y-position above max, handling flat distributions
auto_ypos <- function(y, pad=0.07) {
  y <- y[is.finite(y)]
  if (!length(y)) return(0)
  rng <- range(y, na.rm=TRUE)
  max(rng[2] + pad * diff(rng), rng[2] + max(1e-6, abs(rng[2])*pad))
}

# Boxplot with Wilcoxon + annotation
box_wilcox <- function(data, x, y, x_order=NULL, ylab="", title=NULL,
                       fill_map=NULL, facet=NULL, y_limits=NULL, compare_pairs=NULL) {
  dd <- data
  if (!is.null(x_order)) dd[[x]] <- factor(dd[[x]], levels=x_order)
  p <- ggplot(dd, aes(x=.data[[x]], y=.data[[y]], fill=.data[[x]])) +
    geom_boxplot(outlier.shape = NA) +
    theme_base + lab_xy("", y) + ggtitle(title %||% paste0("Subset: ", subset_id)) +
    scale_x_discrete(guide = guide_axis(angle = 60))
  if (!is.null(fill_map)) {
    p <- p + scale_fill_manual(values = fill_map, guide = "none")
  } else {
    p <- p + guides(fill = "none")
  }
  if (!is.null(y_limits)) p <- p + coord_cartesian(ylim = y_limits)
  
  # comparisons
  if (is.null(compare_pairs)) {
    # default: first two levels if exist
    levs <- levels(dd[[x]])
    if (length(levs) >= 2) compare_pairs <- list(levs[1:2])
  }
  
  if (is.null(facet)) {
    y_pos <- auto_ypos(dd[[y]])
    p <- p + geom_signif(comparisons = compare_pairs,
                         test = "wilcox.test",
                         y_position = y_pos,
                         map_signif_level = TRUE)
  } else {
    # facet-wise annotations (manual=TRUE)
    ann <- dd %>%
      group_by(.data[[facet]]) %>%
      summarise(y_position = auto_ypos(.data[[y]]), .groups="drop")
    # one y_position per facet, one annotation per comparison
    if (!is.null(compare_pairs)) {
      # replicate rows per comparison
      ann <- tidyr::crossing(ann, tibble(start=sapply(compare_pairs, `[`, 1),
                                         end  =sapply(compare_pairs, `[`, 2)))
    } else {
      ann <- ann %>% mutate(start=levels(dd[[x]])[1], end=levels(dd[[x]])[2])
    }
    p <- p +
      geom_signif(
        data = ann,
        aes(xmin=start, xmax=end, y_position=y_position),
        manual = TRUE, inherit.aes = FALSE, map_signif_level = TRUE
      ) +
      facet_wrap(reformulate(facet))
  }
  
  return(p)
}

# Scatter with rho + P
scatter_cor <- function(data, x, y, method="pearson", xlab="", ylab="", title=NULL, xlim=NULL, ylim=NULL) {
  dd <- data %>% filter(is.finite(.data[[x]]), is.finite(.data[[y]]))
  if (!nrow(dd)) return(ggplot() + ggtitle(paste(title,"(no data)")) + theme_void())
  
  ct <- suppressWarnings(cor.test(dd[[x]], dd[[y]], method=method))
  r_value <- unname(ct$estimate)
  p_value <- ct$p.value
  p_lab <- if (p_value < 0.01) "P < 0.01" else if (p_value < 0.05) "P < 0.05" else "P > 0.05"
  
  ggplot(dd, aes(x=.data[[x]], y=.data[[y]])) +
    geom_point(size=1, alpha=.8) +
    theme_base + lab_xy(xlab, ylab) + ggtitle(title) +
    {if(!is.null(xlim)) coord_cartesian(xlim=xlim) else NULL} +
    {if(!is.null(ylim)) coord_cartesian(ylim=ylim) else NULL} +
    annotate("label", x=Inf, y=Inf,
             label=paste0("rho = ", round(r_value,3), "\n", p_lab),
             hjust=1.1, vjust=1.2, size=3)
}

`%||%` <- function(a,b) if (!is.null(a)) a else b

## ===================== Open PDF ============================== ##
pdf(opt$outpdf, width=8.5, height=6.3)

## ===================== (1) SIG vs NON-SIG ==================== ##
# Abundance vs Heritability (scatter + box)
print(scatter_cor(df, "Abundance", "Heritability",
                  method="pearson",
                  xlab="Mean relative abundance (fraction)",
                  ylab="Species heritability (%)",
                  title=NULL))
print(
  box_wilcox(df, x="All_Sig", y="Abundance",
             x_order=c("Non-heritable","Heritable"),
             ylab="Mean relative abundance (fraction)",
             fill_map=c("Non-heritable"="lightblue","Heritable"="lightcoral"),
             y_limits=c(0, quantile(df$Abundance, .99, na.rm=TRUE)*1.1))
)

# Genome size
print(
  box_wilcox(df, x="All_Sig", y="Genome_size",
             x_order=c("Non-heritable","Heritable"),
             ylab="Genome size (AAs)",
             fill_map=c("Non-heritable"="lightblue","Heritable"="lightcoral"))
)

# Genome size by Family (labels cleaned)
print(
  ggplot(df %>% mutate(Family_clean=fct_reorder(Family_clean, Genome_size, .fun=median, na.rm=TRUE)),
         aes(x=Family_clean, y=Genome_size)) +
    geom_boxplot() + theme_base + lab_xy("", "Genome size (AAs)") +
    scale_x_discrete(guide=guide_axis(angle=60))
)

# Number of genes
print(
  box_wilcox(df, x="All_Sig", y="Number_of_genes",
             x_order=c("Non-heritable","Heritable"),
             ylab="Number of genes",
             fill_map=c("Non-heritable"="lightblue","Heritable"="lightcoral"))
)

# Mean gene size
print(
  box_wilcox(df, x="All_Sig", y="Mean_gene_size",
             x_order=c("Non-heritable","Heritable"),
             ylab="Mean gene size (AAs)",
             fill_map=c("Non-heritable"="lightblue","Heritable"="lightcoral"))
)

# Centralities vs All_Sig
for (met in c("Degree","Closeness","Eigenvector","Betweenness")) {
  ylab_txt <- switch(met,
                     Degree="Degree centrality",
                     Closeness="Closeness centrality",
                     Eigenvector="Eigenvector centrality",
                     Betweenness="Betweenness centrality")
  print(
    box_wilcox(df, x="All_Sig", y=met,
               x_order=c("Non-heritable","Heritable"),
               ylab=ylab_txt,
               fill_map=c("Non-heritable"="lightblue","Heritable"="lightcoral"))
  )
}

# Correlations with beta/alpha
print(
  box_wilcox(df, x="All_Sig", y="correlations_beta",
             x_order=c("Non-heritable","Heritable"),
             ylab="Pearson correlation with beta-diversity (PCoA1)",
             y_limits=c(-1,1),
             fill_map=c("Non-heritable"="lightblue","Heritable"="lightcoral"))
)
print(
  box_wilcox(df, x="All_Sig", y="correlations_alpha",
             x_order=c("Non-heritable","Heritable"),
             ylab="Pearson correlation with alpha-diversity",
             y_limits=c(-1,1),
             fill_map=c("Non-heritable"="lightblue","Heritable"="lightcoral"))
)

## = (2) HIGHLY HERITABLE GENERA: Bacteroides vs Prevotella ===== ##
bp <- df %>%
  filter(Genus %in% c("g__Bacteroides","g__Prevotella")) %>%
  mutate(Genus_clean=factor(Genus_clean, levels=c("Bacteroides","Prevotella")))

bp_fill <- c("Bacteroides"="lightblue","Prevotella"="lightgreen")

for (met in c("Genome_size","Number_of_genes","Heritability","Maternal","Cohousing","Variance",
              "correlations_beta","correlations_alpha","correlations_glucose")) {
  ylab_txt <- switch(met,
                     Genome_size="Genome size (AAs)",
                     Number_of_genes="Number of genes",
                     Heritability="Heritability (%)",
                     Maternal="Maternal effects (%)",
                     Cohousing="Cohousing effects (%)",
                     Variance="Variance",
                     correlations_beta="Pearson correlation with beta-diversity (PCoA1)",
                     correlations_alpha="Pearson correlation with alpha-diversity",
                     correlations_glucose="Pearson correlation with glucose levels")
  ylim <- if (grepl("correlations_", met)) c(-1,1) else NULL
  print(
    box_wilcox(bp, x="Genus_clean", y=met,
               x_order=c("Bacteroides","Prevotella"),
               ylab=ylab_txt, y_limits=ylim, fill_map=bp_fill)
  )
}

# Within each genus: Non-heritable vs Heritable, faceted by genus for select metrics
for (met in c("Variance","Maternal","Cohousing","Genome_size","Number_of_genes","Mean_gene_size",
              "Degree","Closeness","Betweenness","Eigenvector",
              "correlations_beta","correlations_alpha","correlations_glucose")) {
  ylab_txt <- switch(met,
                     Genome_size="Genome size (AAs)",
                     Number_of_genes="Number of genes",
                     Mean_gene_size="Mean gene size (AAs)",
                     Degree="Degree centrality",
                     Closeness="Closeness centrality",
                     Betweenness="Betweenness centrality",
                     Eigenvector="Eigenvector centrality",
                     Maternal="Maternal effects (%)",
                     Cohousing="Cohousing effects (%)",
                     Variance="Variance",
                     correlations_beta="Pearson correlation with beta-diversity (PCoA1)",
                     correlations_alpha="Pearson correlation with alpha-diversity",
                     correlations_glucose="Pearson correlation with glucose levels")
  ylim <- if (grepl("correlations_", met)) c(-1,1) else NULL
  print(
    box_wilcox(bp, x="All_Sig", y=met,
               x_order=c("Non-heritable","Heritable"),
               ylab=ylab_txt, y_limits=ylim,
               facet="Genus_clean",
               fill_map=c("Non-heritable"="#1f77b4","Heritable"="#aec7e8"),
               compare_pairs=list(c("Non-heritable","Heritable")))
  )
}

## === (+ Paramuribaculum) Three-way genus comparisons ========= ##
bpp <- df %>%
  filter(Genus %in% c("g__Paramuribaculum","g__Bacteroides","g__Prevotella")) %>%
  mutate(Genus_clean=factor(Genus_clean, levels=c("Paramuribaculum","Bacteroides","Prevotella")))
bpp_fill <- c("Paramuribaculum"="orange","Bacteroides"="lightblue","Prevotella"="lightgreen")
pairs3 <- list(c("Paramuribaculum","Bacteroides"), c("Bacteroides","Prevotella"), c("Paramuribaculum","Prevotella"))

for (met in c("Heritability","Abundance","Genome_size","Number_of_genes","Mean_gene_size","correlations_beta")) {
  ylab_txt <- switch(met,
                     Abundance="Relative abundance (fraction)",
                     Genome_size="Genome size (AAs)",
                     Number_of_genes="Number of genes",
                     Mean_gene_size="Mean gene size (AAs)",
                     Heritability="Species heritability (%)",
                     correlations_beta="Pearson correlation with beta-diversity (PCoA1)")
  ylim <- if (grepl("correlations_", met)) c(-1,1) else NULL
  print(
    box_wilcox(bpp, x="Genus_clean", y=met,
               x_order=levels(bpp$Genus_clean),
               ylab=ylab_txt, y_limits=ylim,
               fill_map=bpp_fill, compare_pairs=pairs3)
  )
}

## ===================== (3) HUBS ============================== ##
hub_order <- c("Non-hub","Cluster Hub","Global Hub","Global and Cluster Hub")
for (met in c("Genome_size","Number_of_genes","Mean_gene_size")) {
  ylab_txt <- switch(met,
                     Genome_size="Genome size (AAs)",
                     Number_of_genes="Number of genes",
                     Mean_gene_size="Mean gene size (AAs)")
  print(
    ggplot(df, aes(x=factor(Hubs, levels=hub_order), y=.data[[met]])) +
      geom_boxplot() + theme_base + lab_xy("", ylab_txt) +
      scale_x_discrete(guide=guide_axis(angle=60)) +
      ggtitle(paste0("Subset: ", subset_id)) +
      geom_signif(comparisons=list(c("Non-hub","Cluster Hub"),
                                   c("Cluster Hub","Global Hub"),
                                   c("Global Hub","Global and Cluster Hub")),
                  map_signif_level=TRUE)
  )
}

## ===================== (4) TAXA slices ======================= ##
# Mean gene size by Phylum / Family / Genus (Bacteroidota/Bacteroidaceae)
print(
  ggplot(df %>% mutate(Phylum=fct_reorder(Phylum, Mean_gene_size, .fun=median, na.rm=TRUE)),
         aes(x=Phylum, y=Mean_gene_size)) +
    geom_boxplot() + theme_base + lab_xy("", "Mean gene size (AAs)") +
    ggtitle(paste0("Subset: ", subset_id)) +
    scale_x_discrete(guide=guide_axis(angle=60))
)
print(
  ggplot(df %>% filter(Phylum=="p__Bacteroidota") %>%
           mutate(Family=fct_reorder(Family, Mean_gene_size, .fun=median, na.rm=TRUE)),
         aes(x=Family, y=Mean_gene_size)) +
    geom_boxplot() + theme_base + lab_xy("", "Mean gene size (AAs)") +
    ggtitle(paste0("Subset: ", subset_id)) +
    scale_x_discrete(guide=guide_axis(angle=60))
)
print(
  ggplot(df %>% filter(Family=="f__Bacteroidaceae") %>%
           mutate(Genus=fct_reorder(Genus, Mean_gene_size, .fun=median, na.rm=TRUE)),
         aes(x=Genus, y=Mean_gene_size)) +
    geom_boxplot() + theme_base + lab_xy("", "Mean gene size (AAs)") +
    ggtitle(paste0("Subset: ", subset_id)) +
    scale_x_discrete(guide=guide_axis(angle=60))
)
print(
  ggplot(df %>% mutate(Cluster=fct_reorder(Cluster, Mean_gene_size, .fun=median, na.rm=TRUE)),
         aes(x=Cluster, y=Mean_gene_size)) +
    geom_boxplot() + theme_base + lab_xy("", "Mean gene size (AAs)") +
    ggtitle(paste0("Subset: ", subset_id)) +
    scale_x_discrete(guide=guide_axis(angle=60))
)

## ===== (5) CENTRALITY vs heritability + by All_Sig =========== ##
for (met in c("Degree","Closeness","Betweenness","Eigenvector")) {
  ylab_txt <- switch(met,
                     Degree="Species degree centrality",
                     Closeness="Species closeness centrality",
                     Betweenness="Species betweenness centrality",
                     Eigenvector="Species eigenvector centrality")
  print(
    scatter_cor(df, x=met, y="Heritability", method="spearman",
                xlab=ylab_txt, ylab="Species heritability (%)", title=NULL)
  )
  print(
    box_wilcox(df, x="All_Sig", y=met,
               x_order=c("Non-heritable","Heritable"),
               ylab=ylab_txt)
  )
}

## ===== (6) Correlations with beta/alpha overall/by groups ===== ##
print(scatter_cor(df, "correlations_beta", "Heritability",
                  method="pearson",
                  xlab="Correlation with beta-diversity (PCoA1)",
                  ylab="Species heritability (%)",
                  ylim=c(0, max(df$Heritability, na.rm=TRUE)*1.1)))
print(
  box_wilcox(df, x="All_Sig", y="correlations_beta",
             x_order=c("Non-heritable","Heritable"),
             ylab="Correlation with beta-diversity (PCoA1)",
             y_limits=c(-1,1))
)

print(scatter_cor(df, "correlations_alpha", "Heritability",
                  method="pearson",
                  xlab="Correlation with alpha-diversity",
                  ylab="Species heritability (%)",
                  ylim=c(0, max(df$Heritability, na.rm=TRUE)*1.1)))
print(
  box_wilcox(df, x="All_Sig", y="correlations_alpha",
             x_order=c("Non-heritable","Heritable"),
             ylab="Correlation with alpha-diversity",
             y_limits=c(-1,1))
)

# beta/alpha inter-relations at species-level
print(scatter_cor(df, "correlations_beta", "correlations_alpha",
                  method="pearson",
                  xlab="Correlation with beta-diversity (PCoA1)",
                  ylab="Correlation with alpha-diversity",
                  xlim=c(-1,1), ylim=c(-1,1)))

# By Phylum / Family(Bacteroidota) / Genus(Bacteroidaceae) / Enterotype / Cluster
by_slices <- list(
  list(var="Phylum", y="correlations_beta", ylab="Correlation with beta-diversity (PCoA1)", ylim=c(-1,1),
       sig_pairs=list(c("p__Bacteroidota","p__Firmicutes_A"))),
  list(var="Phylum", y="correlations_alpha", ylab="Correlation with alpha-diversity", ylim=c(-1,1),
       sig_pairs=list(c("p__Bacteroidota","p__Firmicutes_A")))
)

for (sl in by_slices) {
  v <- sl$var; yv <- sl$y
  ddd <- df
  if (v=="Family") ddd <- df %>% filter(Phylum=="p__Bacteroidota")
  if (v=="Genus")  ddd <- df %>% filter(Family=="f__Bacteroidaceae")
  print(
    ggplot(ddd %>% mutate(!!v := fct_reorder(.data[[v]], .data[[yv]], .fun=median, na.rm=TRUE)),
           aes(x=.data[[v]], y=.data[[yv]])) +
      geom_boxplot() + theme_base + lab_xy("", sl$ylab) +
      scale_x_discrete(guide=guide_axis(angle=60)) +
      {if (!is.null(sl$ylim)) coord_cartesian(ylim=sl$ylim) else NULL} +
      {if (!is.null(sl$sig_pairs)) geom_signif(comparisons=sl$sig_pairs, map_signif_level=TRUE) else NULL}
  )
}

for (v in c("Family","Genus","Cluster")) {
  yv <- "correlations_beta"; ylab_txt <- "Correlation with beta-diversity (PCoA1)"
  ddd <- df
  if (v=="Family") ddd <- df %>% filter(Phylum=="p__Bacteroidota")
  if (v=="Genus")  ddd <- df %>% filter(Family=="f__Bacteroidaceae")
  print(
    ggplot(ddd %>% mutate(!!v := fct_reorder(.data[[v]], .data[[yv]], .fun=median, na.rm=TRUE)),
           aes(x=.data[[v]], y=.data[[yv]])) +
      geom_boxplot() + theme_base + lab_xy("", ylab_txt) +
      scale_x_discrete(guide=guide_axis(angle=60)) +
      coord_cartesian(ylim=c(-1,1))
  )
}

for (v in c("Family","Genus","Enterotype","Cluster")) {
  yv <- "correlations_alpha"; ylab_txt <- "Correlation with alpha-diversity"
  ddd <- df
  if (v=="Family") ddd <- df %>% filter(Phylum=="p__Bacteroidota")
  if (v=="Genus")  ddd <- df %>% filter(Family=="f__Bacteroidaceae")
  print(
    ggplot(ddd %>% mutate(!!v := fct_reorder(.data[[v]], .data[[yv]], .fun=median, na.rm=TRUE)),
           aes(x=.data[[v]], y=.data[[yv]])) +
      geom_boxplot() + theme_base + lab_xy("", ylab_txt) +
      scale_x_discrete(guide=guide_axis(angle=60)) +
      coord_cartesian(ylim=c(-1,1)) +
      {if (v=="Enterotype") geom_signif(comparisons=list(c("Enterotype 1","Enterotype 2")),
                                        map_signif_level=TRUE) else NULL}
  )
}

## ===== (7) Within-genus cluster comparisons ================== ##
# Prevotella: cluster__10 vs cluster__11
prev <- df %>% filter(Genus=="g__Prevotella", Cluster %in% c("cluster__10","cluster__11")) %>%
  mutate(Cluster=factor(Cluster, levels=c("cluster__10","cluster__11")))
for (met in c("Heritability","Abundance","Prevalence","Genome_size","Number_of_genes",
              "Mean_gene_size","correlations_beta","correlations_alpha")) {
  ylab_txt <- switch(met,
                     Abundance="Mean relative abundance (fraction)",
                     Prevalence="Prevalence (%)",
                     Genome_size="Genome size (AAs)",
                     Number_of_genes="Number of genes",
                     Mean_gene_size="Mean gene size (AAs)",
                     Heritability="Species heritability (%)",
                     correlations_beta="Correlations with beta diversity (PCoA1)",
                     correlations_alpha="Correlations with alpha diversity (Rao’s)")
  ylim <- if (grepl("correlations_", met)) c(-1,1) else NULL
  print(
    box_wilcox(prev, x="Cluster", y=met,
               x_order=c("cluster__10","cluster__11"),
               ylab=ylab_txt, y_limits=ylim)
  )
}

# Bacteroides: cluster__10 vs cluster__3  (as per your code)
bact <- df %>% filter(Genus=="g__Bacteroides", Cluster %in% c("cluster__10","cluster__3")) %>%
  mutate(Cluster=factor(Cluster, levels=c("cluster__10","cluster__3")))
for (met in c("Abundance","Prevalence","Genome_size","Number_of_genes",
              "Mean_gene_size","correlations_beta","correlations_alpha")) {
  ylab_txt <- switch(met,
                     Abundance="Mean relative abundance (fraction)",
                     Prevalence="Prevalence (%)",
                     Genome_size="Genome size (AAs)",
                     Number_of_genes="Number of genes",
                     Mean_gene_size="Mean gene size (AAs)",
                     correlations_beta="Correlations with beta diversity (PCoA1)",
                     correlations_alpha="Correlations with alpha diversity (Rao’s)")
  ylim <- if (grepl("correlations_", met)) c(-1,1) else NULL
  print(
    box_wilcox(bact, x="Cluster", y=met,
               x_order=c("cluster__10","cluster__3"),
               ylab=ylab_txt, y_limits=ylim)
  )
}

## ===================== Close PDF ============================= ##
dev.off()

cat("\n✓ Plots written to: ", normalizePath(opt$outpdf, mustWork = FALSE), "\n", sep="")
cat("Subset label:", subset_id, "\n")

# End of script
