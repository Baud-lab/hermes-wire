#!/usr/bin/env Rscript

# Packages
suppressMessages(library(qvalue))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(gap))
suppressMessages(library(argparse))
suppressMessages(library(xfun))
suppressMessages(library(ggsignif))
suppressMessages(library(ape))
suppressMessages(library(tidyr))

# -----------------------------------------------------------------------------
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

# 2. Assign covariates to subsets
assignCovariate <- function(x) {
  subset_id <- x[1]
  covnames <- str_split(x[2], pattern=",", simplify = TRUE) 
  assign(paste('covariates',subset_id,sep='_'), covnames, envir = .GlobalEnv)
}

# 3. Inverse rank normalization
invrank= function(row) {qnorm((rank(row,na.last="keep",ties.method = 'random')-0.5)/sum(!is.na(row)))}

# -----------------------------------------------------------------------------

# Create parser object
parser <- ArgumentParser()
print("Parsing arguments")

# Define desired outputs:
parser$add_argument("-ranks", "--taxonomic_ranks", type="character", help="Taxonomic ranks to be analysed",default="Phylum,Class,Order,Family,Genus,Species")
parser$add_argument("-tax", "--input_taxonomy", type="character", help="File with taxonomy info")
parser$add_argument("-meta", "--input_meta", type="character", help="File with metadata info")
parser$add_argument("-id", "--sample_identifier", type="character", help="Field on the metadata file that identifies the host",default="host_subject_id")
parser$add_argument("-covariates", "--list_of_covariates", type="character", help="File with subsets and covariates info")
parser$add_argument("-met", "--method", type="character", help="Define the method", default="Shallow")
parser$add_argument("-herit", "--heritability", type="character", help="File with heritability info")
parser$add_argument("-net", "--network", type="character", help="File with clusters info")

# Get command line options, if help option encountered - print help and exit:
args <- parser$parse_args()


## Arguments

heritability=args$heritability
network_dir=args$network
method=args$method
method=args$method
sample_identifier<-args$sample_identifier
if (sample_identifier == "") {
  stop("Please inform the sample identifier. It should have exactly the same name of the column where sample IDs are informd on your metadata file.", "\n")
}
input_meta<-args$input_meta
list_of_covariates<-args$list_of_covariates
input_taxonomy<-args$input_taxonomy
ranks<-args$taxonomic_ranks
allowed <- c("Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV")
string_vector <- strsplit(ranks, ",")[[1]]
if (all(string_vector %in% allowed)) {
  print("All ranks are valid")
} else {
  invalid <- string_vector[!string_vector %in% allowed]
  stop("Invalid taxonomic rank(s) found: ", paste(invalid, collapse = ", "),"\n","Ranks should be: Phylum, Class, Order, Family, Genus, Species or ASV" ,"\n")
}



# -----------------------------------------------------------------------------

### Step 1 ####
print("Step 1: Read input files and arguments")

# Load Heritability values
if (heritability != ""){
  load(heritability)
} else {
  stop("No heritability file provided")
}

# Extract metadata info
if (input_meta != "") {
  print("Extracting metadata info")
  metadata<-read_tab_file(input_meta, TRUE)
  if (ncol(metadata)==1){
    metadata<-read_csv_file(input_meta)
  }
  if (length(which(colnames(metadata)==sample_identifier)) == 0) {
    stop("The metadata file has no column with the same name of the sample identifier informed")
  }
} else {
  stop("No metadata file provided")
}

# Read covariates info
if (list_of_covariates != "") {
  print("Reading covariates and covariate types info")
  covs_data<-read_tab_file(list_of_covariates, TRUE)
  if (ncol(covs_data)==1){
    covs_data<-read_csv_file(list_of_covariates)
  }
  subset_ids <- covs_data$subset_id
  if (nrow(covs_data)==0) {
    stop("No subset was provided. If you only worked with one cohort/Study, please include subset 'ALL' on the input file covariates.txt and assign known covariates to it")
  }
  for (i in 1:nrow(covs_data)) {
    if (covs_data$covariates[i] == "") {
      stop(paste0("Subset ",covs_data$subset_id[i]," has no covariate assigned. Please inform at least one known covariate for this subset on the input file covariates.txt and remember that the name of the covariate should correspond with the name of one of the columns of your matadata file. If more than one covariate are assigned for this subset, remember to list them separated by comma with no spaces.")) 
    }
  }
  for (i in 1:nrow(covs_data)) {
    if (nrow(covs_data) > 1 & covs_data$subset_id[i] == "ALL" & !grepl("Study",covs_data$covariates[i])) {
      covs_data$covariates[i]=paste(covs_data$covariates[i],"Study",sep=",")
      print(paste0("Warning: Subset ",covs_data$subset_id[i]," had no covariate 'Study' assigned to it on the input file covariates.txt. Given that other subsets were listed, this covariate was automaticaly included.")) 
    }
  }
  apply(covs_data, 1,  assignCovariate) ## please check: I changed hard coded 1 for subset_id
  covnames <- unique(c(str_split(covs_data$covariates[covs_data$subset_id=="ALL"], pattern=",", simplify = TRUE)))
  total_covs=c()
  for (subset_id in subset_ids[-which(subset_ids=="ALL")]){
    covariate=get(paste("covariates",subset_id,sep="_"))
    total_covs=unique(c(total_covs,covariate))
  }
  missing_covs <- total_covs[!(total_covs %in% covnames)]
  if (length(missing_covs) == 0) {
    print("All elements in the subsets are present in the subset 'ALL'.\n")
  } else {
    print("Warning: Covariates ", paste(missing_covs, collapse = ", "), " are missing in the subset 'ALL' but used on other subsets.", "\n")
  }
} else {
  stop("No covariates file provided")
}

# Read taxonomy
print("Reading taxonomy info")
if (input_taxonomy != "") {
  taxonomy<-read_tab_file(input_taxonomy, FALSE)
  if (grepl("GTDB", input_taxonomy)) {
    colnames(taxonomy)=c("Genome","Taxonomy")
    taxonomy$Taxonomy<-gsub("; ",";",as.character(taxonomy$Taxonomy))
    taxonomy$Taxonomy<-gsub(" ","_",as.character(taxonomy$Taxonomy))
  } else {
    if (grepl("Greengenes2", input_taxonomy)){
      colnames(taxonomy) <- c("Genome","Taxonomy","Confidence","GTDB")
      taxonomy$Taxonomy<-gsub("\\s(?!.*\\s)", "_", as.character(taxonomy$Taxonomy), perl = TRUE)
      taxonomy$Taxonomy<-gsub(" ", "", as.character(taxonomy$Taxonomy), perl = TRUE)
      taxonomy$ASV=paste0("asv__",rownames(taxonomy))
    } else {
      stop("The taxonomy provided should be GTDB for shotgun data or Greengenes2 for 16S")
    }
  }
  taxonomy$Phylum<-sapply(strsplit(taxonomy$Taxonomy, ";"), "[",2)
  taxonomy$Class<-sapply(strsplit(taxonomy$Taxonomy, ";"), "[",3)
  taxonomy$Order<-sapply(strsplit(taxonomy$Taxonomy, ";"), "[",4)
  taxonomy$Family<-sapply(strsplit(taxonomy$Taxonomy, ";"), "[",5)
  taxonomy$Genus<-sapply(strsplit(taxonomy$Taxonomy, ";"), "[",6)
  taxonomy$Species<-sapply(strsplit(taxonomy$Taxonomy, ";"), "[",7)
} else {
  stop("No taxonomy file provided for the main dataset")
}

# Read taxonomic ranks
print("Reading taxonomic ranks")
ranks<-na.omit(c((sapply(strsplit(ranks, ","),"[",1)),
                 (sapply(strsplit(ranks, ","), "[",2)),
                 (sapply(strsplit(ranks, ","), "[",3)),
                 (sapply(strsplit(ranks, ","), "[",4)),
                 (sapply(strsplit(ranks, ","), "[",5)),
                 (sapply(strsplit(ranks, ","), "[",6)),
                 (sapply(strsplit(ranks, ","), "[",7))))


# Read clusters information
if (network_dir != "") {
  print("Reading clusters information")
  network_files <- data.frame(file=list.files(network_dir, pattern=".csv"))
  network_files$subset=gsub(".csv","",network_files$file)
  network_files$subset= sapply(strsplit(network_files$subset, "_"), "[",3)
  for (i in 1:length(network_files$file)){
    clusters=read_csv_file(paste(network_dir,network_files$file[i],sep="/"))
    colnames(clusters)[which(colnames(clusters)=="ASV")]="Genome"
    colnames(clusters)[which(colnames(clusters)=="HDBSCAN")]="Cluster"
    clusters$Cluster=gsub("-1","not.clustered",clusters$Cluster)
    clusters$Cluster=paste("cluster",clusters$Cluster,sep="__")
    if (method=="16S"){
      clusters$ASV=taxonomy$ASV[match(clusters$Genome,taxonomy$Genome)]
    }
    assign(paste('clusters',network_files$subset[i],sep='_'), clusters)
  }
} else {
  stop("No clusters info provided")
}

### End of Step 1 ###

### Step 2 ####
print("Step 2: Make plots")

pdf('heritability_clusters.pdf',bg='white')
for (subset_id in subset_ids){
  filtered_VCs=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),]
  clusters=get(paste("clusters",subset_id,sep="_"))
  filtered_VCs <- filtered_VCs %>% 
    mutate( All_Sig = case_when(
      Significance ==  "Significant (Bonferroni)" | Significance == "Significant (FDR)" ~ "Significant",
      TRUE ~ "Non-significant"))
  if (method=="16S"){
    if ("ASV" %in% ranks){
      print("Comparisons between significant and non-significant ASVs")
      filtered_VCs=filtered_VCs[grepl("ASV",filtered_VCs$rank),]
      filtered_VCs$Genome=taxonomy$ASV[match(filtered_VCs$trait,taxonomy$Genome)]
      filtered_VCs$Phylum=taxonomy$Phylum[match(filtered_VCs$trait,taxonomy$Genome)]
      filtered_VCs$Class=taxonomy$Class[match(filtered_VCs$trait,taxonomy$Genome)]
      filtered_VCs$Order=taxonomy$Order[match(filtered_VCs$trait,taxonomy$Genome)]
      filtered_VCs$Family=taxonomy$Family[match(filtered_VCs$trait,taxonomy$Genome)]
      filtered_VCs$Genus=taxonomy$Genus[match(filtered_VCs$trait,taxonomy$Genome)]
      rank="ASVs"
    }
  } else {
    if ("Species" %in% ranks){
      print("Comparisons between significant and non-significant Species")
      filtered_VCs=filtered_VCs[grepl("Species",filtered_VCs$rank),]
      filtered_VCs$Genome=taxonomy$Genome[match(filtered_VCs$trait,taxonomy$Species)]
      filtered_VCs$Phylum=taxonomy$Phylum[match(filtered_VCs$trait,taxonomy$Species)]
      filtered_VCs$Class=taxonomy$Class[match(filtered_VCs$trait,taxonomy$Species)]
      filtered_VCs$Order=taxonomy$Order[match(filtered_VCs$trait,taxonomy$Species)]
      filtered_VCs$Family=taxonomy$Family[match(filtered_VCs$trait,taxonomy$Species)]
      filtered_VCs$Genus=taxonomy$Genus[match(filtered_VCs$trait,taxonomy$Species)]
      rank="Species"
    }
  }
  # Add cluster and centrality values
  filtered_VCs$Degree_Centrality=clusters$Degree_Centrality[match(filtered_VCs$Genome,clusters$Genome)]
  filtered_VCs$Betweeness_Centrality=clusters$Betweeness_Centrality[match(filtered_VCs$Genome,clusters$Genome)]
  filtered_VCs$Closeness_Centrality=clusters$Closeness_Centrality[match(filtered_VCs$Genome,clusters$Genome)]
  filtered_VCs$PageRank=clusters$PageRank[match(filtered_VCs$Genome,clusters$Genome)]
  filtered_VCs$Cluster=clusters$Cluster[match(filtered_VCs$Genome,clusters$Genome)]
  
  # Calculate the number of significant (FDR) in each cluster
  stats_cluster_all_sig=count(filtered_VCs,All_Sig,Cluster)
  stats_cluster_sep_sig=count(filtered_VCs,Significance,Cluster)
  save(stats_cluster_all_sig, stats_cluster_sep_sig, file="statistics_per_cluster.RData")
  
  # Make the plots
  
  ## Heritability per Cluster
  myplot<-ggplot(filtered_VCs, aes(x=reorder(Cluster,100*var_Ad,FUN = median), y=100*var_Ad)) +
    geom_boxplot() +
    xlab("") +
    ylab("Heritability (%)") +
    ggtitle(paste0("Subset: ",subset_id)) +
    theme_bw() +
    theme(legend.position="bottom",
          axis.line = element_line(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          legend.text=element_text(size=10),
          plot.title = element_text(hjust = 0.5))+
    scale_x_discrete(guide = guide_axis(angle = 45))
  print(myplot)
  
  ## Plot the comparisons between the median heritability of all ASVs on the rank with the own rank heritability
  median_var_Ad <- aggregate(var_Ad ~ ., data = filtered_VCs[, c("var_Ad", "Cluster")], FUN = median)
  median_var_Ad = median_var_Ad[order(median_var_Ad$var_Ad, decreasing = T),]
  colnames(median_var_Ad)[1]="trait"
  if (length(median_var_Ad$trait)>10){
    median_var_Ad=median_var_Ad[1:10,]
  }
  order=median_var_Ad$trait
  filtered_VCs_cluster=all_VCs_full[grepl(subset_id,all_VCs_full$subset_id),]
  filtered_VCs_cluster=filtered_VCs_cluster[grepl("Cluster",filtered_VCs_cluster$rank),]
  filtered_VCs_cluster=filtered_VCs_cluster[filtered_VCs_cluster$trait %in% median_var_Ad$trait,]
  filtered_VCs_cluster$LBL[filtered_VCs_cluster$Significance == "Significant (Bonferroni)" ] <- "**"
  filtered_VCs_cluster$LBL[filtered_VCs_cluster$Significance == "Significant (FDR)"] <- "*"
  filtered_VCs_cluster$LBL[is.na(filtered_VCs_cluster$LBL)] <- ""
  filtered_VCs_cluster$heritability=paste(round(filtered_VCs_cluster$var_Ad*100,digit=2),filtered_VCs_cluster$LBL,sep="")
  filtered_VCs_cluster$median=median_var_Ad$var_Ad[match(filtered_VCs_cluster$trait,median_var_Ad$trait)]
  filtered_VCs_cluster$comp[filtered_VCs_cluster$var_Ad > filtered_VCs_cluster$median] <- "Cluster Heritability > Species Median"
  filtered_VCs_cluster$comp[filtered_VCs_cluster$var_Ad < filtered_VCs_cluster$median] <- "Cluster Heritability < Species Median"
  filtered_VCs_cluster$comp[is.na(filtered_VCs_cluster$comp)] <- ""
  plot_data=filtered_VCs[,c(which(colnames(filtered_VCs)=="var_Ad"),which(colnames(filtered_VCs)=="Cluster"))]
  plot_data=plot_data[plot_data$Cluster %in% median_var_Ad$trait,]
  # Plot
  myplot<-ggplot(plot_data, aes(x = factor(Cluster, levels = order), y = var_Ad*100)) +
    geom_boxplot() +
    geom_label(data = filtered_VCs_cluster,  
               aes(x = factor(trait, levels = order),
                   y = var_Ad*100, 
                   label = heritability, fill = comp),size = 3.5) +
    scale_y_continuous(breaks = seq(0, 100, by = 10)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_line(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position="bottom") +
    guides(fill = guide_legend(override.aes = list(color = NA)), 
           color = "none", 
           shape = "none") +
    scale_fill_manual(values = c("Cluster Heritability > Species Median" = "lightgreen", "Cluster Heritability < Species Median" = "lightcoral"),
                      name = "") +
    labs(x = "", y=paste0("Heritability (%) (Subset: ",subset_id))
  print(myplot)
  
  if (length(filtered_VCs[filtered_VCs$All_Sig=="Significant"])>0 ){
    ## Comparing all types of Significance
    order=c("Non-significant traits","Significant (FDR)","Significant (Bonferroni)")
    ### Degree Centrality
    myplot<-ggplot(filtered_VCs, aes(x = factor(Significance, levels = order), y=100*Degree_Centrality)) +
      geom_boxplot() +
      xlab("") +
      ylab("Degree Centrality (%)") +
      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
      theme_bw() +
      theme(legend.position="bottom",
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text=element_text(size=10),
            plot.title = element_text(hjust = 0.5))+
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      geom_signif(comparisons = list(c("Non-significant traits","Significant (FDR)"),c("Significant (FDR)","Significant (Bonferroni)")),
                  map_signif_level = TRUE)
    print(myplot)
    ### Betweeness Centrality
    myplot<-ggplot(filtered_VCs, aes(x = factor(Significance, levels = order), y=100*Betweeness_Centrality)) +
      geom_boxplot() +
      xlab("") +
      ylab("Betweeness Centrality (%)") +
      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
      theme_bw() +
      theme(legend.position="bottom",
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text=element_text(size=10),
            plot.title = element_text(hjust = 0.5))+
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      geom_signif(comparisons = list(c("Non-significant traits","Significant (FDR)"),c("Significant (FDR)","Significant (Bonferroni)")),
                  map_signif_level = TRUE)
    print(myplot)
    ### Closeness Centrality
    myplot<-ggplot(filtered_VCs, aes(x = factor(Significance, levels = order), y=100*Closeness_Centrality)) +
      geom_boxplot() +
      xlab("") +
      ylab("Closeness Centrality (%)") +
      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
      theme_bw() +
      theme(legend.position="bottom",
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text=element_text(size=10),
            plot.title = element_text(hjust = 0.5))+
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      geom_signif(comparisons = list(c("Non-significant traits","Significant (FDR)"),c("Significant (FDR)","Significant (Bonferroni)")),
                  map_signif_level = TRUE)
    print(myplot)
    ### Page Rank
    myplot<-ggplot(filtered_VCs, aes(x = factor(Significance, levels = order), y=PageRank)) +
      geom_boxplot() +
      xlab("") +
      ylab("Page Rank") +
      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
      theme_bw() +
      theme(legend.position="bottom",
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text=element_text(size=10),
            plot.title = element_text(hjust = 0.5))+
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      geom_signif(comparisons = list(c("Non-significant traits","Significant (FDR)"),c("Significant (FDR)","Significant (Bonferroni)")),
                  map_signif_level = TRUE)
    print(myplot)
    
    
    ## Comparing Significant (FDR) vs Non-Significant
    order=c("Non-significant","Significant")
    ### Degree Centrality
    myplot<-ggplot(filtered_VCs, aes(x = factor(All_Sig, levels = order), y=100*Degree_Centrality)) +
      geom_boxplot() +
      xlab("") +
      ylab("Degree Centrality (%)") +
      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
      theme_bw() +
      theme(legend.position="bottom",
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text=element_text(size=10),
            plot.title = element_text(hjust = 0.5))+
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      geom_signif(comparisons = list(c("Non-significant","Significant")),
                  map_signif_level = TRUE)
    print(myplot)
    ### Betweeness Centrality
    myplot<-ggplot(filtered_VCs, aes(x = factor(All_Sig, levels = order), y=100*Betweeness_Centrality)) +
      geom_boxplot() +
      xlab("") +
      ylab("Betweeness Centrality (%)") +
      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
      theme_bw() +
      theme(legend.position="bottom",
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text=element_text(size=10),
            plot.title = element_text(hjust = 0.5))+
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      geom_signif(comparisons = list(c("Non-significant","Significant")),
                  map_signif_level = TRUE)
    print(myplot)
    ### Closeness Centrality
    myplot<-ggplot(filtered_VCs, aes(x = factor(All_Sig, levels = order), y=100*Closeness_Centrality)) +
      geom_boxplot() +
      xlab("") +
      ylab("Closeness Centrality (%)") +
      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
      theme_bw() +
      theme(legend.position="bottom",
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text=element_text(size=10),
            plot.title = element_text(hjust = 0.5))+
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      geom_signif(comparisons = list(c("Non-significant","Significant")),
                  map_signif_level = TRUE)
    print(myplot)
    ### Page Rank
    myplot<-ggplot(filtered_VCs, aes(x = factor(All_Sig, levels = order), y=PageRank)) +
      geom_boxplot() +
      xlab("") +
      ylab("Page Rank") +
      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
      theme_bw() +
      theme(legend.position="bottom",
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text=element_text(size=10),
            plot.title = element_text(hjust = 0.5))+
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      geom_signif(comparisons = list(c("Non-significant","Significant")),
                  map_signif_level = TRUE)
    print(myplot)
    
    ## Taxonomy: Comparing the two clusters with a higher proportion of species with significant heritability
    top_2_clusters <- filtered_VCs %>%
      group_by(Cluster) %>%
      summarize(percentage_significant = mean(All_Sig == "Significant") * 100) %>%
      top_n(2, percentage_significant)
    filtered_cluster=filtered_VCs[filtered_VCs$Cluster %in% top_2_clusters$Cluster,]
    ## Comparing Significant (FDR) vs Non-Significant (Only top two clusters)
    order=c("Non-significant","Significant")
    ### Degree Centrality
    max_Degree_Centrality=max(filtered_cluster$Degree_Centrality)*100
    myplot<-ggplot(filtered_cluster, aes(x = factor(All_Sig, levels = order), y=100*Degree_Centrality)) +
      geom_boxplot() +
      xlab("") +
      ylab("Degree Centrality (%)") +
      ggtitle(paste0("Subset: ",subset_id)) +
      facet_wrap(~ Cluster, scales = "free") +
      ylim(0, max_Degree_Centrality) +
      theme_bw() +
      theme(legend.position="bottom",
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text=element_text(size=10),
            plot.title = element_text(hjust = 0.5))+
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      geom_signif(comparisons = list(c("Non-significant","Significant")),
                  map_signif_level = TRUE)
    print(myplot)
    ### Betweeness Centrality
    max_Betweeness_Centrality=max(filtered_cluster$Betweeness_Centrality)*100
    myplot<-ggplot(filtered_cluster, aes(x = factor(All_Sig, levels = order), y=100*Betweeness_Centrality)) +
      geom_boxplot() +
      xlab("") +
      ylab("Betweeness (%)") +
      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
      facet_wrap(~ Cluster, scales = "free") +
      ylim(0, max_Betweeness_Centrality) +
      theme_bw() +
      theme(legend.position="bottom",
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text=element_text(size=10),
            plot.title = element_text(hjust = 0.5))+
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      geom_signif(comparisons = list(c("Non-significant","Significant")),
                  map_signif_level = TRUE)
    print(myplot)
    ### Closeness Centrality
    max_Closeness_Centrality=max(filtered_cluster$Closeness_Centrality)*100
    myplot<-ggplot(filtered_cluster, aes(x = factor(All_Sig, levels = order), y=100*Closeness_Centrality)) +
      geom_boxplot() +
      xlab("") +
      ylab("Closeness Centrality (%)") +
      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
      facet_wrap(~ Cluster, scales = "free") +
      ylim(0, max_Closeness_Centrality) +
      theme_bw() +
      theme(legend.position="bottom",
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text=element_text(size=10),
            plot.title = element_text(hjust = 0.5))+
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      geom_signif(comparisons = list(c("Non-significant","Significant")),
                  map_signif_level = TRUE)
    print(myplot)
    
    
    ## Comparing Only Top Two Clusters
    ### Heritability
    myplot<-ggplot(filtered_cluster, aes(x=reorder(Cluster,100*var_Ad,FUN = median), y=100*var_Ad)) +
      geom_boxplot() +
      xlab("") +
      ylab("Heritability (%)") +
      ggtitle(paste0("Subset: ",subset_id)) +
      theme_bw() +
      theme(legend.position="bottom",
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text=element_text(size=10),
            plot.title = element_text(hjust = 0.5))+
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      geom_signif(comparisons = list(c(top_2_clusters$Cluster[1],top_2_clusters$Cluster[2])),
                  map_signif_level = TRUE)
    print(myplot)
    ### Degree Centrality
    myplot<-ggplot(filtered_cluster, aes(x=reorder(Cluster,100*Degree_Centrality,FUN = median), y=100*Degree_Centrality)) +
      geom_boxplot() +
      xlab("") +
      ylab("Degree Centrality (%)") +
      ggtitle(paste0("Subset: ",subset_id)) +
      theme_bw() +
      theme(legend.position="bottom",
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text=element_text(size=10),
            plot.title = element_text(hjust = 0.5))+
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      geom_signif(comparisons = list(c(top_2_clusters$Cluster[1],top_2_clusters$Cluster[2])),
                  map_signif_level = TRUE)
    print(myplot)
    ### Betweeness Centrality
    myplot<-ggplot(filtered_cluster, aes(x=reorder(Cluster,100*Betweeness_Centrality,FUN = median), y=100*Betweeness_Centrality)) +
      geom_boxplot() +
      xlab("") +
      ylab("Betweeness Centrality (%)") +
      ggtitle(paste0("Subset: ",subset_id)) +
      theme_bw() +
      theme(legend.position="bottom",
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text=element_text(size=10),
            plot.title = element_text(hjust = 0.5))+
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      geom_signif(comparisons = list(c(top_2_clusters$Cluster[1],top_2_clusters$Cluster[2])),
                  map_signif_level = TRUE)
    print(myplot)
    
    ### Closeness Centrality
    myplot<-ggplot(filtered_cluster, aes(x=reorder(Cluster,100*Closeness_Centrality,FUN = median), y=100*Closeness_Centrality)) +
      geom_boxplot() +
      xlab("") +
      ylab("Closeness Centrality (%)") +
      ggtitle(paste0("Subset: ",subset_id)) +
      theme_bw() +
      theme(legend.position="bottom",
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text=element_text(size=10),
            plot.title = element_text(hjust = 0.5))+
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      geom_signif(comparisons = list(c(top_2_clusters$Cluster[1],top_2_clusters$Cluster[2])),
                  map_signif_level = TRUE)
    print(myplot)
    
    ## Taxonomy: Comparing the two genera with a higher number of species with significant heritability
    top_2_genera <- filtered_VCs %>%
      filter(All_Sig == "Significant") %>%
      group_by(Genus) %>%
      summarize(num_significant = n()) %>%
      top_n(2, num_significant)
    filtered_genera=filtered_VCs[filtered_VCs$Genus %in% top_2_genera$Genus,]
    
    ## Comparing Significant (FDR) vs Non-Significant (Only top two Genuss)
    order=c("Non-significant","Significant")
    ### Degree Centrality
    max_Degree_Centrality=max(filtered_genera$Degree_Centrality)*100
    myplot<-ggplot(filtered_genera, aes(x = factor(All_Sig, levels = order), y=100*Degree_Centrality)) +
      geom_boxplot() +
      xlab("") +
      ylab("Degree Centrality (%)") +
      ggtitle(paste0("Subset: ",subset_id)) +
      facet_wrap(~ Genus, scales = "free") +
      ylim(0, max_Degree_Centrality) +
      theme_bw() +
      theme(legend.position="bottom",
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text=element_text(size=10),
            plot.title = element_text(hjust = 0.5))+
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      geom_signif(comparisons = list(c("Non-significant","Significant")),
                  map_signif_level = TRUE)
    print(myplot)
    ### Betweeness Centrality
    max_Betweeness_Centrality=max(filtered_genera$Betweeness_Centrality)*100
    myplot<-ggplot(filtered_genera, aes(x = factor(All_Sig, levels = order), y=100*Betweeness_Centrality)) +
      geom_boxplot() +
      xlab("") +
      ylab("Betweeness (%)") +
      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
      facet_wrap(~ Genus, scales = "free") +
      ylim(0, max_Betweeness_Centrality) +
      theme_bw() +
      theme(legend.position="bottom",
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text=element_text(size=10),
            plot.title = element_text(hjust = 0.5))+
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      geom_signif(comparisons = list(c("Non-significant","Significant")),
                  map_signif_level = TRUE)
    print(myplot)
    ### Closeness Centrality
    max_Closeness_Centrality=max(filtered_genera$Closeness_Centrality)*100
    myplot<-ggplot(filtered_genera, aes(x = factor(All_Sig, levels = order), y=100*Closeness_Centrality)) +
      geom_boxplot() +
      xlab("") +
      ylab("Closeness Centrality (%)") +
      ggtitle(paste0("All ",rank,"/Subset: ",subset_id)) +
      facet_wrap(~ Genus, scales = "free") +
      ylim(0, max_Closeness_Centrality) +
      theme_bw() +
      theme(legend.position="bottom",
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text=element_text(size=10),
            plot.title = element_text(hjust = 0.5))+
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      geom_signif(comparisons = list(c("Non-significant","Significant")),
                  map_signif_level = TRUE)
    print(myplot)
    
    
    ## Comparing Only Top Two Genera
    ### Heritability
    myplot<-ggplot(filtered_genera, aes(x=reorder(Genus,100*var_Ad,FUN = median), y=100*var_Ad)) +
      geom_boxplot() +
      xlab("") +
      ylab("Heritability (%)") +
      ggtitle(paste0("Subset: ",subset_id)) +
      theme_bw() +
      theme(legend.position="bottom",
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text=element_text(size=10),
            plot.title = element_text(hjust = 0.5))+
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      geom_signif(comparisons = list(c(top_2_genera$Genus[1],top_2_genera$Genus[2])),
                  map_signif_level = TRUE)
    print(myplot)
    ### Degree Centrality
    myplot<-ggplot(filtered_genera, aes(x=reorder(Genus,100*Degree_Centrality,FUN = median), y=100*Degree_Centrality)) +
      geom_boxplot() +
      xlab("") +
      ylab("Degree Centrality (%)") +
      ggtitle(paste0("Subset: ",subset_id)) +
      theme_bw() +
      theme(legend.position="bottom",
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text=element_text(size=10),
            plot.title = element_text(hjust = 0.5))+
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      geom_signif(comparisons = list(c(top_2_genera$Genus[1],top_2_genera$Genus[2])),
                  map_signif_level = TRUE)
    print(myplot)
    ### Betweeness Centrality
    myplot<-ggplot(filtered_genera, aes(x=reorder(Genus,100*Betweeness_Centrality,FUN = median), y=100*Betweeness_Centrality)) +
      geom_boxplot() +
      xlab("") +
      ylab("Betweeness Centrality (%)") +
      ggtitle(paste0("Subset: ",subset_id)) +
      theme_bw() +
      theme(legend.position="bottom",
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text=element_text(size=10),
            plot.title = element_text(hjust = 0.5))+
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      geom_signif(comparisons = list(c(top_2_genera$Genus[1],top_2_genera$Genus[2])),
                  map_signif_level = TRUE)
    print(myplot)
    
    ### Closeness Centrality
    myplot<-ggplot(filtered_genera, aes(x=reorder(Genus,100*Closeness_Centrality,FUN = median), y=100*Closeness_Centrality)) +
      geom_boxplot() +
      xlab("") +
      ylab("Closeness Centrality (%)") +
      ggtitle(paste0("Subset: ",subset_id)) +
      theme_bw() +
      theme(legend.position="bottom",
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text=element_text(size=10),
            plot.title = element_text(hjust = 0.5))+
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      geom_signif(comparisons = list(c(top_2_genera$Genus[1],top_2_genera$Genus[2])),
                  map_signif_level = TRUE)
    print(myplot)
  }
}
dev.off()

### End of Step 3 ###