#!/usr/bin/env Rscript

# Libraries
suppressMessages(library(argparse))
suppressMessages(library(xfun))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(NetCoMi))
suppressMessages(library(zCompositions))
suppressMessages(library(limma))

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

# -----------------------------------------------------------------------------


# Create parser object
parser <- ArgumentParser()
print("Parsing arguments")

# Define desired outputs:
## Modules
parser$add_argument("-mod_comp", "--module_comparisons", type="character", help="Inclusion of module: Comparisons",default="NO")
## Matrix processing
parser$add_argument("-residues", "--input_residues", type="character", help="File with input matrix: residues")
parser$add_argument("-covariates", "--list_of_covariates", type="character", help="File with subsets and covariates info")
parser$add_argument("-cov_types", "--list_of_covariate_types", type="character", help="File with information of what type of variable the covariate is")
parser$add_argument("-ent", "--enterotypes", type="character", help="File with enterotypes info")
parser$add_argument("-met", "--method", type="character", help="Define the method", default="Shallow")
## Taxonomic
parser$add_argument("-tax", "--input_taxonomy", type="character", help="File with taxonomy info")
parser$add_argument("-ranks", "--taxonomic_ranks", type="character", help="Taxonomic ranks to be analysed",default="Phylum,Class,Order,Family,Genus,Species")
## Network
parser$add_argument("-measure", "--network_measure", type="character", help="Type of method for network measuring",default="spieceasi")
parser$add_argument("-cluster", "--cluster_method", type="character", help="Type of method for cluster definition",default="cluster_fast_greedy")
parser$add_argument("-hub", "--hub_definition", type="character", help="Measurement to be used for hub definition",default="degree")

# Get command line options, if help option encountered - print help and exit:
args <- parser$parse_args()
## Modules
module_comparisons<-args$module_comparisons
if (module_comparisons != "YES" & module_comparisons != "NO") {
  stop("Option for module alpha diversity should be either 'YES' or 'NO' in capital letters.", "\n")
}
## Matrix processing
input_residues<-args$input_residues
list_of_covariates<-args$list_of_covariates
list_of_covariate_types<-args$list_of_covariate_types
enterotypes<-args$enterotypes
sample_identifier<-args$sample_identifier
method=args$method
## Taxonomic
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
## Network
network_measure<-args$network_measure
cluster_method<-args$cluster_method
hub_definition<-args$hub_definition

# -----------------------------------------------------------------------------


### Step 1 ####
print("Step 1: Read input files and arguments")

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
    print("All elements in the subsets are present in the subset 'ALL'.")
  } else {
    print("Warning: Covariates ", paste(missing_covs, collapse = ", "), " are missing in the subset 'ALL' but used on other subsets.", "\n")
  }
  if (list_of_covariate_types != "") {
    covs_types<-read_tab_file(list_of_covariate_types, TRUE)
    if (ncol(covs_types)==1){
      covs_types<-read_csv_file(list_of_covariate_types)
    }
    if (nrow(covs_types)==0) {
      stop("No covariate was provided to the covariate types file. Please, include at least one known covariate on the input file covariates.classification.txt")
    }
    for (i in 1:nrow(covs_types)) {
      if (covs_types$classification[i] == "") {
        stop(paste0("Subset ",covs_types$covariate[i]," has no classification assigned. Please, inform on the input file covariates.classification.txt if the covariate is 'categorical' or 'continuous'.")) 
      }
    }
    for (i in 1:length(covnames)) {
      if (length(which(covs_types$covariate==covnames[i])) == 0) {
        stop(paste0("No classification was assigned to covariate: ",covnames[i],". Please, inform on the input file covariates.classification.txt if this covariate is 'categorical' or 'continuous'"))
      }
    }
    for (subset_id in subset_ids) {
      covariates=get(paste("covariates",subset_id,sep="_"))
      covs=c()
      for (i in 1:length(covariates)){
        cat=covs_types$classification[match(covariates[i],covs_types$covariate)]
        if (cat=="categorical"){
          covs=c(covs,covariates[i])
        }
      }
      assign(paste('covariates',subset_id,sep='_'), covs)
      assign(paste('all_covariates',subset_id,sep='_'), c(covariates))
    }
  }  else {
    stop("No covariate types file provided")
  }
} else {
  stop("No covariates file provided")
}

# Read enterotypes
print("Reading enterotypes info")
if (enterotypes != ""){
  load(enterotypes)
  enterotype_samples=enterotypes_samples_objs[[length(enterotypes_samples_objs)]]
  #enterotype_samples$Enterotype=gsub("enterotype__","Enterotype ",enterotype_samples$Enterotype)
  ents=unique(enterotype_samples$Enterotype)
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
if (! "Species" %in% ranks){
  stop("To define networks, please define it on the Species level.", "\n")
}

# Read residues
if (input_residues != "") {
  print("Reading filtered matrices per subset")
  load(input_residues)
  if (module_comparisons=="YES"){
    permutations <- expand.grid(subset_ids,ranks,methods)
    colnames(permutations)=c("subset_ids","ranks","methods")
    concatenated <- paste(permutations$ranks, permutations$subset_ids, permutations$methods,sep = "_")
  } else {
    permutations <- expand.grid(subset_ids,ranks)
    colnames(permutations)=c("subset_ids","ranks")
    concatenated <- paste(permutations$ranks, permutations$subset_ids, sep = "_")
  }
  for (i in 1:length(concatenated)){
    residuals_qned_counts=residuals_qned_counts_objs[[i]]
    assign(paste('residuals_qned_counts',concatenated[i],sep='_'), residuals_qned_counts)
  }
} else {
  stop("No matrices provided")
}

# Separate samples by enterotypes
print("Separating samples by enterotypes")
for (ent in ents){
  samples=enterotype_samples$sample[enterotype_samples$Enterotype==ent]
  matrix=get(paste('residuals_qned_counts',"Species","ALL",sep="_"))
  matrix=matrix[,colnames(matrix) %in% samples]
  assign(paste('residuals_qned_counts',"Species",ent,sep="_"),matrix)
  print(ent)
}

### End of Step 1 ###

### Step 2: Network estimation ###
print("Step 2: Network estimation")

clusters_objs=list()
stats_clusters_objs=list()
hubs_objs=list()
stats_hubs_objs=list()
centralities_objs=list()
pdf('networks_enterotypes.pdf')
if ("Species" %in% ranks){
  for (ent in ents){
    # Load the matrix
    matrix=t(get(paste('residuals_qned_counts',"Species",ent,sep="_")))
    # Define the network
    net <- netConstruct(data = matrix, 
                        filtTax = "none",
                        filtSamp = "none",
                        measure = network_measure,
                        measurePar = list(nlambda = 10),
                        normMethod = "none", 
                        zeroMethod = "none",
                        sparsMethod = "none", 
                        dissFunc = "signed",
                        verbose = 2,
                        seed = 123456)
    props <- netAnalyze(net, 
                        centrLCC = FALSE,
                        avDissIgnoreInf = TRUE,
                        sPathNorm = FALSE,
                        clustMethod = cluster_method,
                        hubPar = c(hub_definition),
                        hubQuant = 0.9,
                        lnormFit = TRUE,
                        normDeg = FALSE,
                        normBetw = FALSE,
                        normClose = FALSE,
                        normEigen = FALSE)
    net_summary=summary(props)
    nclust <- as.numeric(names(table(props$clustering$clust1)))
    col <- c(topo.colors(nclust), rainbow(length(nclust)))
    # Plot the network
    myplot=plot(props, 
                main = "",
                sameLayout = TRUE, 
                layoutGroup = "union", 
                colorVec = col,
                borderCol = "gray40", 
                nodeSize = "degree", 
                cexNodes = 0.9, 
                nodeSizeSpread = 3, 
                edgeTranspLow = 80, 
                edgeTranspHigh = 50,
                showTitle = TRUE, 
                cexTitle = 2.8,
                mar = c(1,1,3,1), 
                repulsion = 0.9, 
                labels = FALSE, 
                rmSingles = "inboth",
                nodeFilter = "clustMin", 
                nodeFilterPar = 10, 
                nodeTransp = 50, 
                hubTransp = 30)
    print(myplot)
    title(main = paste0("Species network: ",gsub("enterotype__","Enterotype ",ent)), line = 0.5, cex.main = 2.0)
    # Define centralities
    centralities=as.data.frame(props$centralities$degree1)
    colnames(centralities)[1]="Degree"
    centralities$Betweeness=props$centralities$between1
    centralities$Closeness=props$centralities$close1
    centralities$Species=rownames(centralities)
    centralities$Genus=taxonomy$Genus[match(centralities$Species,taxonomy$Species)]
    centralities$Family=taxonomy$Family[match(centralities$Species,taxonomy$Species)]
    centralities$Order=taxonomy$Order[match(centralities$Species,taxonomy$Species)]
    centralities$Class=taxonomy$Class[match(centralities$Species,taxonomy$Species)]
    centralities$Phylum=taxonomy$Phylum[match(centralities$Species,taxonomy$Species)]
    # Define clusters
    clusters=as.data.frame(props$clustering$clust1)
    colnames(clusters)[1]="Cluster"
    clusters$Cluster=paste0("cluster__",clusters$Cluster)
    clusters$Species=rownames(clusters)
    clusters$Genus=taxonomy$Genus[match(clusters$Species,taxonomy$Species)]
    clusters$Family=taxonomy$Family[match(clusters$Species,taxonomy$Species)]
    clusters$Phylum=taxonomy$Phylum[match(clusters$Species,taxonomy$Species)]
    clusters$Degree=centralities$Degree[match(clusters$Species,centralities$Species)]
    clusters$Betweeness=centralities$Betweeness[match(clusters$Species,centralities$Species)]
    clusters$Closeness=centralities$Closeness[match(clusters$Species,centralities$Species)]
    stats_clusters=count(clusters,Cluster,Family)
    # Define hubs
    hubs=as.data.frame(net_summary$hubs)
    colnames(hubs)[1]="Species"
    hubs$Genus=taxonomy$Genus[match(hubs$Species,taxonomy$Species)]
    hubs$Family=taxonomy$Family[match(hubs$Species,taxonomy$Species)]
    hubs$Order=taxonomy$Order[match(hubs$Species,taxonomy$Species)]
    hubs$Class=taxonomy$Class[match(hubs$Species,taxonomy$Species)]
    hubs$Phylum=taxonomy$Phylum[match(hubs$Species,taxonomy$Species)]
    hubs$Cluster=clusters$Cluster[match(hubs$Species,clusters$Species)]
    stats_hubs=count(hubs,Genus)
    stats_hubs$Family=taxonomy$Family[match(stats_hubs$Genus,taxonomy$Genus)]
    stats_hubs$Order=taxonomy$Order[match(stats_hubs$Genus,taxonomy$Genus)]
    stats_hubs$Class=taxonomy$Class[match(stats_hubs$Genus,taxonomy$Genus)]
    stats_hubs$Phylum=taxonomy$Phylum[match(stats_hubs$Genus,taxonomy$Genus)]
    # Define objects to store the data
    clusters_obj <- paste('clusters',subset_id,sep='_')
    assign(clusters_obj, clusters)
    clusters_objs<- append(clusters_objs, list(clusters_obj=clusters))
    stats_clusters_obj <- paste('stats_clusters',subset_id,sep='_')
    assign(stats_clusters_obj, stats_clusters)
    stats_clusters_objs<- append(stats_clusters_objs, list(stats_clusters_obj=stats_clusters))
    hubs_obj <- paste('hubs',subset_id,sep='_')
    assign(hubs_obj, hubs)
    hubs_objs<- append(hubs_objs, list(hubs_obj=hubs))
    stats_hubs_obj <- paste('stats_hubs',subset_id,sep='_')
    assign(stats_hubs_obj, stats_hubs)
    stats_hubs_objs<- append(stats_hubs_objs, list(stats_hubs_obj=stats_hubs))
    centralities_obj <- paste('centralities',subset_id,sep='_')
    assign(centralities_obj, centralities)
    centralities_objs<- append(centralities_objs, list(centralities_obj=centralities))
    print(ent)
  }
} else {
  stop("To define networks, please define it on the Species level.", "\n")
}
dev.off()

# Save results
save(clusters_objs,
     hubs_objs,
     stats_hubs_objs,
     centralities_objs,
     file ="networks_enterotypes.RData")

### End of Step 2 ###