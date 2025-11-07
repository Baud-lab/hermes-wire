#!/usr/bin/env Rscript
###############################################################################
## Network reconstruction & analysis pipeline                                ##
## (clean version – 22 Jun 2025)                                             ##
###############################################################################

# ── Libraries ────────────────────────────────────────────────────────────────
suppressMessages({
  library(argparse);      library(xfun);       library(dplyr)
  library(stringr);       library(NetCoMi);    library(zCompositions)
  library(limma);         library(parallel);   library(WGCNA)
  library(igraph);        library(purrr);      library(LaplacesDemon)
})

# ── 1. Helper functions ──────────────────────────────────────────────────────
read_gzipped <- function(f) if (xfun::file_ext(f) == "gz") gzfile(f) else f

read_tab_file <- function(f, header = FALSE)
  read.delim(read_gzipped(f), header = header, check.names = FALSE)

read_csv_file <- function(f)
  read.csv(read_gzipped(f), stringsAsFactors = FALSE)

check_network <- function(mat) {
  g <- graph_from_adjacency_matrix(mat, weighted = TRUE, mode = "undirected")
  data.frame(
    Density        = sprintf("%.1f%%", mean(mat != 0) * 100),
    NegativeEdges  = sprintf("%.1f%%", sum(mat < 0) / sum(mat != 0) * 100),
    GiantComponent = sprintf("%.1f%%",
                             max(components(g)$csize) / ncol(mat) * 100),
    MedianStrength = sprintf("|r| = %.2f", median(abs(mat[mat != 0])))
  )
}

assignCovariate <- function(x) {
  subset_id <- x[1]
  covnames  <- str_split(x[2], ",", simplify = TRUE)
  assign(paste0("covariates_", subset_id), covnames, envir = .GlobalEnv)
}

# ── Sub-graph plotting helpers ───────────────────────────────────────────────
plot_module <- function(props, mid, modules, node_colors) {
  keep <- modules$Traits[modules$Module == mid]
  myplot <- plot(
    props, rmLoops = TRUE, selNet = 1,
    nodeFilter = "names", nodeFilterPar = keep,
    sameLayout = TRUE, layoutGroup = "union", repulsion = 0.9,
    nodeColor  = "feature", featVecCol = node_colors,
    nodeSize   = "eigenvector", nodeSizeSpread = 3,
    borderCol  = "gray40", labels = TRUE, nodeTransp = 50,
    negDiffCol = TRUE, posCol = "#2E8B57", negCol = "#FF6347",
    edgeWidth  = 2, edgeFilter = "threshold", edgeFilterPar = 0.10,
    edgeTranspLow = 80, edgeTranspHigh = 50,
    highlightHubs = FALSE, hubLabelThresh = 0.5, hubTransp = 30,
    cexHubLabels = 1.2
  )
  title(main = paste0("Module: ", sub("module__", "", mid)))
  print(myplot)
}

plot_cluster <- function(props, cid, clusters_df, node_colors) {
  keep <- clusters_df$Traits[clusters_df$Cluster == cid]
  myplot <- plot(
    props, rmLoops = TRUE, selNet = 1,
    nodeFilter = "names", nodeFilterPar = keep,
    sameLayout = TRUE, layoutGroup = "union", repulsion = 0.9,
    nodeColor  = "feature", featVecCol = node_colors,
    nodeSize   = "eigenvector", nodeSizeSpread = 3,
    borderCol  = "gray40", labels = TRUE, nodeTransp = 50,
    negDiffCol = TRUE, posCol = "#2E8B57", negCol = "#FF6347",
    edgeWidth  = 2, edgeFilter = "threshold", edgeFilterPar = 0.10,
    edgeTranspLow = 80, edgeTranspHigh = 50,
    highlightHubs = FALSE, hubLabelThresh = 0.5, hubTransp = 30,
    cexHubLabels = 1.2
  )
  title(main = paste0("Cluster: ", sub("cluster__", "", cid)))
  print(myplot)
}

# ── Network-topology classifier ──────────────────────────────────────────────
classify_topology <- function(g) {
  m <- list()
  ## basic metrics
  m$n_nodes <- vcount(g);     m$n_edges <- ecount(g)
  m$density <- edge_density(g)
  
  ## centrality
  deg <- igraph::degree(g)
  m$avg_degree                 <- mean(deg)
  m$degree_centralization      <- centr_degree(g)$centralization
  m$betweenness_centralization <- centr_betw(g)$centralization
  m$closeness_centralization   <- centr_clo(g)$centralization
  
  ## connectivity
  m$avg_path_length <- average.path.length(g, directed = FALSE)
  m$diameter        <- diameter(g, directed = FALSE)
  
  ## clustering
  m$transitivity    <- transitivity(g, type = "global")
  m$avg_clustering  <- transitivity(g, type = "average")
  
  ## community structure
  if (vcount(g) > 2) {
    wtc                <- cluster_walktrap(g)
    m$modularity       <- modularity(wtc)
    m$n_communities    <- length(unique(membership(wtc)))
  } else {
    m$modularity <- NA; m$n_communities <- NA
  }
  
  ## hub dominance (top 10 %)
  if (vcount(g) >= 10) {
    top_hubs      <- names(sort(deg, TRUE))[1:ceiling(0.1 * vcount(g))]
    hub_strength  <- sum(strength(g, vids = top_hubs))
    m$hub_dominance <- hub_strength / sum(strength(g))
  } else m$hub_dominance <- NA
  
  ## qualitative class
  m$class <- "Other topology"
  if (!is.na(m$density)  && m$density  > 0.5)                 m$class <- "Dense core"
  if (!is.na(m$modularity) && m$modularity > 0.5 &&
      !is.na(m$n_communities) && m$n_communities >= 3)        m$class <- "Modular structure"
  if (!is.na(m$hub_dominance) && m$hub_dominance > 0.4)       m$class <- "Hub-dominated"
  if (!is.na(m$avg_path_length) && m$avg_path_length < 2)      m$class <- "Small-world like"
  if (!is.na(m$transitivity) && m$transitivity > 0.3)          m$class <- "Cliquish structure"
  
  m
}

# ── 2. CLI argument parsing ─────────────────────────────────────────────────
parser <- ArgumentParser()
message("Parsing arguments")

## Desired outputs
parser$add_argument("-mod_ana", "--module_analysis",
                    default = "Pooled",
                    help    = "Analyse modules? (Pooled/…)")

parser$add_argument("-mod_comp", "--module_comparisons",
                    default = "NO",
                    help    = "Compare modules? YES/NO")

parser$add_argument("-mod_clust", "--module_cluster",
                    default = "NO",
                    help    = "Define clusters? YES/NO")

## Matrix processing
parser$add_argument("-residues", "--input_residues", type = "character")
parser$add_argument("-covariates", "--list_of_covariates", type = "character")
parser$add_argument("-cov_types", "--list_of_covariate_types", type = "character")
parser$add_argument("-phe", "--pheno", default = "microbiome_taxonomic")
parser$add_argument("-met", "--method", default = "Shallow")

## Taxonomic
parser$add_argument("-tax", "--input_taxonomy", type = "character")
parser$add_argument("-tax_ranks", "--taxonomic_ranks",
                    default = "Phylum,Class,Order,Family,Genus,Species")

## Functional
parser$add_argument("-func_ranks", "--functional_ranks",
                    default = "EC1,EC2,EC3,EC4")
parser$add_argument("-samp_filt", "--tax_sample_filter", default = "NO")
parser$add_argument("-samps", "--tax_defined_samples_and_enterotypes",
                    type = "character")

## Enterotypes / clusters / network
parser$add_argument("-ent", "--enterotypess", type = "character")
parser$add_argument("-umap", "--umap",           type = "character")
parser$add_argument("-umap_res", "--umap_res",   type = "character")
parser$add_argument("-measure", "--network_measure", default = "bicor")
parser$add_argument("-outl", "--outliers", type = "double", default = 0.05)
parser$add_argument("-thresh", "--threshold", type = "double", default = 0.25)
parser$add_argument("-alpha", "--alpha", type = "double", default = 0.05)
parser$add_argument("-cluster", "--cluster_method",
                    default = "cluster_walktrap")
parser$add_argument("-hub", "--hub_definition",
                    default = "betweenness,eigenvector")

args <- parser$parse_args()

# ── 3. Pre-processing & basic checks ────────────────────────────────────────
module_analysis    <- args$module_analysis
module_comparisons <- args$module_comparisons
if (!module_comparisons %in% c("YES", "NO"))
  stop("module_comparisons must be 'YES' or 'NO'.")
module_cluster <- args$module_cluster

input_residues           <- args$input_residues
list_of_covariates       <- args$list_of_covariates
list_of_covariate_types  <- args$list_of_covariate_types
method                   <- args$method
methods                  <- factor(method)       # keep plural for expand.grid
pheno                    <- args$pheno

input_taxonomy  <- args$input_taxonomy
tax_ranks       <- args$taxonomic_ranks
func_ranks      <- args$functional_ranks
samp_filt       <- args$tax_sample_filter
samps           <- args$tax_defined_samples_and_enterotypes

umap_dir        <- args$umap
umap_res_dir    <- args$umap_res

network_measure <- args$network_measure
outliers        <- args$outliers
threshold       <- args$threshold
alpha           <- args$alpha
cluster_method  <- args$cluster_method
hub_def         <- args$hub_definition
hub_definition  <- na.omit(unlist(strsplit(hub_def, ",")))


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

if (module_analysis=="Pooled"){
  subset_ids="ALL"
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

# Read ranks
print("Reading ranks")
if (pheno=="microbiome_taxonomic"){
  ranks<-na.omit(c((sapply(strsplit(tax_ranks, ","),"[",1)),
                   (sapply(strsplit(tax_ranks, ","), "[",2)),
                   (sapply(strsplit(tax_ranks, ","), "[",3)),
                   (sapply(strsplit(tax_ranks, ","), "[",4)),
                   (sapply(strsplit(tax_ranks, ","), "[",5)),
                   (sapply(strsplit(tax_ranks, ","), "[",6)),
                   (sapply(strsplit(tax_ranks, ","), "[",7))))
} else {
  ranks<-na.omit(c((sapply(strsplit(func_ranks, ","),"[",1)),
                   (sapply(strsplit(func_ranks, ","), "[",2)),
                   (sapply(strsplit(func_ranks, ","), "[",3)),
                   (sapply(strsplit(func_ranks, ","), "[",4))))
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

# Enterotypes
load(args$enterotypes)

### End of Step 1 ###

### Step 2: Network estimation ###
print("Step 2: Network estimation")

if (pheno=="microbiome_taxonomic" && module_cluster=="YES"){
  # Read clusters information
  if (umap_dir != "" & umap_res_dir != "") {
    print("Reading clusters information")
    
    umap_files <- data.frame(file=list.files(umap_dir, pattern=".csv"))
    umap_files$subset=gsub(".csv","",umap_files$file)
    umap_files$subset= sapply(strsplit(umap_files$subset, "_"), "[",4)
    
    umap_res_files <- data.frame(file=list.files(umap_res_dir, pattern=".csv"))
    umap_res_files$subset=gsub(".csv","",umap_res_files$file)
  } else {
    stop("No clusters info provided")
  }
}


print("Starting to build the network")
if (("Species" %in% ranks) || ("EC4" %in% ranks)){
  for (subset_id in subset_ids){
    print(subset_id)
    ## Load the matrix
    if (pheno=="microbiome_taxonomic"){
      matrix=t(get(paste('residuals_qned_counts',"Species",subset_id,sep="_")))
    } else {
      matrix=t(get(paste('residuals_qned_counts',"EC4",subset_id,sep="_")))
    }
    
    if (pheno=="microbiome_taxonomic" && module_cluster=="YES"){
      clusters=read_csv_file(paste(umap_dir,umap_files$file[umap_files$subset==subset_id],sep="/"))
      samp_umap=read_csv_file(paste(umap_res_dir,umap_res_files$file[umap_res_files$subset==subset_id],sep="/"))
      if (pheno=="microbiome_taxonomic"){
        clusters$Traits=taxonomy$Species[match(samp_umap$Genome_id,taxonomy$Genome)]
      } else {
        clusters$Traits=samp_umap$Genome_id
      }
      clusters=clusters[clusters$Traits %in% colnames(matrix),]
      clusters$Cluster=gsub("-1","noise",clusters$Cluster)
      clusters$Cluster=paste("cluster",clusters$Cluster,sep="__")
      clusters_original=clusters
    }
    
    
    
    #################################################################################### Building full Network   

    ###############################################################################
    # 2.  Build the network with bicor
    ###############################################################################
    
    # Compute raw bicor matrix
    print("Compute raw bicor matrix")
    raw_bicor <- WGCNA::bicor(matrix, maxPOutliers = 0.2)
    # Check for negative values
    cat("Raw bicor matrix summary:\n")
    summary(c(raw_bicor))
    cat("Negative correlations:", sum(raw_bicor < 0, na.rm = TRUE), "\n")
    
    print("Build the network with bicor")
    ncores <- parallel::detectCores() - 1
    net_w <- netConstruct(
      data        = matrix,
      dataType    = "counts",
      filtTax          = "none",
      filtSamp         = "none",
      measure     = network_measure,
      measurePar  = list(maxPOutliers = outliers),
      sparsMethod = "threshold",
      thresh      = threshold,
      alpha       = alpha,
      adjust      = "adaptBH",
      dissFunc    = "signed",
      zeroMethod  = "none",
      normMethod  = "none",
      cores       = ncores,
      seed        = 123
    )
    
    # For single network
    check_network_res=check_network(net_w$assoMat)
    print(check_network_res)
    
    print("Imput zeroes to diagonals and define negative associations")
    for (nm in grep("^adjaMat", names(net_w), value = TRUE)) {
      diag(net_w[[nm]]) <- 0
    }
    cat("Association matrix summary:\n")
    summary(c(net_w$assoMat1))
    cat("Negative associations:", sum(net_w$assoMat1 < 0, na.rm = TRUE), "\n")
    
    print("Making adjustments for manual estimation of betweenness centrality")
    g <- graph_from_adjacency_matrix(net_w$adjaMat1, mode = "undirected",
                                     weighted = TRUE, diag = FALSE)
    E(g)$weight <- 1 / E(g)$weight
    #net_w$graph   <- g            # add it back for downstream functions
    
    print("Make the unweighted version")
    net_uw <- net_w
    net_uw$parameters$weighted <- FALSE
    
    # keep the thresholded matrix but turn it into 0/1
    bin_adj <- (net_w$adjaMat1 != 0) * 1     # TRUE→1, FALSE→0
    diag(bin_adj) <- 0                       # just in case
    
    g_uw <- graph_from_adjacency_matrix(
      bin_adj,
      mode     = "upper",           # or "plus"; see below
      weighted = NULL,
      diag     = FALSE)
    
    comp <- components(g_uw)                         # or clusters()
    lcc_id <- which.max(comp$csize)                  # largest component
    g_lcc  <- induced_subgraph(g_uw, which(comp$membership == lcc_id))
    cat("Edges seen by igraph:", ecount(g_lcc), "\n")        # should now be 207 584
    cat("Is graph connected? ", is.connected(g_lcc), "\n")   # optional
    metrics_lcc <- classify_topology(g_lcc)          # run on the LCC only
    print(metrics_lcc)
    
    
    print("Saving")
    save(net_w,
         net_uw,
         check_network_res,
         file ="networks.RData")
    
    #################################################################################### Extracting information from the full network  
    
    ###############################################################################
    # 3.  Analyse communities and hubs
    ###############################################################################
    print("Extract props: Weighted")
    ### –– 3. Analyse communities and hubs ––––––––––––––––––
    props_w <- netAnalyze(net_w, 
                          centrLCC = FALSE,
                          avDissIgnoreInf = TRUE,
                          sPathNorm = FALSE,
                          clustMethod = cluster_method,
                          hubPar = c(hub_definition),
                          hubQuant = 0.95,
                          lnormFit = TRUE,
                          normDeg = FALSE,
                          normBetw = FALSE,
                          normClose = FALSE,
                          normEigen = FALSE)
    
    # find every slot whose name starts with "assoMat"
    ## -------- 1.  locate every association matrix slot ------------------------
    asso_slots <- grep("^assoMat", names(props_w$input), value = TRUE)
    ## -------- 2.  drop self-loops (diag = 0)  ---------------------------------
    for (s in asso_slots) {
      mat <- props_w$input[[s]]
      if (is.null(mat)) {
        message(sprintf("slot %s is empty – skipped", s))
        next
      }
      if (!is.matrix(mat) && !inherits(mat, "Matrix"))
        stop(sprintf("slot %s is not a matrix-like object", s))
      mat <- as.matrix(mat)     # keeps dimnames
      diag(mat) <- 0            # remove self-edges
      props_w$input[[s]] <- mat   # write back
    }
    
    print("Extract props: Unweighted")
    ### –– 3. Analyse communities and hubs ––––––––––––––––––
    props_uw <- netAnalyze(net_uw, 
                           centrLCC = FALSE,
                           avDissIgnoreInf = TRUE,
                           sPathNorm = FALSE,
                           clustMethod = cluster_method,
                           hubPar = c(hub_definition),
                           hubQuant = 0.95,
                           lnormFit = TRUE,
                           normDeg = FALSE,
                           normBetw = FALSE,
                           normClose = FALSE,
                           normEigen = FALSE)
    
    # find every slot whose name starts with "assoMat"
    ## -------- 1.  locate every association matrix slot ------------------------
    asso_slots <- grep("^assoMat", names(props_uw$input), value = TRUE)
    ## -------- 2.  drop self-loops (diag = 0)  ---------------------------------
    for (s in asso_slots) {
      mat <- props_uw$input[[s]]
      if (is.null(mat)) {
        message(sprintf("slot %s is empty – skipped", s))
        next
      }
      if (!is.matrix(mat) && !inherits(mat, "Matrix"))
        stop(sprintf("slot %s is not a matrix-like object", s))
      mat <- as.matrix(mat)     # keeps dimnames
      diag(mat) <- 0            # remove self-edges
      props_uw$input[[s]] <- mat   # write back
    }
    
    print("Saving")
    save(net_w,
         net_uw,
         props_w,
         props_uw,
         check_network_res,
         file ="networks.RData")
    
    
    ###############################################################################
    # 4.  Tidy outputs: one fast join instead of many match() calls
    ###############################################################################
    print("Define centrality values")
    bw <- betweenness(g, weights = E(g)$weight)
    
    centralities <- props_w$centralities
    centr_vecs <- centralities[!sapply(centralities, is.null)]
    centralities   <- as.data.frame(do.call(cbind, centr_vecs),
                                    check.names = FALSE)
    colnames(centralities)=c("degree","betweenness","closeness","eigenvector")
    centralities$betweenness=bw
    n <- nrow(centralities)
    centralities$betweenness <- centralities$betweenness / ((n - 1) * (n - 2) / 2)
    #if (pheno=="microbiome_taxonomic"){
    centralities   <- tibble::rownames_to_column(centralities, "Traits")
    #centralities   <- dplyr::left_join(centralities, taxonomy, by = "Traits")
    #} else {
    #  centralities   <- tibble::rownames_to_column(centralities, "Traits")
    #}
    
    print("Saving")
    save(net_w,
         net_uw,
         props_w,
         props_uw,
         centralities,
         check_network_res,
         file ="networks.RData")
    
    print("Define modules within the network")
    modules <- props_w$clustering$clust1 %>%
      as.data.frame()
    colnames(modules)="Module"
    modules <- modules %>%
      tibble::rownames_to_column("Traits") %>%
      mutate(Module = paste0("module__", Module)) %>%
      left_join(centralities, by = "Traits")
    
    print("Saving")
    save(net_w,
         net_uw,
         props_w,
         props_uw,
         centralities,
         modules,
         check_network_res,
         file ="networks.RData")
    
    if (pheno=="microbiome_taxonomic" && module_cluster=="YES") {
      print("Map clusters to the network")
      clusters <- clusters_original %>%
        left_join(centralities, by = "Traits")
      
      print("Saving")
      save(net_w,
           net_uw,
           props_w,
           props_uw,
           centralities,
           modules,
           clusters,
           check_network_res,
           file ="networks.RData")
    }
    
    
    print("Estimate hubs for modules")
    # Merge centrality and module data
    combined_data <- modules
    
    # 1. Global Hub Classification (whole network) - CORRECTED
    global_hubs <- combined_data %>%
      mutate(
        # Calculate z-scores for each centrality measure
        across(c(hub_definition), 
               ~ as.numeric(scale(.)), 
               .names = "z_{.col}"),
        # Create composite centrality score (average of z-scores)
        global_composite = rowMeans(
          cbind(z_betweenness, z_eigenvector),  ########################### hardcoded !!!!!!!
          na.rm = TRUE
        ),
        # Classify as global hub if in top 10% of composite scores
        global_hub = global_composite >= quantile(global_composite, 0.90, na.rm = TRUE)
      )
    
    module_hubs <- global_hubs %>%
      group_by(Module) %>%
      mutate(module_size = n()) %>% 
      mutate(across(all_of(hub_definition),
                    ~ if (module_size[1] > 1) scale(.)[, 1] else 0,
                    .names = "modz_{.col}")) %>%
      rowwise() %>%
      mutate(module_composite = mean(c_across(starts_with("modz_")), na.rm = TRUE),
             module_hub       = if_else(module_size == 1,
                                        TRUE,
                                        module_composite >=
                                          quantile(module_composite, 0.80,
                                                   na.rm = TRUE))) %>%
      ungroup()
    
    # 3. Final Hub Classification
    hubs_modules <- module_hubs %>%
      dplyr::select(Traits, Module, global_hub, module_hub) %>%
      mutate(
        hub_type = case_when(
          global_hub & module_hub ~ "Global and Module Hub",
          global_hub ~ "Global Hub",
          module_hub ~ "Module Hub",
          TRUE ~ "Non-hub"
        )
      )
    
    if (pheno=="microbiome_taxonomic" && module_cluster=="YES"){
      print("Saving")
      save(net_w,
           net_uw,
           props_w,
           props_uw,
           centralities,
           modules,
           clusters,
           hubs_modules,
           check_network_res,
           file ="networks.RData")
    } else {
      save(net_w,
           net_uw,
           props_w,
           props_uw,
           centralities,
           modules,
           hubs_modules,
           check_network_res,
           file ="networks.RData")
    }
    
    if (pheno=="microbiome_taxonomic" && module_cluster=="YES"){
      print("Estimate hubs for clusters")
      # Merge centrality and clusters data
      global_hubs$Module=clusters$Cluster[match(global_hubs$Traits,clusters$Traits)]
      
      clusters_hubs <- global_hubs %>%
        group_by(Module) %>%
        mutate(clusters_size = n()) %>%
        mutate(across(all_of(hub_definition),
                      ~ if (clusters_size[1] > 1) scale(.)[, 1] else 0,
                      .names = "modz_{.col}")) %>%
        rowwise() %>%
        mutate(clusters_composite = mean(c_across(starts_with("modz_")), na.rm = TRUE),
               clusters_hub       = if_else(clusters_size == 1,
                                            TRUE,
                                            clusters_composite >=
                                              quantile(clusters_composite, 0.80,
                                                       na.rm = TRUE))) %>%
        ungroup()
      
      # 3. Final Cluster Classification
      hubs_clusters <- clusters_hubs %>%
        dplyr::select(Traits, Module, global_hub, clusters_hub) %>%
        mutate(
          hub_type = case_when(
            global_hub & clusters_hub ~ "Global and Cluster Hub",
            global_hub ~ "Global Hub",
            clusters_hub ~ "Cluster Hub",
            TRUE ~ "Non-hub"
          )
        )
      colnames(hubs_clusters)[2]="Cluster"
      
      print("Saving")
      save(net_w,
           net_uw,
           props_w,
           props_uw,
           centralities,
           modules,
           clusters,
           hubs_modules,
           hubs_clusters,
           check_network_res,
           file ="networks.RData")
    }
    
    
    #################################################################################### Print plots for full network
    
    print("Select colors for the modules")
    ext_ids <- sort(unique(modules$Module))
    pal_ext <- topo.colors(length(ext_ids))
    names(pal_ext) <- ext_ids
    node_cols <- setNames(pal_ext[modules$Module[match(colnames(matrix),modules$Traits)]],colnames(matrix))
    node_cols[is.na(node_cols)] <- "#DDDDDD"

    # Verify color vector matches node names
    all_nodes <- colnames(props_w$input$adjaMat1)
    cat("Color vector valid:", all(all_nodes %in% names(node_cols)), "\n")
    cat("Missing nodes:", setdiff(all_nodes, names(node_cols)), "\n")
    
    print("Printing the modules plots")
    topology_results=list()
    pdf("networks_per_module.pdf")
    for (mid in ext_ids) {
      print("Defining topologies per module")
      Traits_in_module <- modules$Traits[modules$Module == mid]
      if (length(Traits_in_module) < 2) {
        cat("Skipping module", mid, "- too few nodes\n")
        next
      }
      # Create subgraph
      subg <- induced_subgraph(g_uw, vids = Traits_in_module)
      comp <- components(subg)                         # or clusters()
      lcc_id <- which.max(comp$csize)                  # largest component
      subg_lcc  <- induced_subgraph(subg, which(comp$membership == lcc_id))
      # Classify topology
      metrics <- classify_topology(subg_lcc)
      # Store results
      topology_results[[mid]] <- data.frame(
        Module = mid,
        Nodes = metrics$n_nodes,
        Edges = metrics$n_edges,
        Density = round(metrics$density, 4),
        Avg_Degree = round(metrics$avg_degree, 2),
        Avg_Path_Length = round(metrics$avg_path_length, 2),
        Transitivity = round(metrics$transitivity, 3),
        Modularity = round(ifelse(is.na(metrics$modularity), 0, metrics$modularity), 3),
        Communities = ifelse(is.na(metrics$n_communities), 1, metrics$n_communities),
        Hub_Dominance = round(ifelse(is.na(metrics$hub_dominance), 0, metrics$hub_dominance), 3),
        Classification = metrics$class,
        stringsAsFactors = FALSE
      )
      print(paste0("Module topology: ",topology_results[[mid]]$Classification))
      cat("\n")
      print("Printing plots per module")
      plot_module(props_w, mid, modules, node_cols)
      print(mid)
    }
    dev.off()
    
    # Combine results into data frame
    print("Combine topology results")
    topology_modules <- do.call(rbind, topology_results)
    rownames(topology_modules) <- NULL

    #################################################################################### Print plots for full network (per cluster)
    
    if (pheno=="microbiome_taxonomic" && module_cluster=="YES"){
      print("Select colors for the clusters")
      ext_ids <- sort(unique(clusters$Cluster))
      pal_ext <- topo.colors(length(ext_ids))
      names(pal_ext) <- ext_ids
      node_cols <- setNames(pal_ext[clusters$Cluster[match(colnames(matrix),clusters$Traits)]],colnames(matrix))
      node_cols[is.na(node_cols)] <- "#DDDDDD"
      
      # Verify color vector matches node names
      all_nodes <- colnames(props_w$input$adjaMat1)
      cat("Color vector valid:", all(all_nodes %in% names(node_cols)), "\n")
      cat("Missing nodes:", setdiff(all_nodes, names(node_cols)), "\n")
      
      print("Defining topologies and printing plots per cluster")
      topology_results <- list()
      pdf("networks_per_cluster.pdf")
      for (cid in ext_ids){
        
        print("Defining topologies per cluster")
        Traits_in_cluster <- clusters$Traits[clusters$Cluster == cid]
        if (length(Traits_in_cluster) < 2) {
          cat("Skipping cluster", cid, "- too few nodes\n")
          next
        }
        # Create subgraph
        subg <- induced_subgraph(g_uw, vids = Traits_in_cluster)
        comp <- components(subg)                         # or clusters()
        lcc_id <- which.max(comp$csize)                  # largest component
        subg_lcc  <- induced_subgraph(subg, which(comp$membership == lcc_id))
        # Classify topology
        metrics <- classify_topology(subg_lcc)
        # Store results
        topology_results[[cid]] <- data.frame(
          Cluster = cid,
          Nodes = metrics$n_nodes,
          Edges = metrics$n_edges,
          Density = round(metrics$density, 4),
          Avg_Degree = round(metrics$avg_degree, 2),
          Avg_Path_Length = round(metrics$avg_path_length, 2),
          Transitivity = round(metrics$transitivity, 3),
          Modularity = round(ifelse(is.na(metrics$modularity), 0, metrics$modularity), 3),
          Communities = ifelse(is.na(metrics$n_communities), 1, metrics$n_communities),
          Hub_Dominance = round(ifelse(is.na(metrics$hub_dominance), 0, metrics$hub_dominance), 3),
          Classification = metrics$class,
          stringsAsFactors = FALSE
        )
        print(paste0("Cluster topology: ",topology_results[[cid]]$Classification))
        cat("\n")
        print("Printing plots per cluster")
        plot_cluster(props_w, cid, clusters, node_cols)
        
        print(cid)
      }
      dev.off()
      
      # Combine results into data frame
      print("Combine topology results")
      topology_clusters <- do.call(rbind, topology_results)
      rownames(topology_clusters) <- NULL
      
      save(net_w,
           net_uw,
           props_w,
           props_uw,
           centralities,
           modules,
           clusters,
           hubs_modules,
           hubs_clusters,
           topology_modules,
           topology_clusters,
           check_network_res,
           file ="networks.RData")
      
      #################################################################################### Print plots for full network (per cluster pair)
      
      #if (pheno=="microbiome_taxonomic" && module_cluster=="YES"){
      #  if (subset_id=="ALL"){
      #    print("Printing pairs of clusters")
      #    clusts=list(
      #      p1=c("cluster__3","cluster__0"),
      #      p2=c("cluster__3","cluster__8"),
      #      p3=c("cluster__3","cluster__4"),
      #      p4=c("cluster__3","cluster__15"),
      #      p5=c("cluster__3","cluster__10"),
      #      p6=c("cluster__4","cluster__15"),
      #      p7=c("cluster__4","cluster__8"),
      #      p8=c("cluster__8","cluster__15"),
      #      p9=c("cluster__8","cluster__10"),
      #      p10=c("cluster__1","cluster__2"),
      #      p11=c("cluster__1","cluster__3"),
      #      p12=c("cluster__2","cluster__3")
      #    )
      #    
      #    pdf("networks_per_cluster_pair.pdf")
      #    for(i in 1:length(clusts)){
      #      
      #      print("Defining topologies per cluster pair")
      #      Traits_in_cluster <- clusters$Traits[clusters$Cluster %in% clusts[[i]]]
      #      if (length(Traits_in_cluster) < 2) {
      #        cat("Skipping cluster", cid, "- too few nodes\n")
      #        next
      #      }
      #      # Create subgraph
      #      subg <- induced_subgraph(g_uw, vids = Traits_in_cluster)
      #      comp <- components(subg)                         # or clusters()
      #      lcc_id <- which.max(comp$csize)                  # largest component
      #      subg_lcc  <- induced_subgraph(subg, which(comp$membership == lcc_id))
      #      # Classify topology
      #      metrics <- classify_topology(subg_lcc)
      #      # Store results
      #      topology_results[[cid]] <- data.frame(
      #        Cluster = cid,
      #        Nodes = metrics$n_nodes,
      #        Edges = metrics$n_edges,
      #        Density = round(metrics$density, 4),
      #        Avg_Degree = round(metrics$avg_degree, 2),
      #        Avg_Path_Length = round(metrics$avg_path_length, 2),
      #        Transitivity = round(metrics$transitivity, 3),
      #        Modularity = round(ifelse(is.na(metrics$modularity), 0, metrics$modularity), 3),
      #        Communities = ifelse(is.na(metrics$n_communities), 1, metrics$n_communities),
      #        Hub_Dominance = round(ifelse(is.na(metrics$hub_dominance), 0, metrics$hub_dominance), 3),
      #        Classification = metrics$class,
      #        stringsAsFactors = FALSE
      #      )
      #      
      #      print("Printing plots per cluster")
      #      keep <- Traits_in_cluster
      #      myplot=plot(  props_w,
      #                    rmLoops = TRUE,
      #                    selNet = 1,
      #                    ## layout & general look
      #                    sameLayout      = TRUE,
      #                    layoutGroup     = "union",
      #                    repulsion       = 0.9,
      #                    mar             = c(1, 1, 3, 1),
      #                    showTitle       = TRUE,
      #                    cexTitle        = 2.8,
      #                    
      #                    ## node aesthetics
      #                    nodeColor       = "feature",      # << tell it to use our own vector
      #                    featVecCol        = node_cols,  # << exactly one entry per node
      #                    nodeSize        = "eigenvector",
      #                    nodeSizeSpread  = 3,
      #                    nodeTransp = 50, 
      #                    hubTransp = 30,
      #                    borderCol       = "gray40",
      #                    labels          = FALSE,
      #                    rmSingles       = "inboth",
      #                    nodeFilter     = "names",          # keep only these node names …
      #                    nodeFilterPar  = keep,       # … that belong to module_id
      #                    
      #                    ## *** edge aesthetics ****************************************************
      #                    ##      – different colours for positive / negative
      #                    negDiffCol = TRUE,
      #                    posCol = "#2E8B57",   # or any single colour
      #                    negCol = "#FF6347",
      #                    edgeFilter     = "threshold",
      #                    edgeFilterPar  = 0.10,
      #                    edgeWidth       = 2,                          # thicker lines
      #                    edgeTranspLow = 80, 
      #                    edgeTranspHigh = 50,
      #      )
      #      print(myplot)
      #      title(main = paste0("Clusters: ",paste(gsub("cluster__","",clusts[[i]]),collapse=", ")), line = 0.5, cex.main = 2.8)
      #      print(paste0("Clusters: ",paste(gsub("cluster__","",clusts[[i]]),collapse=", ")))
      #    }
      #    dev.off()
      #    # Combine results into data frame
      #    print("Combine topology results")
      #    topology_clusters_pairs <- do.call(rbind, topology_results)
      #    rownames(topology_clusters_pairs) <- NULL
      #  }
      #} else {topology_clusters_pairs=c() }
    }
    
    print("Defining topologies of the full network")
    metrics <- classify_topology(g_lcc)
    # Store results
    topology_full_network <- data.frame(
      Nodes = metrics$n_nodes,
      Edges = metrics$n_edges,
      Density = round(metrics$density, 4),
      Avg_Degree = round(metrics$avg_degree, 2),
      Avg_Path_Length = round(metrics$avg_path_length, 2),
      Transitivity = round(metrics$transitivity, 3),
      Modularity = round(ifelse(is.na(metrics$modularity), 0, metrics$modularity), 3),
      Communities = ifelse(is.na(metrics$n_communities), 1, metrics$n_communities),
      Hub_Dominance = round(ifelse(is.na(metrics$hub_dominance), 0, metrics$hub_dominance), 3),
      Classification = metrics$class,
      stringsAsFactors = FALSE
    )
    print(paste0("Network topology: ",topology_full_network$Classification))
    
    print("Identify LCCs")
    lcc_nodes  <- props_w$lccNames1
    
    if (pheno=="microbiome_taxonomic" && module_cluster=="YES"){
      print("Saving")
      save(net_w,
           net_uw,
           props_w,
           props_uw,
           centralities,
           modules,
           clusters,
           hubs_modules,
           hubs_clusters,
           topology_modules,
           topology_clusters,
           #topology_clusters_pairs,
           topology_full_network,
           lcc_nodes,
           check_network_res,
           file ="networks.RData")
    } else {
      print("Saving")
      save(net_w,
           net_uw,
           props_w,
           props_uw,
           centralities,
           modules,
           hubs_modules,
           topology_modules,
           topology_full_network,
           lcc_nodes,
           check_network_res,
           file ="networks.RData")
    }
    
    
    #################################################################################### Print plots for full network (entire community)
    
    # Plot the full network
    print("Plot the network")
    pdf('networks.pdf')
    myplot=plot(  props_w,
                  rmLoops = TRUE,
                  selNet = 1,
                  ## layout & general look
                  sameLayout      = TRUE,
                  layoutGroup     = "union",
                  repulsion       = 0.9,
                  mar             = c(1, 1, 3, 1),
                  showTitle       = TRUE,
                  cexTitle        = 2.8,
                  
                  ## node aesthetics
                  nodeColor       = "feature",      # << tell it to use our own vector
                  featVecCol        = node_cols,  # << exactly one entry per node
                  nodeSize        = "eigenvector",
                  nodeSizeSpread  = 3,
                  borderCol       = "gray40",
                  labels          = FALSE,
                  rmSingles       = "inboth",
                  nodeFilter      = "clustMin",
                  nodeFilterPar   = 10,
                  nodeTransp = 50, 
                  hubTransp = 30,
                  
                  ## *** edge aesthetics ****************************************************
                  ##      – different colours for positive / negative
                  negDiffCol = TRUE,
                  posCol = "#2E8B57",   # or any single colour
                  negCol = "#FF6347",
                  edgeFilter     = "threshold",
                  edgeFilterPar  = 0.10,
                  edgeWidth       = 2,                          # thicker lines
                  edgeTranspLow = 80, 
                  edgeTranspHigh = 50,
    )
    print(myplot)
    title(main = paste0("Network: ",subset_id), line = 0.5, cex.main = 2.8)
    dev.off()
    
  }
  
  #################################################################################### Define the networks per enterotype
  
  print("Load enterotypes")
  for (i in 1:length(subset_ids)){
    if (samp_filt=="YES" && samps !="" && pheno!="microbiome_taxonomic"){
      ent_samples=read.csv(samps,sep=",")
    } else if (samp_filt=="YES" && samps =="" && pheno!="microbiome_taxonomic"){
      stop("If samples will be filtered based on previous taxonomic profiling, please inform the path to the file with such information")
    } else {
      ent_samples=enterotypes_samples_objs[[i]]
    }
    assign(paste('ent_samples',subset_ids[i],sep='_'), ent_samples)
  }
  ent_names=sort(as.character(unique(ent_samples$Enterotype)))
  
  
  print("Compare the networks of the two enterotypes")
  print("Make the list")
  mat_list=list()
  for (ent in ent_names){
    idx          <- match(ent_samples$sample[ent_samples$Enterotype == ent], rownames(matrix))
    ent_matrix   <- matrix[idx, ]
    mat_list[[ent]] <- ent_matrix
    # Compute raw bicor matrix
    print("Compute raw bicor matrix")
    raw_bicor <- WGCNA::bicor(ent_matrix, maxPOutliers = 0.2)
    # Check for negative values
    cat("Raw bicor matrix summary:\n")
    summary(c(raw_bicor))
    cat("Negative correlations:", sum(raw_bicor < 0, na.rm = TRUE), "\n")
  }

  ###############################################################################
  # 2.  Build the network with bicor
  ###############################################################################
  
  print("Build the network with bicor")
  ncores <- parallel::detectCores() - 1
  net_w <- netConstruct(
    data        = mat_list[[1]],
    data2        = mat_list[[2]],
    dataType    = "counts",
    filtTax          = "none",
    filtSamp         = "none",
    measure     = network_measure,
    measurePar  = list(maxPOutliers = outliers),
    sparsMethod = "threshold",
    thresh      = threshold,
    alpha       = alpha,
    adjust      = "adaptBH",
    dissFunc    = "signed",
    zeroMethod  = "none",
    normMethod  = "none",
    cores       = ncores,
    seed        = 123
  )
  
  print("Imput zeroes to diagonals and define negative associations")
  for (nm in grep("^adjaMat", names(net_w), value = TRUE)) {
    diag(net_w[[nm]]) <- 0
  }
  cat("Association matrix summary:\n")
  summary(c(net_w$assoMat1))
  cat("Negative associations:", sum(net_w$assoMat1 < 0, na.rm = TRUE), "\n")
  
  print("Making adjustments for manual estimation of betweenness centrality")
  g <- graph_from_adjacency_matrix(net_w$adjaMat1, mode = "undirected",
                                   weighted = TRUE, diag = FALSE)
  E(g)$weight <- 1 / E(g)$weight
  #net_w$graph   <- g            # add it back for downstream functions
  
  print("Make the unweighted version")
  net_uw <- net_w
  net_uw$parameters$weighted <- FALSE
  
  # keep the thresholded matrix but turn it into 0/1
  bin_adj <- (net_w$adjaMat1 != 0) * 1     # TRUE→1, FALSE→0
  diag(bin_adj) <- 0                       # just in case
  
  g_uw <- graph_from_adjacency_matrix(
    bin_adj,
    mode     = "upper",           # or "plus"; see below
    weighted = NULL,
    diag     = FALSE)
  
  comp <- components(g_uw)                         # or clusters()
  lcc_id <- which.max(comp$csize)                  # largest component
  g_lcc  <- induced_subgraph(g_uw, which(comp$membership == lcc_id))
  cat("Edges seen by igraph:", ecount(g_lcc), "\n")        # should now be 207 584
  cat("Is graph connected? ", is.connected(g_lcc), "\n")   # optional
  metrics_lcc <- classify_topology(g_lcc)          # run on the LCC only
  print(metrics_lcc)
  
  print("Saving")
  save(net_w,
       net_uw,
       file ="networks_comp_ent.RData")
  
  #################################################################################### Extracting information from the full network  
  
  ###############################################################################
  # 3.  Analyse communities and hubs
  ###############################################################################
  print("Extract props: Weighted")
  ### –– 3. Analyse communities and hubs ––––––––––––––––––
  props_w <- netAnalyze(net_w, 
                        centrLCC = FALSE,
                        avDissIgnoreInf = TRUE,
                        sPathNorm = FALSE,
                        clustMethod = cluster_method,
                        hubPar = c(hub_definition),
                        hubQuant = 0.95,
                        lnormFit = TRUE,
                        normDeg = FALSE,
                        normBetw = FALSE,
                        normClose = FALSE,
                        normEigen = FALSE)
  
  # find every slot whose name starts with "assoMat"
  ## -------- 1.  locate every association matrix slot ------------------------
  asso_slots <- grep("^assoMat", names(props_w$input), value = TRUE)
  ## -------- 2.  drop self-loops (diag = 0)  ---------------------------------
  for (s in asso_slots) {
    mat <- props_w$input[[s]]
    if (is.null(mat)) {
      message(sprintf("slot %s is empty – skipped", s))
      next
    }
    if (!is.matrix(mat) && !inherits(mat, "Matrix"))
      stop(sprintf("slot %s is not a matrix-like object", s))
    mat <- as.matrix(mat)     # keeps dimnames
    diag(mat) <- 0            # remove self-edges
    props_w$input[[s]] <- mat   # write back
  }
  
  print("Extract props: Unweighted")
  ### –– 3. Analyse communities and hubs ––––––––––––––––––
  props_uw <- netAnalyze(net_uw, 
                         centrLCC = FALSE,
                         avDissIgnoreInf = TRUE,
                         sPathNorm = FALSE,
                         clustMethod = cluster_method,
                         hubPar = c(hub_definition),
                         hubQuant = 0.95,
                         lnormFit = TRUE,
                         normDeg = FALSE,
                         normBetw = FALSE,
                         normClose = FALSE,
                         normEigen = FALSE)
  
  # find every slot whose name starts with "assoMat"
  ## -------- 1.  locate every association matrix slot ------------------------
  asso_slots <- grep("^assoMat", names(props_uw$input), value = TRUE)
  ## -------- 2.  drop self-loops (diag = 0)  ---------------------------------
  for (s in asso_slots) {
    mat <- props_uw$input[[s]]
    if (is.null(mat)) {
      message(sprintf("slot %s is empty – skipped", s))
      next
    }
    if (!is.matrix(mat) && !inherits(mat, "Matrix"))
      stop(sprintf("slot %s is not a matrix-like object", s))
    mat <- as.matrix(mat)     # keeps dimnames
    diag(mat) <- 0            # remove self-edges
    props_uw$input[[s]] <- mat   # write back
  }
  
  print("Saving")
  save(net_w,
       net_uw,
       props_w,
       props_uw,
       file ="networks_comp_ent.RData")
  
  print("Compare networks: Weighted")
  cmp_w <- netCompare(
    x = props_w,
    permTest = F,
    adjust = "adaptBH",
    cores = ncores,
    seed = 123
  )
  
  print("Compare networks: Unweighted")
  cmp_uw <- netCompare(
    x = props_uw,
    permTest = F,
    adjust = "adaptBH",
    cores = ncores,
    seed = 123
  )
  
  print("Saving")
  save(net_w,
       net_uw,
       props_w,
       props_uw,
       cmp_w,
       cmp_uw,
       file ="networks_comp_ent.RData")
  
  # Global properties
  print("Extract global properties")
  glob      <- stack(cmp_uw$diffGlobal)           # whole graph
  glob_lcc  <- stack(cmp_uw$diffGlobalLCC)        # just the largest connected component
  print("Summary of all global differences")
  print(glob)                           # easy data-frame view
  print("Summary of largest connected component (LCC)")
  print(glob_lcc)  
  
  # Differential centrality
  print("Extract differential centrality")
  diff_centrality <- cmp_w$diffCentr
  #print(diff_centrality)
  
  # LCCs
  print("Identify LCCs")
  g1 <- graph_from_adjacency_matrix(cmp_w$adjaMatrices$adja1,
                                    mode = "undirected", weighted = TRUE)
  which_comp <- components(g1)$membership
  lcc_nodes  <- names(which_comp[which_comp == which.max(table(which_comp))])
  
  print("Extract Jaccard results")
  jaccard_results <- data.frame(
    Measure = c("degree", "betweenness", "closeness", "eigenvector", "hub"),
    Jaccard = c(
      cmp_w$jaccDeg["jacc"],
      cmp_w$jaccBetw["jacc"],
      cmp_w$jaccClose["jacc"],
      cmp_w$jaccEigen["jacc"],
      cmp_w$jaccHub["jacc"]
    ),
    p_greater = c(
      cmp_w$jaccDeg["p.greater"],
      cmp_w$jaccBetw["p.greater"],
      cmp_w$jaccClose["p.greater"],
      cmp_w$jaccEigen["p.greater"],
      cmp_w$jaccHub["p.greater"]
    ),
    p_less = c(
      cmp_w$jaccDeg["p.less"],
      cmp_w$jaccBetw["p.less"],
      cmp_w$jaccClose["p.less"],
      cmp_w$jaccEigen["p.less"],
      cmp_w$jaccHub["p.less"]
    ))
  print(jaccard_results)
  
  print("Rand Index: clustering similarity")
  rand_results <- data.frame(
    Scope = c("Whole Network", "LCC"),
    Rand_Index = c(cmp_w$randInd["value"], cmp_w$randIndLCC["value"]),
    p_value = c(cmp_w$randInd["pval"], cmp_w$randIndLCC["pval"])
  )
  print(rand_results)
  
  print("Graphlet correlation distance")
  gcd_values <- data.frame(
    Scope = c("Whole Network", "LCC"),
    GCD = c(cmp_uw$gcd$gcd, cmp_uw$gcdLCC$gcd)
  )
  print(gcd_values)
  
  # Compute differential network
  print("Compute differential network")
  diff_net <- diffnet(
    net_w,                      # <-- one argument, a microNet object
    diffMethod = "fisher",    # Fisher’s z-test
    adjust     = "BH",        # Benjamini–Hochberg
    cores      = ncores,
    seed       = 123
  )
  
  print("Saving")
  save(net_w,
       net_uw,
       props_w,
       props_uw,
       cmp_w,
       cmp_uw,
       glob,
       glob_lcc,
       diff_centrality,
       lcc_nodes,
       jaccard_results,
       rand_results,
       gcd_values,
       diff_net,
       file="networks_comp_ent.RData"
  )
  
  # Plot networks
  print("Plot networks")
  pdf("networks_comp_ent.pdf")
  myplot=plot(  props_w,
                rmLoops = TRUE,
                selNet = 1,
                ## layout & general look
                sameLayout      = TRUE,
                layoutGroup     = "union",
                repulsion       = 0.9,
                mar             = c(1, 1, 3, 1),
                showTitle       = TRUE,
                cexTitle        = 2.8,
                
                ## node aesthetics
                nodeColor       = "feature",      # << tell it to use our own vector
                featVecCol        = node_cols,  # << exactly one entry per node
                nodeSize        = "eigenvector",
                nodeSizeSpread  = 3,
                borderCol       = "gray40",
                labels          = FALSE,
                rmSingles       = "inboth",
                nodeFilter      = "clustMin",
                nodeFilterPar   = 10,
                nodeTransp = 50, 
                highlightHubs = FALSE,
                hubTransp = 30,
                groupNames = ent_names,
                
                ## *** edge aesthetics ****************************************************
                ##      – different colours for positive / negative
                negDiffCol = TRUE,
                posCol = "#2E8B57",   # or any single colour
                negCol = "#FF6347",
                edgeFilter     = "threshold",
                edgeFilterPar  = 0.10,
                edgeWidth       = 2,                          # thicker lines
                edgeTranspLow = 80, 
                edgeTranspHigh = 50,
  )
  print(myplot)
  dev.off()
  
} else {
  stop("To define networks, please define it on the lowest level (Species or EC4).", "\n")
}

# ── Wrap-up ──────────────────────────────────────────────────────────────────
disableWGCNAThreads()