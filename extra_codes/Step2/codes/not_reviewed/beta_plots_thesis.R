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

# 3. Clean ids
clean_ids<-function(oldids,metadata) {
  raw_ids<-str_split(oldids, pattern="_", simplify = TRUE)
  for (i in 1:ncol(raw_ids)){
    w=which(metadata[[sample_identifier]] == raw_ids[1,i])
    if (length(w) != 0){
      ids<-sapply(strsplit(oldids, "_"), "[",i)
    }
  }
  return(ids)
}

input_meta="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/metadata_Shallow.txt"
list_of_covariates="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/covariates_Shallow.txt"
list_of_covariate_types="/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/input/covariates_classification.txt"
sample_identifier= "host_subject_id"
module_analysis="Pooled"

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

print("Load enterotypes")
load("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Tax_Mat_Net/Enterotypes/enterotypes.RData")
number=1
ent_samples=enterotypes_samples_objs[[number]]
ent_samples$Enterotype=gsub("enterotype__","Enterotype ",ent_samples$Enterotype)
metadata$Enterotype=ent_samples$Enterotype[match(metadata$host_subject_id,ent_samples$sample)]
covariates=c(covariates,"Enterotype")

cov_colors=list()
for (i in 1:length(covariates)){
  metadata[[covariates[i]]]=factor(metadata[[covariates[i]]])
  total_covs = sort(unique(metadata[[covariates[i]]]))
  if (covariates[i]=="Batch"){
    colours = metafolio::gg_color_hue(n =  length(total_covs), c = 100, l = 65)
    shuffled_indices <- sample(length(colours))
    colours <- colours[shuffled_indices]
  } else if (covariates[i]=="Sex"){
    colours = c("darkgreen","lightgreen")
  } else if (covariates[i]=="Study"){
    colours = c("darkblue","lightblue")
  } else {colours = c("darkred","orange")}
  names(colours)=total_covs
  cov_colors[[covariates[i]]] <- colours
}

methods=c("Tax","Func")
steps=c("With fixed effects","Without fixed effects")


pdf("/users/abaud/fmorillo/paper_figures/beta_plots_thesis.pdf")
permanova.total=data.frame()
for (method in methods){
  if (method=="Func"){
    text="EC4"
  } else {text="Species"}
  for (step in steps){
    if (step==steps[1]){
      load(paste0("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_",method,"_Mat_Net/Diversity/beta_diversity.RData"))
      gp.dist=gp.dist_objs[[4]]
    } else {
      load(paste0("/nfs/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_",method,"_Mat_Net/Diversity/residuals_qned_counts_beta.RData"))
      gp.dist=gp.dist_objs[[2]]
    }
    samples= attr(gp.dist, "Labels")
    this_metadata=metadata[metadata[[sample_identifier]] %in% samples,]
    rownames(this_metadata)=this_metadata[[sample_identifier]]
    this_metadata=this_metadata[samples,]
    gp.pcoa <- ape::pcoa(gp.dist)
    gp.pc <- data.frame(gp.pcoa$vectors)
    gp.pc <- gp.pc[rownames(this_metadata),]
    for (covariate in covariates){
      print(paste(step,text,covariate,sep="/"))
      formula <- as.formula(paste("gp.dist ~", covariate))
      permanova <- adonis2(formula, data = this_metadata, permutations = 999)
      permanova$step=step
      permanova$method=method
      permanova$covariate=covariate
      permanova.total=rbind(permanova.total,permanova)
      R2 <- round(permanova[1,3],digits=2)
      p_value <- round(permanova[1,5],digits=3)
      myplot<-ggplot(gp.pc, aes(x=Axis.1, y=Axis.2, colour = as.factor(this_metadata[[covariate]]))) +
        geom_point() + 
        ggtitle(paste0("Principal Coordinates by ",covariate,": ",step,"/",text)) +
        scale_colour_manual(values = cov_colors[[covariate]]) +
        xlab("PCo1") +
        ylab("PCo2") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, size = 14),
              axis.line = element_line(),
              plot.margin = margin(5.5,5.5, 5.5, 6.5, "pt"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              legend.text=element_text(size=12),
              legend.title = element_blank(),
              legend.position="bottom") +
        annotate("label", 
                 x = max(gp.pc$Axis.1), 
                 y = max(gp.pc$Axis.2), 
                 label = paste0("p = ",p_value),
                 hjust = 1, vjust = 1,
                 size = 5,
                 label.size = 0.5,
                 label.r = unit(0.2, "lines"),
                 label.padding = unit(0.25, "lines"),
                 fill = "white",
                 color = "black")
      print(myplot)
    }
  }
}
dev.off()

save(permanova.total,file="~/paper_figuers/permanova_beta_all.RData")
