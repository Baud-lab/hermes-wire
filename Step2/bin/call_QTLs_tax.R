#!/usr/bin/env Rscript

# Packages

suppressMessages(library(rhdf5))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(argparse))
suppressMessages(library(ggplot2))
suppressMessages(library(xfun))
suppressMessages(library(qqman))

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

# 2. Assign covariates to subsets
assignCovariate <- function(x) {
  subset_id <- x[1]
  covnames <- str_split(x[2], pattern=",", simplify = TRUE) 
  assign(paste('covariates',subset_id,sep='_'), covnames, envir = .GlobalEnv)
}

# -----------------------------------------------------------------------------


#Create parser object
parser <- ArgumentParser()
print("Parsing arguments")

#Define desired outputs:
#GLOBAL FEATURES:
parser$add_argument("-odir", "--out_dir", type="character", help="Output directory")
parser$add_argument("-model", "--model_type", type="character", help="Name of the model",default="uni")
parser$add_argument("-pheno", "--phenotype_version", type="character", help="Phenotype version",default="bacterial_taxa")
parser$add_argument("-kinloco", "--kinship_version_loco", type="character", help="Name of the kinship model for GRM LOCO",default="pruned_dosages")
parser$add_argument("-eff", "--effect", type="character", help="Effect to be checked")
parser$add_argument("-sig_l", "--significance_logp", type="double", help="Significance level",default=5.8)
parser$add_argument("-min_l", "--min_logp", type="double", help="Minimum threshold for -logP values",default=4)
parser$add_argument("-covariates", "--list_of_covariates", type="character", help="File with subsets and covariates info")
parser$add_argument("-ranks", "--taxonomic_ranks", type="character", help="Taxonomic ranks to be analysed",default="Phylum,Class,Order,Family,Genus,Species")
parser$add_argument("-res", "--residuals", type="character", help="R object with all matrices containing residuals")


#Get command line options, if help option encountered - print help and exit:
args <- parser$parse_args()

## Arguments
OUT_DIR=args$out_dir
MODEL=args$model_type
EFF=args$effect
GRM_V=args$kinship_version_loco
PHENO_V=args$phenotype_version 
suggestive_l=args$significance_logp
min_l=args$min_logp
ranks<-args$taxonomic_ranks
allowed <- c("Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV")
string_vector <- strsplit(ranks, ",")[[1]]
if (all(string_vector %in% allowed)) {
  print("All ranks are valid")
} else {
  invalid <- string_vector[!string_vector %in% allowed]
  stop("Invalid taxonomic rank(s) found: ", paste(invalid, collapse = ", "),"\n","Ranks should be: Phylum, Class, Order, Family, Genus, Species or ASV" ,"\n")
}
residuals=args$residuals
list_of_covariates<-args$list_of_covariates

if (EFF == 'DGE,cageEffect') effect = 'DGE_cageEffect' else if (EFF == 'cageEffect') effect = 'cageEffect'

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
    print("All elements in the subsets are present in the subset 'ALL'.\n")
  } else {
    print("Warning: Covariates ", paste(missing_covs, collapse = ", "), " are missing in the subset 'ALL' but used on other subsets.", "\n")
  }
} else {
  stop("No covariates file provided")
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


# Read residuals
if (residuals ==""){
  stop("Please, provide a RData file with the residuals if you are not using module 'Matrix Processing'.")
} else {
  print("Reading residuals")
  load(residuals)
  permutations <- expand.grid(subset_ids,ranks)
  concatenated <- paste(permutations$Var2, permutations$Var1, sep = "_")
  for (i in 1:length(concatenated)){
    residuals_qned_counts=residuals_qned_counts_objs[[i]]
    assign(paste('residuals_qned_counts',concatenated[i],sep='_'), residuals_qned_counts)
  }
  #sig_table=data.frame(subset = subset_ids, study_wide_significance = 0)
  for (subset_id in subset_ids){
    all_residuals_counts=data.frame()
    for (rank in ranks){
      residuals_qned_counts=get(paste("residuals_qned_counts",rank,subset_id,sep="_"))
      all_residuals_counts=rbind(all_residuals_counts,residuals_qned_counts)
    }
    assign(paste('all_residuals_counts',subset_id,sep='_'), all_residuals_counts)
  }
  for (subset_id in subset_ids){
    all_residuals_counts=get(paste('all_residuals_counts',subset_id,sep='_'))
    pca=prcomp(all_residuals_counts,center=TRUE,scale=TRUE)$sdev
    total_var = sum(pca^2)
    for (k in 1:length(pca)){
      if (sum((pca[1:k])^2) >= 0.95*total_var) {
        sig_l=-log10(10^(-suggestive_l)/k)
        #sig_table$study-wide_significance[sig_table$subset == subset_id] <- sig_l
        assign(paste('sig_l',subset_id,sep='_'), sig_l)
        cat(subset_id,": K is: ",k," / Significance level is: ",sig_l,"\n")
        break
      }
    }
  }
}

### End of Step 1 ###

### Step 2 ####
print("Step 2: Build QTLs table")

## Load model "Full"
pvalues_dir = paste0(OUT_DIR,"/",MODEL,'variate/',PHENO_V,"/pvalues_LOCO/",GRM_V,"_",effect,"/")
files = list.files(pvalues_dir)
files = files[grep('h5',files)]
length(files)

alpha = 10^(-min_l)

my_f = function(file) {
  pheno_name = sub('.h5','',file,fixed=T)
  print(pheno_name)
  
  my_h5 = try(h5read(paste(pvalues_dir,'/',pheno_name,'.h5',sep=''),'/'))
  if (inherits(my_h5,'try-error')) return('Error')
  
  my_pvalues = c()
  for (chr in 1:20) {
    #if (!paste('pvalues_chr',chr,sep='') %in% names(my_h5)) return('Error')
    my_pvalues = c(my_pvalues, my_h5[[paste('pvalues_chr',chr,sep='')]])
  }
  my_h5$pvalues = my_pvalues
  my_h5 = data.frame(my_h5[['chr']],my_h5[['pos']],my_h5[['pvalues']])
  colnames(my_h5) = c('chr','pos','pvalues')
  my_h5$id = paste(my_h5[,'chr'],my_h5[,'pos'],sep='_')
  
  # handle missing (-999) P values
  missing_chrs = unique(my_h5[my_h5$pvalues == (-999),'chr'])
  if (all(my_h5$chr %in% missing_chrs)) {
    print(paste('No non-missing P values for file', file))
    return(NULL)
  } else {
    if (length(missing_chrs)>0){
      print(paste("P values for chromosomes(s)", missing_chrs, "are missing"))
    }
  }
  my_h5 = my_h5[!my_h5$chr %in% missing_chrs,]
  
  my_h5$logP=-log10(my_h5[,'pvalues'])
  my_h5=my_h5[order(my_h5[,'pvalues']),]
  
  
  if (my_h5[1,'pvalues']>alpha) return(NULL)
  
  peaks_marker =c()
  peaks_pos=c()
  peaks_chr=c()
  peaks_pvalue=c()
  peaks_logP=c()
  peaks_ci_starts=c()
  peaks_ci_stops=c()	
  over = F
  while (over==F) {
    qtl_chr=my_h5[1,'chr']
    qtl_pos=my_h5[1,'pos']
    qtl_marker=my_h5[1,'id']
    qtl_pvalue=my_h5[1,'pvalues']
    qtl_logP=my_h5[1,'logP']
    
    one_side_window1=1500000
    
    qtl_ci_start=qtl_pos-one_side_window1
    qtl_ci_stop=qtl_pos+one_side_window1
    
    peaks_ci_starts=c(peaks_ci_starts,qtl_ci_start)
    peaks_ci_stops=c(peaks_ci_stops,qtl_ci_stop)		
    
    peaks_marker=c(peaks_marker,qtl_marker)
    peaks_pos=c(peaks_pos,qtl_pos)
    peaks_chr=c(peaks_chr,qtl_chr)
    peaks_pvalue=c(peaks_pvalue,qtl_pvalue)
    peaks_logP=c(peaks_logP,qtl_logP)
    
    remove=which(my_h5[,'chr']==qtl_chr & my_h5[,'pos'] >=qtl_ci_start & my_h5[,'pos']<= qtl_ci_stop)
    my_h5=my_h5[-remove,]		
    my_h5=my_h5[order(my_h5[,'logP'],decreasing=T),]
    if (dim(my_h5)[1]==0 || my_h5[1,'pvalues']>alpha) over=T
  }
  
  all_peaks=data.frame(measure=pheno_name,marker=peaks_marker,chr=peaks_chr,pos=peaks_pos,pvalue = peaks_pvalue, logP=peaks_logP,ci_starts=peaks_ci_starts,ci_stops=peaks_ci_stops,stringsAsFactors=F)
  all_peaks$merge=0
  all_new=NULL
  for (chr in unique(all_peaks[,'chr'])) {
    soub=all_peaks[which(all_peaks[,'chr']==chr),]
    changed=T
    while(changed){
      changed=F
      if (dim(soub)[1]!=1) {
        for (k in (1:(dim(soub)[1]-1))) {
          marker1=soub[k,'marker']
          for (l in ((k+1):dim(soub)[1])) {
            marker2=soub[l,'marker']
            if (soub[k,'ci_stops']>=soub[l,'ci_starts'] & soub[k,'ci_starts']<=soub[l,'ci_stops']) {
              changed=T
              if (soub[k,'merge']==0 & soub[l,'merge']==0) {
                soub[k,'merge']=k
                soub[l,'merge']=k
              } else if (soub[k,'merge']!=0 | soub[l,'merge']!=0) {
                soub[k,'merge']=soub[k,'merge']
                soub[l,'merge']=soub[k,'merge']
              }
            }
          }
        }
        
        merge_values=unique(soub[,'merge'])
        if (all(merge_values!=0)) {
          now=NULL
          for (merge_value in merge_values) {
            soubsoub=soub[soub[,'merge']==merge_value,]
            w=which.max(soubsoub[,'logP'])
            add=data.frame(measure=pheno_name,marker='merged_peak',chr=chr,pos=soubsoub[w,'pos'],pvalue = soubsoub[w,'pvalue'], logP=soubsoub[w,'logP'],ci_starts=min(soubsoub[,'ci_starts']),ci_stops=max(soubsoub[,'ci_stops']),merge=0,stringsAsFactors=F)
            now=rbind(now,add)
          }
          soub=now
        } else {
          now=soub[soub[,'merge']==0,]
          for (merge_value in merge_values[-which(merge_values==0)]) {
            soubsoub=soub[soub[,'merge']==merge_value,]
            w=which.max(soubsoub[,'logP'])
            add=data.frame(measure=pheno_name,marker='merged_peak',chr=chr,pos=soubsoub[w,'pos'],pvalue = soubsoub[w,'pvalue'], logP=soubsoub[w,'logP'],ci_starts=min(soubsoub[,'ci_starts']),ci_stops=max(soubsoub[,'ci_stops']),merge=0,stringsAsFactors=F)
            now=rbind(now,add)
          }
          soub=now
        }
      }
    }
    all_new=rbind(all_new,soub)
  }
  return(all_new)
}

res = lapply(files, my_f)
g = grep('Error', res)
print(g)

all_QTLs = do.call('rbind',res)
all_QTLs$logP=as.numeric(all_QTLs$logP)
all_QTLs = all_QTLs[order(all_QTLs$logP, decreasing = T),]
all_QTLs=all_QTLs[!grepl("Error",all_QTLs$measure),]
all_QTLs$subset_id <- sapply(strsplit(all_QTLs$measure, "_"), tail, 1)
all_QTLs$rank <- sapply(strsplit(all_QTLs$measure, "_"), head, 1)
for (subset_id in subset_ids){
  all_QTLs$measure=gsub(paste("\\_",subset_id,"$",sep=""),"",all_QTLs$measure)
}
for (i in 1:nrow(all_QTLs)){
  if (all_QTLs$rank[i]=="p"){
    all_QTLs$rank[i]="Phylum"
  } else{
    if (all_QTLs$rank[i]=="c"){
      all_QTLs$rank[i]="Class"
    } else{
      if (all_QTLs$rank[i]=="o"){
        all_QTLs$rank[i]="Order"
      } else{
        if (all_QTLs$rank[i]=="f"){
          all_QTLs$rank[i]="Family"
        } else{
          if (all_QTLs$rank[i]=="g"){
            all_QTLs$rank[i]="Genus"
          } else{
            if (all_QTLs$rank[i]=="s"){
              all_QTLs$rank[i]="Species"
            } else{
              if (all_QTLs$rank[i]=="alpha"){
                all_QTLs$rank[i]="Alpha Diversity"
              } else{
                if (all_QTLs$rank[i]=="beta"){
                  all_QTLs$rank[i]="Beta Diversity"
                } else {
                  if (all_QTLs$rank[i]=="asv"){
                    all_QTLs$rank[i]="ASV"
                  } else {
                    if (all_QTLs$rank[i]=="cluster"){
                      all_QTLs$rank[i]="Clusters"
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
all_QTLs$snp=paste0("chr",all_QTLs$chr,":",all_QTLs$pos)

### End of Step 2 ###

### Step 3 ####
print("Step 3: Make plots")

# 1. QTL Association plot
print("QTL Association plot")

## Define maximum values to be used for the axes of all plots
max_logP = 0
for (subset_id in subset_ids){
  all_QTLs_plot=all_QTLs[grepl(subset_id,all_QTLs$subset_id),]
  max_logP_sub = max(all_QTLs_plot[,'logP'], na.rm = T)
  max_logP = max(max_logP,max_logP_sub)
}

## Structure and make the plot
all_QTLs_tmp=data.frame()
pdf('gwas_qtl_association.pdf',bg='white', width = 12, height = 7)
for (subset_id in subset_ids){
  ## Calculate Padj
  sig_l=get(paste("sig_l",subset_id,sep="_"))
  QTLs=all_QTLs[grepl(subset_id,all_QTLs$subset_id),]
  QTLs <- QTLs %>% 
    mutate( Significance = case_when(logP > sig_l & subset_id == subset_id ~ "Significant", TRUE ~ "Not significant"))
  sig_snps=QTLs$snp[grepl("Significant",QTLs$Significance)]
  ## Make QTL association plot
  manhattan(x=QTLs,
            chr="chr",
            bp="pos",
            p="pvalue",
            snp="snp",
            suggestiveline =suggestive_l,
            genomewideline =sig_l,
            ylim=c(3.5,(max_logP+2)),
            main = paste("Study: ",subset_id),
            cex.main = 1.5)
  all_QTLs_tmp=rbind(all_QTLs_tmp,QTLs)
  print(subset_id)
}
dev.off()
all_QTLs=all_QTLs_tmp
all_QTLs = all_QTLs[order(all_QTLs$logP, decreasing = T),]

# Save the final data frame
save(all_QTLs,file = 'All_QTLs.RData')

# 2. Overall statistics
print("Overall statistics plot")

pdf('gwas_statistics.pdf',bg='white', width = 8, height = 6)
par(mar = c(2,2,4,2))

## 1. Per subsets
statistics=data.frame(subset=subset_ids)
### Samples
statistics$samples=all_QTLs$sample_size[match(statistics[[1]],all_QTLs$subset_id)]
### Traits
traits=count(all_QTLs,subset_id)
statistics$traits=traits$n[match(statistics[[1]],traits$subset_id)]
### Significant
sig_hits=all_QTLs[!all_QTLs$Significance=="Not significant",]
stats_sig_herit=sig_hits %>%                            
  group_by(subset_id) %>%
  summarize(min = min(logP*100), median=median(logP*100), max = max(logP*100))
sig_hits=count(sig_hits, subset_id)
statistics$sig_hits=sig_hits$n[match(statistics[[1]],sig_hits$subset_id)]
statistics$percentage=100*statistics$sig_hits/statistics$traits
statistics <- replace(statistics, is.na(statistics), 0)
assign("statistics_subsets",statistics)
### Make bar plot
statistics_plot <- tidyr::pivot_longer(statistics_subsets, cols = c(traits, sig_hits), names_to = "traits", values_to = "value")
max_value=max(statistics_plot$value)
myplot <- ggplot(statistics_plot, aes(x = subset, y = value, fill = traits)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13)) +
  labs(x = "", y = "Number of hits") +
  coord_cartesian(ylim = c(0, max_value+100)) +
  scale_fill_manual(values = c("traits" = "azure4", "sig_hits" = "red"),
                    labels = c("Significant","All hits"),
                    name = "") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  geom_text(aes(label = value), position = position_dodge(width = 0.9), vjust = -0.5)
print(myplot)

## 2. Per rank
order=c()
if ("Alpha Diversity" %in% unique(all_QTLs$rank)){
  order=c(order,"Alpha Diversity")
}
if ("Beta Diversity" %in% unique(all_QTLs$rank)){
  order=c(order,"Beta Diversity")
}
if ("Clusters" %in% unique(all_QTLs$rank)){
  order=c(order,"Clusters")
}
order=c(order,ranks)
for (subset_id in subset_ids){
  statistics=data.frame(rank=order)
  ### Traits
  filtered_VCs=all_QTLs[grepl(subset_id,all_QTLs$subset_id),]
  measure=count(filtered_VCs,rank)
  statistics$traits=measure$n[match(statistics[[1]],measure$rank)]
  ### Significant
  sig_hits=filtered_VCs[!filtered_VCs$Significance=="Not significant",]
  stats_sig_herit=sig_hits %>%                            
    group_by(rank) %>%
    summarize(min = min(logP*100), median=median(logP*100), max = max(logP*100))
  sig_hits=count(sig_hits, rank)
  statistics$sig_hits=sig_hits$n[match(statistics[[1]],sig_hits$rank)]
  statistics$percentage=100*statistics$sig_hits/statistics$traits
  statistics$min_herit=stats_sig_herit$min[match(statistics[[1]],stats_sig_herit$rank)]
  statistics$median_herit=stats_sig_herit$median[match(statistics[[1]],stats_sig_herit$rank)]
  statistics$max_herit=stats_sig_herit$max[match(statistics[[1]],stats_sig_herit$rank)]
  ### Make bar plot
  statistics_plot <- tidyr::pivot_longer(statistics, cols = c(traits, sig_hits), names_to = "traits", values_to = "value")
  statistics_plot$value[is.na(statistics_plot$value)] <- 0
  max_value=max(statistics_plot$value)
  myplot <- ggplot(statistics_plot, aes(x = factor(rank, levels = order), y = value, fill = traits)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(),
          axis.line.y = element_line(),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 13)) +
    coord_cartesian(ylim = c(0, max_value+100)) +
    labs(x= "", y = paste0("Number of loci associations (Subset: ",subset_id,")")) +
    scale_fill_manual(values = c("traits" = "azure4", "sig_hits" = "red"),
                      labels = c("Significant","All loci associations"),
                      name = "") +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    geom_text(aes(label = value), position = position_dodge(width = 0.9), vjust = -0.5)
  print(myplot)
  assign(paste("statistics_ranks",subset_id,sep="_"),statistics)
  print(subset_id)
}
dev.off()
