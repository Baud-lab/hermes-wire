#!/usr/bin/env nextflow

// Copyright (c) 2025, Centre for Genomic Regulation (CRG) and the authors.

/*
===========================================================
 HERMES-WIRE made by CRG

 @authors
 Felipe Morillo Sanz Dias <felipe.morillo.es@gmail.com>
 Luca Cozzuto <luca.cozzuto@crg.eu>
 Helene Tonnele <helene.tonnele@crg.eu>
 Amelie Baud <amelie.baud@crg.eu>
=========================================================== 
*/

/*
 Input parameters: read pairs, genome fasta file...
*/

params.help            = false
params.resume          = false


log.info """

====================================================

1. MODULE SELECTION

module analysis             : ${params.module_analysis}
module profile              : ${params.module_profile}
module comparisons          : ${params.module_comparisons}
module matrix processing    : ${params.module_matrix_processing}
module cluster              : ${params.module_cluster}
module network              : ${params.module_network}
module enterotype           : ${params.module_enterotype}
module alpha diversity      : ${params.module_alpha_diversity}
module beta diversity       : ${params.module_beta_diversity}
module heritability         : ${params.module_heritability}
module gwas                 : ${params.module_gwas}

----------------------------------------------------------------------

2.1 PROFILE: TAXONOMIC

taxonomy                    : ${params.taxonomic_file_taxonomy}
phylotree                   : ${params.taxonomic_file_phylotree}
taxonomic ranks             : ${params.taxonomic_parameter_ranks}

----------------------------------------------------------------------

2.2 PROFILE: FUNCTIONAL

functional ranks           : ${params.functional_parameter_ranks}
sample filter              : ${params.functional_tax_sample_filter}
samples and enterotypes    : ${params.functional_tax_defined_samples_and_enterotypes}

----------------------------------------------------------------------

3. MODULE: HARMONIZATION

2nd dataset table           : ${params.harmonization_file_secondary_counts_table}
2nd dataset taxonomy        : ${params.harmonization_file_secondary_taxonomy}
2nd dataset min.init depth  : ${params.harmonization_parameter_secondary_min_initial_depth}
2nd dataset min.final depth : ${params.harmonization_parameter_secondary_min_final_depth}
methods to be compared      : ${params.harmonization_parameter_methods}
taxa harmonization          : ${params.harmonization_parameter_taxa_harmonization}
hybrid matrix               : ${params.harmonization_parameter_hybrid_matrix}
proportion of main          : ${params.harmonization_parameter_hybrid_matrix_proportion_of_main}

----------------------------------------------------------------------

4. MODULE: MATRIX PROCESSING

counts table                : ${params.processing_file_counts_table}
remove ids                  : ${params.processing_file_remove_ids}
metadata                    : ${params.processing_file_metadata}
library                     : ${params.processing_file_library}
covariates                  : ${params.processing_file_covariates}
covariate types             : ${params.processing_file_cov_types}
sample identifier           : ${params.processing_parameter_sample_identifier}
host identifier             : ${params.processing_parameter_host_identifier}
minimum initial depth       : ${params.processing_parameter_min_initial_depth}
minimum abundance           : ${params.processing_parameter_min_abundance}
minimum final depth         : ${params.processing_parameter_min_final_depth}
minimum prevalence          : ${params.processing_parameter_min_prevalence}
sample material             : ${params.processing_parameter_sample_material}
method                      : ${params.processing_parameter_method}

----------------------------------------------------------------------

5. MODULE: CLUSTERS

number of neighbors         : ${params.cluster_parameter_number_of_neighbors}
minimum distance            : ${params.cluster_parameter_minimum_distance}
number of components        : ${params.cluster_parameter_number_of_components}
metric umap                 : ${params.cluster_parameter_metric_umap}
metric hdb                  : ${params.cluster_parameter_metric_hdb}
minimum cluster size        : ${params.cluster_parameter_minimum_cluster_size}
minimum samples             : ${params.cluster_parameter_minimum_samples}
epsilon                     : ${params.cluster_parameter_epsilon}
target DBCV                 : ${params.cluster_parameter_target_DBCV}
target ARI                  : ${params.cluster_parameter_target_ARI}
taxa harmonization          : ${params.cluster_parameter_taxa_harmonization}

----------------------------------------------------------------------

7. MODULE: ENTEROTYPES

enterotype rank name        : ${params.enterotypes_parameter_rank_name}
enterotype kmax             : ${params.enterotypes_parameter_kmax}
enterotype nstart           : ${params.enterotypes_parameter_nstart}

----------------------------------------------------------------------

8. MODULE: NETWORKS

network measure             : ${params.network_parameter_measure}
outliers                    : ${params.network_parameter_outliers}
threshold                   : ${params.network_parameter_threshold}
alpha                       : ${params.network_parameter_alpha}
cluster method              : ${params.network_parameter_cluster_method}
hub definition              : ${params.network_parameter_hub_definition}

----------------------------------------------------------------------

9. MODULES: ALPHA & BETA DIVERSITY

diversity type              : ${params.diversity_parameter_type}
alpha hill number           : ${params.diversity_parameter_alpha_hill_number}
alpha nboot                 : ${params.diversity_parameter_alpha_nboot}
alpha confidence            : ${params.diversity_parameter_alpha_conf}

----------------------------------------------------------------------

10. MODULES: HOST GENETIC EFFECTS ANALYSES

h5 file with traits         : ${params.genetic_file_h5}
h5 file with genetic info   : ${params.genetic_file_gen_info}
effects (Null)              : ${params.genetic_parameter_effects_null}
effects (Full)              : ${params.genetic_parameter_effects_full}

----------------------------------------------------------------------

10.1. MODULE: HERITABILITY

prev                        : ${params.heritability_file_prev}
catalogue_stats             : ${params.heritability_file_catalogue_stats}
kinship                     : ${params.heritability_parameter_kinship}
significance_fdr            : ${params.heritability_parameter_significance_fdr}
significance_bonferroni     : ${params.heritability_parameter_significance_bonferroni}

----------------------------------------------------------------------

10.2. MODULE: GWAS

residuals                   : ${params.gwas_file_residuals}
dosage_field                : ${params.gwas_parameter_dosage_field}
kinship_loco                : ${params.gwas_parameter_kinship_loco}
significance_logp           : ${params.gwas_parameter_significance_logp}
minimum logp                : ${params.gwas_parameter_min_logp}

----------------------------------------------------------------------

11. OUTPUT FOLDER

output (output folder)      : ${params.output_folder}

----------------------------------------------------------------------
"""


if (params.help) {
    log.info 'This is the MicroBiome Profiles\'s pipeline'
    log.info 'Enjoy it'
    log.info '\n'
    exit 1
}

/*
* Define output sub-folders
*/
outputMatrix_Processing     = "${params.output_folder}/Matrix_Processing"

/*
* Defining modules
*/
def local_modules_file = file("${projectDir}/local_modules.nf")
def dummy = file("${workflow.projectDir}/dummy.txt")

// Define internal parameters for the genetic analysis
cage = "all"
dam="all"
model="uni"
effectNull = "${params.genetic_parameter_effects_null}"
effectFull = "${params.genetic_parameter_effects_full},${effectNull}"
if (effectFull.contains("IGE")) {
  effectFull = "${params.genetic_parameter_effects_full},IEE,${effectNull}"
}

// Module Samples Comparison
include { samples_harmonization }  from "${local_modules_file}"
// Module Matrix Processing
include { matrix_processing_filtering; matrix_processing_beta;  matrix_processing_alpha; matrix_processing_clusters_umap; matrix_processing_clusters_sparcc; matrix_processing_clusters_network; matrix_processing_enterotypes; matrix_processing_network; matrix_processing_network_enterotypes; matrix_processing_collapse_clusters; matrix_processing_data_transformation; matrix_processing_add_beta; matrix_processing_add_alpha; matrix_processing_h5} from "${local_modules_file}"
// Modules Heritability and GWAS
include { get_row_from_matrix } from "${local_modules_file}"
include { run_exp_bivar as run_exp_bivar_full }  from "${local_modules_file}" addParams(effect: effectFull)
include { run_exp_bivar as run_exp_bivar_null }  from "${local_modules_file}" addParams(effect: effectNull)
include { concat_files as concat_full }  from "${local_modules_file}" 
include { concat_files as concat_null }  from "${local_modules_file}" 
include { mine_results;  heritability_process; heritability_cluster_analysis}  from "${local_modules_file}"
include { exp_varianceDecomp_LOCO as exp_varianceDecomp_LOCO_full }  from "${local_modules_file}" addParams(effect: effectFull)
include { exp_varianceDecomp_LOCO as exp_varianceDecomp_LOCO_null }  from "${local_modules_file}" addParams(effect: effectNull)
include { map_DGE_LOCO_GPKronSum as map_DGE_LOCO_GPKronSum_full } from "${local_modules_file}" addParams(effect: effectFull)
include { create_dir_structure } from "${local_modules_file}" addParams(effect: effectFull)
include { call_QTLs; genotypes } from "${local_modules_file}"

/*
* Getting data from parameters 
*/

//Modules
module_analysis= params.module_analysis
if (module_analysis !="Pooled" && module_analysis !="Cross-study"){
  println "WARNING: The options for 'Module Analysis' should be 'Pooled' or 'Cross-study'. By default, if none of these options is informed it is treated as 'Pooled'."
  module_analysis="Pooled"
}
module_profile= params.module_profile
if (module_profile !="Taxonomic" && module_profile !="Functional"){
  error "ERROR: Inform if the profile you are analysing is 'Taxonomic' or 'Functional'.\n"
} else {
  if (module_profile =="Taxonomic") {
    pheno="microbiome_taxonomic"
  } else {
    pheno="microbiome_functional"
  }
}
module_comparisons= params.module_comparisons
if (module_comparisons!="NO" && module_comparisons!="YES") {
  println "WARNING: The options for 'Module Harmonization' should be 'YES' or 'NO'. By default, if none of these options is informed it is treated as 'NO'."
  module_comparisons="NO"
}
module_cluster= params.module_cluster
if (module_cluster!="NO" && module_cluster!="YES") {
  println "WARNING: The options for 'Module Cluster' should be 'YES' or 'NO'. By default, if none of these options is informed it is treated as 'NO'."
  module_cluster="NO"
}
module_enterotype= params.module_enterotype
if (module_enterotype!="NO" && module_enterotype!="YES") {
  println "WARNING: The options for 'Module Enterotype' should be 'YES' or 'NO'. By default, if none of these options is informed it is treated as 'NO'."
  module_enterotype="NO"
}
module_network= params.module_network
if (module_network!="NO" && module_network!="YES") {
  println "WARNING: The options for 'Module Network' should be 'YES' or 'NO'. By default, if none of these options is informed it is treated as 'NO'."
  module_network="NO"
}
module_matrix_processing= params.module_matrix_processing
if (module_matrix_processing!="NO" && module_matrix_processing!="YES") {
  println "WARNING: The options for 'Module Matrix Processing' should be 'YES' or 'NO'. By default, if none of these options is informed it is treated as 'NO'."
  module_matrix_processing="NO"
}
if (module_matrix_processing=="NO"){
  module_alpha_diversity="NO"
  module_beta_diversity="NO"
  println "As module 'Matrix Processing' was not selected, modules 'Alpha Diversity' and 'Beta Diversity' were also disabled"
} else {
  module_alpha_diversity= params.module_alpha_diversity
  if (module_alpha_diversity!="NO" && module_alpha_diversity!="YES") {
    println "WARNING: The options for 'Module Alpha Diversity' should be 'YES' or 'NO'. By default, if none of these options is informed it is treated as 'NO'."
    module_alpha_diversity="NO"
  }
  module_beta_diversity= params.module_beta_diversity
  if (module_beta_diversity!="NO" && module_beta_diversity!="YES") {
    println "WARNING: The options for 'Module Beta Diversity' should be 'YES' or 'NO'. By default, if none of these options is informed it is treated as 'NO'."
    module_beta_diversity="NO"
  }
}
module_heritability= params.module_heritability
if (module_heritability!="NO" && module_heritability!="YES") {
  println "WARNING: The options for 'Module Heritability' should be 'YES' or 'NO'. By default, if none of these options is informed it is treated as 'NO'."
  module_heritability="NO"
}
module_gwas= params.module_gwas
if (module_gwas!="NO" && module_gwas!="YES") {
  println "WARNING: The options for 'Module GWAS' should be 'YES' or 'NO'. By default, if none of these options is informed it is treated as 'NO'."
  module_gwas="NO"
}

// Harmonization for comparisons with another dataset (i.e. 16S data)
input_matrix_harm=file(params.harmonization_file_secondary_counts_table)
taxonomy_harm = file(params.harmonization_file_secondary_taxonomy)
min_initial_depth_harm=Channel.of(params.harmonization_parameter_secondary_min_initial_depth)
min_final_depth_harm=Channel.of(params.harmonization_parameter_secondary_min_final_depth)
methods=Channel.of(params.harmonization_parameter_methods)
taxa_harmonization=params.harmonization_parameter_taxa_harmonization
hybrid=Channel.of(params.harmonization_parameter_hybrid_matrix)
proportion=Channel.of(params.harmonization_parameter_hybrid_matrix_proportion_of_main)

// Matrix processing
if (module_profile == "Taxonomic") {
  input = file(params.processing_file_counts_table) 
} else {
  file1 = file("${params.processing_file_counts_table}/EC1.out")
  file2 = file("${params.processing_file_counts_table}/EC2.out")
  file3 = file("${params.processing_file_counts_table}/EC3.out")
  file4 = file("${params.processing_file_counts_table}/EC4.out")
  input = [file1, file2, file3, file4]
  //input = [file1]
}
remove_ids = file(params.processing_file_remove_ids)
metadata = file(params.processing_file_metadata)
library = file(params.processing_file_library)
covariates = file(params.processing_file_covariates)
cov_types = file(params.processing_file_cov_types)
sample_identifier=Channel.of(params.processing_parameter_sample_identifier)
host_identifier=Channel.of(params.processing_parameter_host_identifier)
min_initial_depth=Channel.of(params.processing_parameter_min_initial_depth)
min_abundance=Channel.of(params.processing_parameter_min_abundance)
min_final_depth=Channel.of(params.processing_parameter_min_final_depth)
min_prevalence=Channel.of(params.processing_parameter_min_prevalence)
sample_material=params.processing_parameter_sample_material
if ((module_comparisons=="YES" || module_heritability=="YES" || module_gwas=="YES") && sample_material=="") {
  error "ERROR: If modules 'Comparisons', 'Heritability' or 'GWAS' were selected the 'Sample Material' should be provided.\n"
}
method=Channel.of(params.processing_parameter_method)

// Taxonomic
taxonomy = file(params.taxonomic_file_taxonomy)
phylotree = file(params.taxonomic_file_phylotree)
taxonomic_ranks=Channel.of(params.taxonomic_parameter_ranks)

// Functional
functional_ranks=Channel.of(params.functional_parameter_ranks)
tax_sample_filter=Channel.of(params.functional_tax_sample_filter)
tax_defined_samples_and_enterotypes = file(params.functional_tax_defined_samples_and_enterotypes)

// Clusters
cluster_neighbors=Channel.of(params.cluster_parameter_number_of_neighbors)
cluster_distance=Channel.of(params.cluster_parameter_minimum_distance)
cluster_components=Channel.of(params.cluster_parameter_number_of_components)
cluster_metric_umap=Channel.of(params.cluster_parameter_metric_umap)
cluster_metric_hdb=Channel.of(params.cluster_parameter_metric_hdb)
cluster_size=Channel.of(params.cluster_parameter_minimum_cluster_size)
cluster_samples=Channel.of(params.cluster_parameter_minimum_samples)
cluster_epsilon=Channel.of(params.cluster_parameter_epsilon)
cluster_target_DBCV= Channel.of(params.cluster_parameter_target_DBCV)
cluster_target_ARI=Channel.of(params.cluster_parameter_target_ARI)
cluster_harm=Channel.of(params.cluster_parameter_taxa_harmonization)

// Enterotypes
enterotype_rank=Channel.of(params.enterotypes_parameter_rank_name)
enterotype_kmax=Channel.of(params.enterotypes_parameter_kmax)
enterotype_nstart=Channel.of(params.enterotypes_parameter_nstart)

// Network
network_measure=Channel.of(params.network_parameter_measure)
outliers=Channel.of(params.network_parameter_outliers)
threshold=Channel.of(params.network_parameter_threshold)
alpha=Channel.of(params.network_parameter_alpha)
cluster_method=Channel.of(params.network_parameter_cluster_method)
hub_definition=Channel.of(params.network_parameter_hub_definition)

// Alpha and Beta Diversity
diversity_type=Channel.of(params.diversity_parameter_type)
alpha_hill_number=Channel.of(params.diversity_parameter_alpha_hill_number)
alpha_nboot=Channel.of(params.diversity_parameter_alpha_nboot)
alpha_conf=Channel.of(params.diversity_parameter_alpha_conf)

// Overall sample Genetic Effects modules
h5 = params.genetic_file_h5
if (h5==""){
  matrix=Channel.empty()
}
gene_info = file(params.genetic_file_gen_info)

// Heritability
prev = params.heritability_file_prev
catalogue_stats = file(params.heritability_file_catalogue_stats)
kinship = Channel.of(params.heritability_parameter_kinship)
significance_fdr = Channel.of(params.heritability_parameter_significance_fdr) //Check
significance_bonferroni = Channel.of(params.heritability_parameter_significance_bonferroni) //Check

// GWAS
residuals = params.gwas_file_residuals
dosage_field = Channel.of(params.gwas_parameter_dosage_field) //Check
kinship_loco = Channel.of(params.gwas_parameter_kinship_loco) //Check
significance_logp = Channel.of(params.gwas_parameter_significance_logp) //Check
min_logp = Channel.of(params.gwas_parameter_min_logp) //Check

//Empty variables
harmonized_samples=Channel.empty()
matrix_1st_filter=Channel.empty()
collapsed_ranks=Channel.empty()
data_transformation=Channel.empty()
counts_data=Channel.empty()
matrices_residues=Channel.empty()
final_counts_data=Channel.empty()
beta=Channel.empty()
alpha=Channel.empty()
input_cluster=Channel.empty()
add_beta=Channel.empty()
add_alpha=Channel.empty()
filtered_notprev=Channel.empty()
filtered_prev=Channel.empty()
final_filtered_prev=Channel.empty()
residuals=Channel.empty()
final_residuals=Channel.empty()
umap=Channel.empty()
sparcc=Channel.empty()
network=Channel.empty()
umap_path=Channel.empty()
umap_all = Channel.value(dummy)
network_path=Channel.empty()
network_all=Channel.empty()
ent_samples=Channel.empty()
heritability_results=Channel.empty()

workflow {
	// 1. Select samples and/or taxa for comparisons between two different datasets (i.e. Shallow shotgun data vs 16S)
  if (module_comparisons == "YES") {
    input_harmonized = methods.map {[
        it,
        input,
        input_matrix_harm,
        remove_ids,
        metadata,
        library,
        covariates,
        cov_types,
        module_analysis,
        params.processing_parameter_sample_identifier,
        params.processing_parameter_method,
        params.processing_parameter_min_initial_depth,
        params.processing_parameter_min_abundance,
        params.processing_parameter_min_final_depth,
        params.harmonization_parameter_secondary_min_initial_depth,
        params.harmonization_parameter_secondary_min_final_depth,
        params.processing_parameter_min_prevalence,
        params.processing_parameter_sample_material,
        taxonomy,
        taxonomy_harm,
        params.taxonomic_parameter_ranks,
        params.harmonization_parameter_hybrid_matrix,
        params.harmonization_parameter_hybrid_matrix_proportion_of_main
    ]}
    harmonized_samples = samples_harmonization(input_harmonized)
    controls=harmonized_samples.map{[it[0]]}
    filtered_notprev=harmonized_samples.map{[it[1]]}
    if (taxa_harmonization=="YES"){
          filtered_prev=harmonized_samples.map{[it[5]]}
        } else {
          filtered_prev=harmonized_samples.map{[it[4]]}
        }
  }
  
  // 2. Process the matrix
  if (module_matrix_processing == "YES") {
      if (module_comparisons == "NO"){
        // 2.1 Remove low abundance taxa and low depth samples
        input_data = method.map {
          [
              it,
              module_heritability,
              module_gwas,
              module_analysis,
              input,
              remove_ids,
              metadata,
              library,
              covariates,
              cov_types,
              params.processing_parameter_sample_identifier,
              params.processing_parameter_min_initial_depth,
              params.processing_parameter_min_abundance,
              params.processing_parameter_min_final_depth,
              params.processing_parameter_min_prevalence,
              params.processing_parameter_sample_material,
              pheno,
              taxonomy,
              params.taxonomic_parameter_ranks,
              params.functional_parameter_ranks,
              params.functional_tax_sample_filter,
              tax_defined_samples_and_enterotypes
          ]
        }
        matrix_1st_filter = matrix_processing_filtering(input_data)
        controls=matrix_1st_filter.map{[it[0]]}
        filtered_notprev=matrix_1st_filter.map{[it[1]]}
        filtered_prev=matrix_1st_filter.map{[it[2]]}
      }
    
    if (module_heritability == "YES" || module_gwas == "YES" || module_network == "YES" || module_enterotype == "YES" ||  module_cluster == "YES" || module_alpha_diversity == "YES" ||  module_beta_diversity == "YES"){
      // 2.2 Data transformation
    input_data_transformation=filtered_notprev.combine(filtered_prev).map{[
          it[0],
          it[1],
          module_comparisons,
          module_analysis,
          params.harmonization_parameter_methods,
          params.processing_parameter_method,
          metadata,
          library,
          covariates,
          cov_types,
          params.processing_parameter_sample_identifier,
          phylotree,
          taxonomy,
          params.taxonomic_parameter_ranks,
          params.functional_parameter_ranks,
          pheno
        ]}
    data_transformation=matrix_processing_data_transformation(input_data_transformation)
    residuals=data_transformation.map{[it[1]]}
    counts_data=data_transformation.map{[it[2]]}
    matrices_residues=data_transformation.map{[it[4]]}
    }
    
    if (module_cluster == "YES") {
      // 2.3 Define clusters and calculate network measurements
      input_umap=data_transformation.map{[it[3]]}.flatten().map{[
          it,
          params.cluster_parameter_number_of_neighbors,
          params.cluster_parameter_minimum_distance,
          params.cluster_parameter_number_of_components,
          params.cluster_parameter_metric_umap,
          params.cluster_parameter_metric_hdb,
          params.cluster_parameter_epsilon,
          params.cluster_parameter_minimum_cluster_size,
          params.cluster_parameter_target_DBCV,
          params.cluster_parameter_target_ARI,
          params.cluster_parameter_minimum_samples
      ]}
      //input_sparcc=matrix_1st_filter.map{[it[3]]}.flatten()
      umap=matrix_processing_clusters_umap(input_umap)
      umap_path=umap.map{[it[1]]}
      umap_all=umap_path.flatten().collect()
      //sparcc=matrix_processing_clusters_sparcc(input_sparcc)
      //input_network= umap.join(sparcc.groupTuple())
      //network=matrix_processing_clusters_network(input_network)
      //network_path=network.map{[it[1]]}
      //network_all=network_path.flatten().collect()
      input_collapse=filtered_notprev.combine(filtered_prev).combine(umap_all.map{[it]}).map{[
          it[0],
          it[1],
          it[2],
          module_comparisons,
          module_analysis,
          params.cluster_parameter_taxa_harmonization,
          params.processing_parameter_method,
          params.harmonization_parameter_methods,
          metadata,
          library,
          covariates,
          cov_types,
          phylotree,
          taxonomy,
          params.taxonomic_parameter_ranks,
          params.functional_parameter_ranks,
          pheno,
          params.processing_parameter_sample_identifier,
          params.processing_parameter_min_final_depth,
          params.processing_parameter_min_prevalence
      ]}
      clusters=matrix_processing_collapse_clusters(input_collapse)
      final_counts_data=clusters.map{it[0]}
      final_filtered_prev=clusters.map{it[1]}
      final_residuals=clusters.map{[it[3]]}
      final_filtered_notprev=clusters.map{it[7]}
    } else {
      final_filtered_prev=filtered_prev
      final_filtered_notprev=filtered_notprev
      final_residuals=residuals
      final_counts_data=counts_data
    }
    
    if (module_enterotype=="YES") {
      // 2.4 Define enterotypes
      input_enterotype=final_residuals.combine(final_filtered_prev).combine(final_counts_data).map{[
        it[0],
        it[1],
        it[2],
        module_comparisons,
        module_analysis,
        params.processing_parameter_method,
        metadata,
        library,
        covariates,
        cov_types,
        params.processing_parameter_sample_identifier,
        taxonomy,
        params.taxonomic_parameter_ranks,
        params.functional_parameter_ranks,
        pheno,
        params.enterotypes_parameter_rank_name,
        params.enterotypes_parameter_kmax,
        params.enterotypes_parameter_nstart
      ]}
      enterotype=matrix_processing_enterotypes(input_enterotype)
      ent_samples=enterotype.map{[it[0]]}
    }
    
    if (module_network=="YES") {
      // 2.5 Define networks
      input_network=final_residuals.combine(ent_samples).combine(matrices_residues).combine(umap_all.map{[it]}).map{[
        it[0],
        it[1],
        it[2],
        it[3],
        module_comparisons,
        module_cluster,
        module_analysis,
        params.processing_parameter_method,
        covariates,
        cov_types,
        taxonomy,
        params.taxonomic_parameter_ranks,
        params.functional_parameter_ranks,
        pheno,
        params.network_parameter_measure,
        params.network_parameter_outliers,
        params.network_parameter_threshold,
        params.network_parameter_alpha,
        params.network_parameter_cluster_method,
        params.network_parameter_hub_definition,
        params.functional_tax_sample_filter,
        tax_defined_samples_and_enterotypes
      ]}
      network=matrix_processing_network(input_network)
    }
    
    //if (module_network=="YES" && module_enterotype=="YES") {
    //  // 2.6 Define networks for enterotypes
    //  input_network_enterotypes=final_residuals.combine(ent_samples).map{[
    //    it[0],
    //    it[1],
    //    module_comparisons,
    //    params.processing_parameter_method,
    //    covariates,
    //    cov_types,
    //    taxonomy,
    //    params.taxonomic_parameter_ranks,
    //    params.functional_parameter_ranks,
    //    pheno,
    //    params.network_parameter_measure,
    //    params.network_parameter_cluster_method,
    //    params.network_parameter_hub_definition
    //  ]}
    //  network_enterotypes=matrix_processing_network_enterotypes(input_network_enterotypes)
    //}
    
    if (module_beta_diversity=="YES") {
      // 2.7 Calculate beta diversity
      input_beta=controls.combine(final_filtered_prev).map{[
        it[0],
        it[1],
        module_comparisons,
        module_analysis,
        params.processing_parameter_method,
        metadata,
        library,
        covariates,
        cov_types,
        params.processing_parameter_sample_identifier,
        phylotree,
        taxonomy,
        params.taxonomic_parameter_ranks,
        params.functional_parameter_ranks,
        pheno,
        params.diversity_parameter_type
      ]}
      beta=matrix_processing_beta(input_beta)
    }
    
    if (module_alpha_diversity=="YES") {
      // 2.8 Calculate alpha diversity
      input_alpha=final_filtered_prev.map{[
        it,
        module_comparisons,
        module_analysis,
        params.processing_parameter_method,
        metadata,
        library,
        covariates,
        cov_types,
        params.processing_parameter_sample_identifier,
        phylotree,
        taxonomy,
        params.taxonomic_parameter_ranks,
        params.functional_parameter_ranks,
        pheno,
        params.diversity_parameter_type,
        params.diversity_parameter_alpha_hill_number,
        params.diversity_parameter_alpha_nboot,
        params.diversity_parameter_alpha_conf
      ]}
      alpha=matrix_processing_alpha(input_alpha)
    }
    
    if ((module_heritability == "YES" || module_gwas == "YES")){
      //2.9 Add alpha and beta diversity to transformed data
      // Alpha and Beta
      if (module_beta_diversity=="YES" && module_alpha_diversity=="YES"){
        // 2.9.1 Add beta diversity to transformed data
        input_add_beta=final_counts_data.combine(beta.map{[it[0]]}).map{[
          it[0],
          it[1],
          params.processing_parameter_method,
        ]}
        add_beta=matrix_processing_add_beta(input_add_beta)
        // 2.9.2 Add alpha diversity to transformed data with beta
        input_add_alpha=add_beta.combine(alpha.map{[it[0]]}).map{[
          it[0],
          it[1],
          params.processing_parameter_method,
        ]}
        final_counts_data=matrix_processing_add_alpha(input_add_alpha)
      } else{
          // Beta
          if (module_beta_diversity=="YES" && module_alpha_diversity=="NO"){
            // 2.9.1 Add beta diversity to transformed data
            input_add_beta=final_counts_data.combine(beta.map{[it[0]]}).map{[
              it[0],
              it[1],
              params.processing_parameter_method,
            ]}
            final_counts_data=matrix_processing_add_beta(input_add_beta)
          } else {
            // Alpha
            if (module_beta_diversity=="NO" && module_alpha_diversity=="YES") {
              // 2.9.1 Add alpha diversity to transformed data without beta
              input_add_alpha=final_counts_data.combine(alpha.map{[it[0]]}).map{[
                it[0],
                it[1],
                params.processing_parameter_method,
              ]}
              final_counts_data=matrix_processing_add_alpha(input_add_alpha)
            }
          }
      }
      // 2.10 Make H5 file
      input_h5 = final_counts_data.combine(final_residuals).map {
          [
              it[0],
              it[1],
              params.heritability_parameter_kinship,
              params.gwas_parameter_kinship_loco,
              module_heritability,
              module_gwas,
              module_analysis,
              metadata,
              covariates,
              cov_types,
              params.processing_parameter_sample_identifier,
              params.processing_parameter_host_identifier,
              params.processing_parameter_method,
              params.taxonomic_parameter_ranks,
              params.functional_parameter_ranks,
              pheno,
              cage,
              dam,
              gene_info,
              params.gwas_parameter_dosage_field
          ]
      }
      matrix=matrix_processing_h5(input_h5)
    }
  } else {
      if ((module_heritability == "YES" || module_gwas == "YES") && h5==""){
        error "ERROR: Inform the path to a H5 file for Host Genetic Analises if you don't intend to run the 'Matrix Processing' step.\n"
      } else {
        if (module_heritability == "YES" || module_gwas == "YES"){
         matrix = Channel.fromPath(h5) 
        }
      }
      if ( module_heritability =="YES" && prev==""){
        error "ERROR: Inform the path to a RData file with the matrices after filtering traits by prevalence if you don't intend to run the 'Matrix Processing' step.\n"
      } else {
        final_filtered_prev = Channel.fromPath(prev)
      }      
      if ( module_gwas =="YES" && residuals==""){
        error "ERROR: Inform the path to a RData file with the matrices after CLR transformation if you don't intend to run the 'Matrix Processing' step.\n"
      } else {
        final_residuals = Channel.fromPath(params.gwas_file_residuals)
      }      
    }
  
  if (module_heritability == "YES" || module_gwas == "YES") {
    
    // 3. Get the number of rows from the filtered matrix
    max_rows_ch = get_row_from_matrix(Channel.of(pheno).combine(matrix)).first().toInteger()
    rows = max_rows_ch.flatMap { value -> (1..value) }
    
    if (module_heritability == "YES") {
      
      // 4. Calculate heritability for both the null (with no DGE) and full (with DGE) models
      data = rows.combine(matrix).map{ [ it[0], it[1], pheno, cage, dam, params.heritability_parameter_kinship ] }
      full_res = run_exp_bivar_full(data)
      null_res = run_exp_bivar_null(data)
      cat_full = concat_full("output_full.txt", full_res.collect()) // Check the new scripts
      cat_null = concat_null("output_null.txt", null_res.collect()) // Check the new scripts
      cat_files = cat_full.map {
          [effectFull, params.heritability_parameter_kinship, model, pheno, it]
      }.mix(cat_null.map {
          [effectNull, params.heritability_parameter_kinship, model, pheno, it]
      })
      mined_res = mine_results(cat_files)
      mined_res.branch { effect, rdata ->
          efull: effect == effectFull
          enull: effect == effectNull
      }.set { h5files }
      
      // 5. Calculate p-values and define significance for Heritability
      pval_data = kinship.combine(Channel.of(model)).combine(Channel.of(pheno)).combine(h5files.enull).combine(h5files.efull).combine(Channel.of(module_analysis)).combine(Channel.of(module_comparisons)).combine(Channel.of(module_cluster)).combine(method).combine(methods).combine(significance_fdr).combine(significance_bonferroni).combine(min_prevalence).combine(taxonomic_ranks).combine(functional_ranks).combine(final_filtered_prev).combine(Channel.of([catalogue_stats,taxonomy, phylotree, covariates]))
      heritability_results=heritability_process(pval_data)
      //if (module_cluster=="YES"){
      //  // 6. Analyse heritable traits from the point of view of network analyses
      //  cluster_herit_input = heritability_results.map{[it[0]]}.combine(network_all.map{[it]}).map{[
      //    it[0],
      //    it[1],
      //    params.processing_parameter_method,
      //    metadata,
      //    covariates,
      //    params.processing_parameter_sample_identifier,
      //    taxonomy,
      //    params.taxonomic_parameter_ranks
      //  ]}
      //  heritability_cluster_analysis(cluster_herit_input) 
      //}
      
    }
    
    if (module_gwas == "YES") {
      
      // 7. Process GWAS
      data2 = rows.combine(matrix).map{ [ it[0], it[1], pheno, cage, dam, params.gwas_parameter_kinship_loco ] }
      full_res2 = exp_varianceDecomp_LOCO_full(data2).collect()
      data_folder = full_res2.map{ [pheno, params.gwas_parameter_kinship_loco, it ] }
      fold_struct = create_dir_structure(data_folder)
      data_map = data2.combine(fold_struct)
      outdir = map_DGE_LOCO_GPKronSum_full(data_map)
      dir_for_call = outdir.collect().flatten().first()
      data_for_call = dir_for_call.map{
      	[  effectFull,
      	params.gwas_parameter_kinship_loco,
      	model,
      	pheno,
      	params.gwas_parameter_significance_logp,
      	params.gwas_parameter_min_logp,
      	it ] }.combine(Channel.of([covariates])).combine(taxonomic_ranks).combine(functional_ranks).combine(final_residuals) 
      
      // 8. Call QTLs and calculate significance for loci association
      outDir = call_QTLs(data_for_call).collect().flatten().first()
      
      //// 9. Make boxplots of relative abundances by genotypes for significant associations between traits and SNPs
      //geno_input = outDir.map{[it[0]]}.combine(final_residuals).combine(matrix).map{[
      //  it[0],
      //  it[1],
      //  it[2],
      //  module_cluster,
      //  metadata,
      //  covariates,
      //  params.processing_parameter_sample_identifier,
      //  params.processing_parameter_host_identifier,
      //  params.processing_parameter_method,
      //  params.taxonomic_parameter_ranks
      //]}
      //geno_input.view()
      //genotypes_out = genotypes(geno_input)
    }
  }
  
}


/*
* When is finished
*/
workflow.onComplete {
    println "Pipeline MicroBiome Profiling is completed!"
    println "Started at  $workflow.start" 
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'FAILED' }"
}
