/*
* Define output sub-folders
*/
outputHarmonization     = "${params.output_folder}/Harmonization"
outputMatrix_Processing     = "${params.output_folder}/Matrix_Processing"
outputCluster_Analyses     = "${params.output_folder}/Cluster_Analyses"
outputEnterotypes     = "${params.output_folder}/Enterotypes"
outputDiversity     = "${params.output_folder}/Diversity"
outputNetwork     = "${params.output_folder}/Network"
outputHerit   = "${params.output_folder}/Heritability"
outputGWAS   = "${params.output_folder}/GWAS"

// Define internal parameter for the effects evaluated on the heritability analysis
params.effect = ""

// Processes

// 1. Select samples for comparisons between two different datasets (i.e. Shallow shotgun data vs 16S)
process samples_harmonization {

    tag ("Sample harmonization: ${methods}")
    publishDir(outputHarmonization, mode:'copy') 
    memory='50G'
    cpus='1'
    time='0.5h'
	
    container 'morillofe/microbiomeprofile:0.7'
    
    input:
    tuple val(methods),
    path(input),
    path(input_harm),
    path(remove_ids),
    path(metadata),
    path(library),
    path(covariates),
    path(cov_types),
    val(module_analysis),
    val(sample_identifier),
    val(method),
    val(min_initial_depth),
    val(min_abundance),
    val(min_final_depth),
    val(min_initial_depth_harm),
    val(min_final_depth_harm),
    val(min_prevalence),
    val(sample_material),
    path(taxonomy),
    path(taxonomy_harm),
    val(taxonomic_ranks),
    val(hybrid),
    val(proportion)
    
    output:
    tuple path("controls.RData"),
    path("filtered_traits_matrix_noprev.RData"),
    path("alpha_table_main.RData"),
    path("alpha_table_secondary.RData"),
    path("harmonized_samples.RData"),
    path("harmonized_samples_and_taxa.RData"),
    path("depth_funnel_focal_method.pdf"),
    path("asymptotic_richness_vs_depth.pdf"),
    path("phylum_stacked_barplots.pdf")
    
	
    script:
    """
    samples_harmonization.R \
            -mets ${methods} \
            -met ${method} \
            -input ${input} \
            -input_harm ${input_harm} \
            -exclude ${remove_ids} \
            -meta ${metadata} \
            -lib ${library} \
            -covariates ${covariates} \
            -cov_types ${cov_types} \
            -mod_ana ${module_analysis} \
            -id ${sample_identifier} \
            -mid ${min_initial_depth} \
            -ma ${min_abundance} \
            -mfd ${min_final_depth} \
            -mid_harm ${min_initial_depth_harm} \
            -mfd_harm ${min_final_depth_harm} \
            -mp ${min_prevalence} \
            -mat ${sample_material} \
            -tax ${taxonomy} \
            -taxo_harm ${taxonomy_harm} \
            -ranks ${taxonomic_ranks} \
            -hyb ${hybrid} \
            -prop ${proportion}
    """
}

// 2. Matrix processing - Filtering low abundance and prevalence bacteria and low depth samples
process matrix_processing_filtering {

    tag ("Filter taxa and samples: ${method}")
    publishDir(outputMatrix_Processing, mode:'copy') 
    memory='10G'
    cpus='1'
    time='0.5h'
	
    container 'morillofe/microbiomeprofile:0.7'
    
    input:
    tuple val(method),
    val(module_heritability),
    val(module_gwas),
    val(module_analysis),
    path(input),
    path(remove_ids),
    path(metadata),
    path(library),
    path(covariates),
    path(cov_types),
    val(sample_identifier),
    val(min_initial_depth),
    val(min_abundance),
    val(min_final_depth),
    val(min_prevalence),
    val(sample_material),
    val(pheno),
    path(taxonomy),
    val(taxonomic_ranks),
    val(functional_ranks),
    val(tax_sample_filter),
    path(tax_defined_samples_and_enterotypes)

    output:
    tuple path("controls.RData"),
    path("filtered_traits_matrix_noprev.RData"),
    path("filtered_prev.RData"),
    path("matrices_abundance_filtered/*.csv"),
    path("*.pdf")
	
    script:
    """
    if [ ${pheno} = "microbiome_taxonomic" ]; then
        matrix_processing_filtering.R \
            -met ${method} \
            -mod_herit ${module_heritability} \
            -mod_gwas ${module_gwas} \
            -mod_ana ${module_analysis} \
            -input ${input} \
            -exclude ${remove_ids} \
            -meta ${metadata} \
            -lib ${library} \
            -covariates ${covariates} \
            -cov_types ${cov_types} \
            -id ${sample_identifier} \
            -mid ${min_initial_depth} \
            -ma ${min_abundance} \
            -mfd ${min_final_depth} \
            -mp ${min_prevalence} \
            -mat ${sample_material} \
            -tax ${taxonomy} \
            -ranks  ${taxonomic_ranks}
    else
        matrix_processing_filtering_functional.R \
            -met ${method} \
            -mod_herit ${module_heritability} \
            -mod_gwas ${module_gwas} \
            -mod_ana ${module_analysis} \
            -ec1 ${input[0]} \
            -ec2 ${input[1]} \
            -ec3 ${input[2]} \
            -ec4 ${input[3]} \
            -exclude ${remove_ids} \
            -meta ${metadata} \
            -lib ${library} \
            -covariates ${covariates} \
            -cov_types ${cov_types} \
            -id ${sample_identifier} \
            -mid ${min_initial_depth} \
            -ma ${min_abundance} \
            -mfd ${min_final_depth} \
            -mp ${min_prevalence} \
            -mat ${sample_material} \
            -ranks  ${functional_ranks} \
            -samp_filt ${tax_sample_filter} \
            -samps ${tax_defined_samples_and_enterotypes}
    fi
    mkdir matrices_abundance_filtered
    mv *.csv matrices_abundance_filtered/
    """
}

process matrix_processing_data_transformation {

    tag ("Data transformation: ${method}")
    publishDir(outputMatrix_Processing, mode:'copy') 
    memory='1.5G'
    cpus='1'
    time='0.5h'
	
    container 'morillofe/microbiomeprofile:0.7'
    
    input:
    tuple path(filtered_notprev),
    path(filtered_prev),
    val(module_comparisons),
    val(module_analysis),
    val(methods),
    val(method),
    path(metadata),
    path(library),
    path(covariates),
    path(cov_types),
    val(sample_identifier),
    path(phylotree),
    path(taxonomy),
    val(taxonomic_ranks),
    val(functional_ranks),
    val(pheno)

    output:
    tuple path("filtered_clr_counts.RData"),
    path("residuals_qned_counts.RData"),
    path("count_data.RData"),
    path("matrices_residues/*.csv"),
    path("matrices_residues")
	
    script:
    """
   matrix_processing_data_transformation.R \
        -mets ${methods} \
        -met ${method} \
        -phe ${pheno} \
        -mod_comp ${module_comparisons} \
        -mod_ana ${module_analysis} \
        -prev ${filtered_prev} \
        -not_prev ${filtered_notprev} \
        -meta ${metadata} \
        -lib ${library} \
        -covariates ${covariates} \
        -cov_types ${cov_types} \
        -id ${sample_identifier} \
        -phy ${phylotree} \
        -tax ${taxonomy} \
        -tax_ranks  ${taxonomic_ranks} \
        -func_ranks  ${functional_ranks}
    mkdir matrices_residues
    mv *.csv matrices_residues/
    """
}

process matrix_processing_clusters_umap {

    tag ("Define clusters (UMAP-HDBSCAN): ${input_cluster.baseName}")
    publishDir(outputCluster_Analyses, mode:'copy') 
    memory='1G'
    cpus='1'
    time='0.5h'
	
    container 'morillofe/micnet:1.0.12'
    
    input:
    tuple path(input_cluster),
    val(cluster_neighbors),
    val(cluster_distance),
    val(cluster_components),
    val(cluster_metric_umap),
    val(cluster_metric_hdb),
    val(cluster_epsilon),
    val(cluster_size),
    val(target_DBCV),
    val(target_ARI),
    val(cluster_samples)

    output:
    tuple val("${input_cluster.baseName}"),
    path("Cluster_UMAP_HDBSCAN_${input_cluster.baseName}.csv"),
    path("UMAP_HDBSCAN_Clusters_${input_cluster.baseName}.pdf")
	
    script:
    def subset = input_cluster.baseName
    """
    micnet.py --tool umap_hdbscan --file_path ${input_cluster} \
    --n_neighbors ${cluster_neighbors} \
    --min_dist ${cluster_distance} \
    --n_components ${cluster_components} \
    --metric_umap ${cluster_metric_umap} \
    --metric_hdb ${cluster_metric_hdb} \
    --cluster_epsilon ${cluster_epsilon} \
    --min_cluster_size ${cluster_size} \
    --target_dbcv ${target_DBCV} \
    --target_ari ${target_ARI} \
    --min_sample ${cluster_samples}
    mv Cluster_UMAP_HDBSCAN.csv Cluster_UMAP_HDBSCAN_${subset}.csv
    mv UMAP_HDBSCAN_Clusters.pdf UMAP_HDBSCAN_Clusters_${subset}.pdf
    """
}

process matrix_processing_clusters_sparcc {

    tag ("Define clusters (SPARCC): ${input_cluster.baseName}")
    publishDir(outputCluster_Analyses, mode:'copy') 
    memory='10G'
    cpus='1'
    time='1.5h'
	
    container 'morillofe/micnet:1.0.12'
    
    input:
    path(input_cluster)

    output:
    tuple val("${input_cluster.baseName}"), path("Cluster_SparCC_${input_cluster.baseName}.csv")
	
    script:
    def subset = input_cluster.baseName
    """
    micnet.py --tool sparcc --file_path ${input_cluster} \
    --n_iteractions 20 \
    --x_iteractions 10 \
    --threshold 0.1 \
    --num_simulate_data 5
    mv Cluster_SparCC.csv Cluster_SparCC_${subset}.csv
    """
}

process matrix_processing_clusters_network {

    tag ("Define clusters (Network): ${subset}")
    publishDir(outputCluster_Analyses, mode:'copy') 
    memory='10G'
    cpus='1'
    time='6h'
	
    container 'morillofe/micnet:1.0.12'
    
    input:
    tuple val(subset), path(umap), path(sparcc) 

    output:
    tuple val("${subset}"), path("Cluster_Centralities_${subset}.csv"),
    path("Network_${subset}")
	
    script:
    """
    mkdir Network_${subset}
    micnet.py --tool network --sparcc_path ${sparcc} --umap_path ${umap} \
    --nboot_top 20 \
    --bins_top 20 \
    --prem_perc 0.1
    mv Cluster_Centralities.csv Cluster_Centralities_${subset}.csv
    mv *.txt *.pdf Network_${subset}
    """
}

process matrix_processing_collapse_clusters {

    tag ("Collapse clusters: ${method}")
    publishDir(outputCluster_Analyses, mode:'copy') 
    memory='1.5G'
    cpus='1'
    time='0.5h'
	
    container 'morillofe/microbiomeprofile:0.7'
    
    input:
    tuple path(filtered_notprev),
    path(filtered_prev),
    path(umap_all),
    val(module_comparisons),
    val(module_analysis),
    val(cluster_harm),
    val(method),
    val(methods),
    path(metadata),
    path(library),
    path(covariates),
    path(cov_types),
    path(phylotree),
    path(taxonomy),
    val(taxonomic_ranks),
    val(functional_ranks),
    val(pheno),
    val(sample_identifier),
    val(min_final_depth),
    val(min_prevalence)

    output:
    tuple path("count_data_new.RData"),
    path("filtered_prev_new.RData"),
    path("filtered_clr_counts_new.RData"),
    path("residuals_qned_counts_new.RData"),
    path("filtered_prev_clusters.RData"),
    path("filtered_clr_counts_clusters.RData"),
    path("residuals_qned_counts_clusters.RData"),
    path("filtered_traits_matrix_noprev_new.RData"),
    path("cluster_stacked_barplots.pdf"),
    path("stacked_barplots_new.pdf")
	
    script:
    """
    mkdir umap_dir; cp ${umap_all} umap_dir/
    matrix_processing_collapse_clusters.R \
        -prev ${filtered_prev} \
        -not_prev ${filtered_notprev} \
        -umap umap_dir \
        -mod_comp ${module_comparisons} \
        -mod_ana ${module_analysis} \
        -harm ${cluster_harm} \
        -met ${method} \
        -phe ${pheno} \
        -mets ${methods} \
        -meta ${metadata} \
        -lib ${library} \
        -covariates ${covariates} \
        -cov_types ${cov_types} \
        -phy ${phylotree} \
        -tax ${taxonomy} \
        -tax_ranks  ${taxonomic_ranks} \
        -func_ranks  ${functional_ranks} \
        -id ${sample_identifier} \
        -mfd ${min_final_depth} \
        -mp ${min_prevalence}
    """
}

process matrix_processing_enterotypes {

    tag ("Enterotypes: ${method}")
    publishDir(outputEnterotypes, mode:'copy') 
    memory='1.5G'
    cpus='1'
    time='0.5h'
	
    container 'morillofe/microbiomeprofile:0.7'
    
    input:
    tuple path(final_residuals),
    path(final_prev),
    path(previous_counts),
    val(module_comparisons),
    val(module_analysis),
    val(method),
    path(metadata),
    path(library),
    path(covariates),
    path(cov_types),
    val(sample_identifier),
    path(taxonomy),
    val(taxonomic_ranks),
    val(functional_ranks),
    val(pheno),
    val(ent_rank),
    val(ent_kmax),
    val(ent_nstart)
    

    output:
    tuple path("enterotypes.RData"),
    path("enterotypes.pdf")
	
    script:
    """
    matrix_processing_enterotypes.R \
        -met ${method} \
        -phe ${pheno} \
        -mod_comp ${module_comparisons} \
        -mod_ana ${module_analysis} \
        -input ${final_prev} \
        -residues ${final_residuals} \
        -counts ${previous_counts} \
        -meta ${metadata} \
        -lib ${library} \
        -covariates ${covariates} \
        -cov_types ${cov_types} \
        -id ${sample_identifier} \
        -tax ${taxonomy} \
        -tax_ranks ${taxonomic_ranks} \
        -func_ranks ${functional_ranks} \
        -ent_rank ${ent_rank} \
        -ent_kmax ${ent_kmax} \
        -ent_nstart ${ent_nstart}
    """
}

process matrix_processing_network {

    tag ("Network: ${method}")
    publishDir(outputNetwork, mode:'copy') 
    memory='10G'
    cpus='1'
    time='48h'
	
    container 'morillofe/netcomi:1.3'
    
    input:
    tuple path(final_residuals),
    path(ent_samples),
    path(matrices_residues),
    path(umap_all),
    val(module_comparisons),
    val(module_cluster),
    val(module_analysis),
    val(method),
    path(covariates),
    path(cov_types),
    path(taxonomy),
    val(taxonomic_ranks),
    val(functional_ranks),
    val(pheno),
    val(network_measure),
    val(network_outliers),
    val(network_threshold),
    val(network_alpha),
    val(cluster_method),
    val(hub_definition),
    val(tax_sample_filter),
    path(tax_defined_samples_and_enterotypes)
    

    output:
    tuple path("*.pdf"),
    path("*.RData")
	
    script:
    """
     mkdir umap_dir; cp ${umap_all} umap_dir/
    matrix_processing_network.R \
        -met ${method} \
        -phe ${pheno} \
        -umap umap_dir \
        -umap_res ${matrices_residues} \
        -mod_comp ${module_comparisons} \
        -mod_clust ${module_cluster} \
        -mod_ana ${module_analysis} \
        -residues ${final_residuals} \
        -ent ${ent_samples} \
        -covariates ${covariates} \
        -cov_types ${cov_types} \
        -tax ${taxonomy} \
        -tax_ranks ${taxonomic_ranks} \
        -func_ranks ${functional_ranks} \
        -measure ${network_measure} \
        -outl ${network_outliers} \
        -thresh ${network_threshold} \
        -alpha ${network_alpha} \
        -cluster ${cluster_method} \
        -hub ${hub_definition} \
        -samp_filt ${tax_sample_filter} \
        -samps ${tax_defined_samples_and_enterotypes}
    """
}

process matrix_processing_network_enterotypes {

    tag ("Network: ${method}")
    publishDir(outputNetwork, mode:'copy') 
    memory='10G'
    cpus='1'
    time='48h'
	
    container 'morillofe/netcomi:1.2'
    
    input:
    tuple path(final_residuals),
    path(enterotypes),
    val(module_comparisons),
    val(module_analysis),
    val(method),
    path(covariates),
    path(cov_types),
    path(taxonomy),
    val(taxonomic_ranks),
    val(functional_ranks),
    val(pheno),
    val(network_measure),
    val(cluster_method),
    val(hub_definition)
    

    output:
    tuple path("networks_enterotypes.pdf"),
    path("networks_enterotypes.RData")
	
    script:
    """
    matrix_processing_network_enterotypes.R \
        -met ${method} \
        -mod_comp ${module_comparisons} \
        -mod_ana ${module_analysis} \
        -residues ${final_residuals} \
        -covariates ${covariates} \
        -cov_types ${cov_types} \
        -ent ${enterotypes} \
        -tax ${taxonomy} \
        -ranks ${taxonomic_ranks} \
        -measure ${network_measure} \
        -cluster ${cluster_method} \
        -hub ${hub_definition}
    """
}

process matrix_processing_beta {

    tag ("Beta diversity: ${method}")
    publishDir(outputDiversity, mode:'copy') 
    memory='3G'
    cpus='1'
    time='0.12h'
	
    container 'morillofe/microbiomeprofile:0.7'
    
    input:
    tuple path(controls),
    path(filtered_prev),
    val(module_comparisons),
    val(module_analysis),
    val(method),
    path(metadata),
    path(library),
    path(covariates),
    path(cov_types),
    val(sample_identifier),
    path(phylotree),
    path(taxonomy),
    val(taxonomic_ranks),
    val(functional_ranks),
    val(pheno),
    val(diversity_type)
    

    output:
    tuple path("count_data_beta.RData"),
    path("residuals_qned_counts_beta.RData"),
    path("beta_diversity.RData"),
    path("permanova_beta_diversity.RData"),
    path("beta_diversity.pdf")
	
    script:
    """
    matrix_processing_beta.R \
        -met ${method} \
        -phe ${pheno} \
        -mod_comp ${module_comparisons} \
        -mod_ana ${module_analysis} \
        -cont ${controls} \
        -input ${filtered_prev} \
        -meta ${metadata} \
        -lib ${library} \
        -covariates ${covariates} \
        -cov_types ${cov_types} \
        -id ${sample_identifier} \
        -phy ${phylotree} \
        -tax ${taxonomy} \
        -tax_ranks ${taxonomic_ranks} \
        -func_ranks ${functional_ranks} \
        -divt  ${diversity_type} \
    """
}

process matrix_processing_alpha {

    tag ("Alpha diversity: ${method}")
    publishDir(outputDiversity, mode:'copy') 
    memory='10G'
    cpus='1'
    time='48h'
	
    container 'morillofe/microbiomeprofile:0.7'
    
    input:
    tuple path(filtered_prev),
    val(module_comparisons),
    val(module_analysis),
    val(method),
    path(metadata),
    path(library),
    path(covariates),
    path(cov_types),
    val(sample_identifier),
    path(phylotree),
    path(taxonomy),
    val(taxonomic_ranks),
    val(functional_ranks),
    val(pheno),
    val(diversity_type),
    val(alpha_hill_number),
    val(alpha_nboot),
    val(alpha_conf)

    output:
    tuple path("count_data_alpha.RData"),
    path("residuals_qned_counts_alpha.RData"),
    path("alpha_diversity.RData"),
    path("anova_alpha_diversity.RData"),
    path("alpha_diversity.pdf")
	
    script:
    """
    matrix_processing_alpha.R \
        -met ${method} \
        -phe ${pheno} \
        -input ${filtered_prev} \
        -meta ${metadata} \
        -lib ${library} \
        -mod_comp ${module_comparisons} \
        -mod_ana ${module_analysis} \
        -covariates ${covariates} \
        -cov_types ${cov_types} \
        -id ${sample_identifier} \
        -phy ${phylotree} \
        -tax ${taxonomy} \
        -tax_ranks ${taxonomic_ranks} \
        -func_ranks ${functional_ranks} \
        -divt  ${diversity_type} \
        -alpha_h  ${alpha_hill_number} \
        -alpha_n  ${alpha_nboot} \
        -alpha_c  ${alpha_conf} \
    """
}

process matrix_processing_add_beta {

    tag ("Add beta diversity values: ${method}")
    publishDir(outputDiversity, mode:'copy') 
    memory='1G'
    cpus='1'
    time='0.1h'
	
    container 'morillofe/microbiomeprofile:0.7'
    
    input:
    tuple path(input),
    path(beta),
    val(method)

    output:
    path("count_data_with_beta.RData")
	
    script:
    """
   matrix_processing_add_beta.R \
        -met ${method} \
        -input ${input} \
        -b ${beta} \
    """
}

process matrix_processing_add_alpha {

    tag ("Add alpha diversity values: ${method}")
    publishDir(outputDiversity, mode:'copy') 
    memory='1G'
    cpus='1'
    time='0.1h'
	
    container 'morillofe/microbiomeprofile:0.7'
    
    input:
    tuple path(input),
    path(alpha),
    val(method)

    output:
    path("count_data_with_alpha.RData")
	
    script:
    """
   matrix_processing_add_alpha.R \
        -met ${method} \
        -input ${input} \
        -a ${alpha} \
    """
}

process matrix_processing_h5 {

    tag ("H5 file: ${kinship}")
    publishDir(outputMatrix_Processing, mode:'copy') 
    memory='5G'
    cpus='1'
    time='0.2h'
	
    container 'morillofe/microbiomeprofile:0.7'
    
    input:
    tuple path(counts_data),
    path(final_residuals),
    val(kinship),
    val(kloco),
    val(module_heritability),
    val(module_gwas),
    val(module_analysis),
    path(metadata),
    path(covariates),
    path(cov_types),
    val(sample_identifier),
    val(host_identifier),
    val(method),
    val(taxonomic_ranks),
    val(functional_ranks),
    val(pheno),
    val(cage),
    val(dam),
    path(gene_info),
    val(dosage_field)
	
    output:
    path("${kinship}.h5")
	
    script:
    """
    cp ${gene_info} ${kinship}.h5
    matrix_processing_h5.R \
        -mod_herit ${module_heritability} \
        -mod_gwas ${module_gwas} \
        -mod_ana ${module_analysis} \
        -input ${counts_data} \
        -meta ${metadata} \
        -covariates ${covariates} \
        -cov_types ${cov_types} \
        -res ${final_residuals} \
        -sid ${sample_identifier} \
        -hid ${host_identifier} \
        -met ${method} \
        -pheno ${pheno} \
        -tax_ranks ${taxonomic_ranks} \
        -func_ranks ${functional_ranks} \
        -cage ${cage} \
        -dam ${dam} \
        -kin ${kinship} \
        -dosf ${dosage_field} \
        -kinloco ${kloco} \
        -out ${kinship}.h5
    """
}

// 3. Get the number of rows from the filtered matrix
process get_row_from_matrix {

	tag ("Get rows from matrix: ${matrix}")
	memory='1G'
  cpus='1'
  time='0.05h'
	
	container 'morillofe/microbiomeprofile:0.7'
    
    input:
    tuple val(pheno), path(matrix)	
	
    output:
	stdout
	
	script:
	"""
	if [ ! -f ${matrix} ]; then
        error "Error: H5 file with microbiome traits not found. Please, select a H5 file if you are not using module 'Matrix Processing'"
  else
        get_row_from_mat.R -h5 ${matrix} -pheno ${pheno}
  fi
	
	"""

}

// 4. Calculate heritability for both the null (with no DGE) and full (with DGE) models
process run_exp_bivar {

    tag ("Feature number: ${row_num}")	
    container 'biocorecrg/microbiomepython:0.1'
    label('ignore')
    //cache true
    
    input:
    tuple val(row_num), file(matrix), val(pheno), val(cage), val(dam), val(kinship)	

    output:
	  path("*.txt")
	
 	
    script:
    """
    exp_bivar12.py \
  	${matrix} \
  	${pheno} \
  	None \
  	${cage} \
  	${dam} \
  	${kinship} \
  	--out ./ \
  	-p ${row_num} \
  	-e ${params.effect}	
  	cp ./*/*/*/*.txt .
  	"""

}

process concat_files {

    memory='5G'
    cpus='1'
    time='0.5h'
	
    input:
	  val(outfile)
	  path(input_files)
	
    output:
	  path(outfile)
		
    script:

    """
	  cat ${input_files} >> ${outfile}
    """

}


process mine_results {

    tag ("Aggregate results: ${effect}")
		memory='1G'
    cpus='1'
    time='0.2h'
    
    container 'morillofe/microbiomeprofile:0.7'
    
    input:
    tuple val(effect), val(kinship), val(model), val(pheno), path(catfile)
	
    output:
	  tuple val(effect), path("*.RData")
	
    script:
    """
		mine_results.R  \
		-input ${catfile} \
		-kin ${kinship} \
		-model ${model} \
		-pheno ${pheno} \
		-eff ${effect}
    """

}

// 5. Calculate p-values and define significance for Heritability
process heritability_process {
    
    memory='1G'
    cpus='1'
    time='0.5h'    
    publishDir(outputHerit, mode:'copy')
    container 'morillofe/microbiomeprofile:0.7'

    input:
    tuple val(kin),
    val(model),
    val(pheno),
    val(effN),
    path(inputN),
    val(effF),
    path(inputF),
    val(module_analysis),
    val(module_comparisons),
    val(module_cluster),
    val(method),
    val(methods),
    val(significance_fdr),
    val(significance_bonferroni),
    val(min_prevalence),
    val(taxonomic_ranks),
    val(functional_ranks),
    path(filtered_prev),
    path(catalogue_stats),
    path(taxonomy),
    path(phylotree),
    path(covariates)

    output:
    tuple path("heritability.RData"),
    path("common_traits.RData"),
    //path("anova_rank_effect_on_heritability.RData"),
    //path("statistics_per_rank.RData"),
    path("*.pdf")

    script:
    """
    heritability.R \
    	  -inputN ${inputN} \
    	  -inputF ${inputF} \
    	  -kin ${kin} \
    	  -model ${model} \
    	  -pheno ${pheno} \
    	  -effN ${effN} \
    	  -effF ${effF} \
    	  -mod_ana ${module_analysis} \
    	  -mod_comp ${module_comparisons} \
    	  -mod_clust ${module_cluster} \
    	  -met ${method} \
    	  -mets ${methods} \
    	  -sig_f ${significance_fdr} \
    	  -sig_b ${significance_bonferroni} \
    	  -mp ${min_prevalence} \
    	  -tax_ranks ${taxonomic_ranks} \
    	  -func_ranks ${functional_ranks} \
    	  -prev ${filtered_prev} \
    	  -stats ${catalogue_stats} \
    	  -tax ${taxonomy} \
    	  -phy ${phylotree} \
    	  -covariates ${covariates}
    """

}

// 6. Analyse heritable traits from the point of view of network analyses
process heritability_cluster_analysis {
    
    memory='1G'
    cpus='1'
    time='0.5h'    
    publishDir(outputHerit, mode:'copy')

    container 'morillofe/microbiomeprofile:0.7'

    input:
    tuple path(heritability),
    path(network),
    val(method),
    path(metadata),
    path(covariates),
    val(sample_identifier),
    path(taxonomy),
    val(taxonomic_ranks)
    
    output:
    path("*.pdf")

    script:
    """
    mkdir network_dir; mv ${network} network_dir/
    heritability_cluster_analysis.R \
    	  -herit ${heritability} \
    	  -net network_dir \
    	  -met ${method} \
    	  -meta ${metadata} \
    	  -id ${sample_identifier} \
    	  -ranks ${taxonomic_ranks} \
    	  -tax ${taxonomy} \
    	  -covariates ${covariates}
    """

}

// 7. Process GWAS analyses
process exp_varianceDecomp_LOCO {

    tag ("Feature number: ${row_num}")	
    container 'biocorecrg/microbiomepython:0.2'
    label('ignore')
    //cache true

    input:
    tuple val(row_num), file(matrix), val(pheno), val(cage), val(dam), val(grm_v)	

    output:
    path("*__*")
	
    script:
    """
  	exp_varianceDecomp_LOCO.py \
  	${matrix} \
  	${pheno} \
  	None \
  	${cage} \
  	${grm_v} \
  	--out ./ \
  	-p ${row_num} \
  	-e ${params.effect}	
  	cp -r ./*/*/*/*/* .
  	"""

}

process create_dir_structure {

    tag ("Directory name: univariate")	
   
    input:
    tuple val(pheno), val(kinsloco), path(infiles)	

    output:
    path("univariate")
	
    script:
    def outdir = "univariate/${pheno}/null_covars_LOCO/${kinsloco}" + "_" + "${params.effect}".replace(",", "_")
    
    """
    mkdir -p ${outdir}
    mv ${infiles} ${outdir}
	  """

}

process map_DGE_LOCO_GPKronSum {

    tag ("Feature number: ${row_num}")	
    container 'biocorecrg/microbiomepython:0.2'
    label('ignore')
    memory = '35G'
    
    input:
    tuple val(row_num), file(matrix), val(pheno), val(cage), val(dam), val(grm_v), path(micro_fold)	

    output:
    path("univariate")
	
    script:
    """
  	map_DGE_LOCO_GPKronSum.py \
  	${matrix} \
  	${pheno} \
  	None \
  	${cage} \
  	${grm_v} \
  	--out ./ \
  	-p ${row_num} \
  	-e ${params.effect}	
  	"""

}

// 8. Call QTLs and calculate significance for loci association
process call_QTLs {
    memory='10G'
    cpus='1'
    time='48h'
    publishDir(outputGWAS, mode:'copy')
    tag ("Variance decomposition: ${effect}")
    container 'morillofe/microbiomeprofile:0.7'
    
    input:
    tuple val(effect),
    val(kinship),
    val(model),
    val(pheno),
    val(significance_logp),
    val(min_logp),
    path(inputdir),
    path(covariates),
    val(taxonomic_ranks),
    val(functional_ranks),
    path(final_residuals)
	
    output:
	  tuple path("All_QTLs.RData"), path("gwas_qtl_association.pdf.pdf"), path("gwas_statistics.pdf")
	
    script:
    """
  	if [ ${pheno} = "microbiome_taxonomic" ]; then
      	call_QTLs_tax.R -odir ./ \
          	-model ${model} \
          	-pheno ${pheno} \
          	-kinloco ${kinship} \
          	-eff ${effect} \
          	-sig_l ${significance_logp} \
          	-min_l ${min_logp} \
          	-covariates ${covariates} \
            -ranks ${taxonomic_ranks} \
            -res ${final_residuals}
    else
      	call_QTLs_func.R -odir ./ \
          	-model ${model} \
          	-pheno ${pheno} \
          	-kinloco ${kinship} \
          	-eff ${effect} \
          	-sig_l ${significance_logp} \
          	-min_l ${min_logp} \
          	-covariates ${covariates} \
            -ranks ${functional_ranks} \
            -res ${final_residuals}
    fi
    """

}

// 9. Generate Residuals boxplots for each significant trait x SNP association
process genotypes {
    publishDir(outputGWAS, mode:'copy')
    tag ("Mapped genotypes: ${effect}")
    container 'morillofe/microbiomeprofile:0.7'
    
    input:
    tuple path(qtls),
    path(filtered_residuals),
    path(h5),
    path(metadata),
    path(covariates),
    val(sample_identifier),
    val(host_identifier),
    val(method),
    val(taxonomic_ranks)

    output:
	  path("genotypes_boxplot.pdf")
	
    script:
    """
    genotypes.R \
      	-meta ${metadata} \
      	-covariates ${covariates} \
      	-res ${filtered_residuals} \
        -sid ${sample_identifier} \
        -hid ${host_identifier} \
      	-method ${method} \
      	-ranks ${taxonomic_ranks} \
        -qtls ${qtls} \
        -h5 ${h5}
    """

}
