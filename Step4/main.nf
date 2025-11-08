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

/* ========================== PARAMS (defaults + echo) ========================= */
params.help   = params.help   ?: false
params.resume = params.resume ?: false

params.outdir                    = params.outdir                    ?: "./results"
params.h5_path                   = params.h5_path                   ?: null
params.genes_xlsx                = params.genes_xlsx                ?: null
params.metadata_tsv              = params.metadata_tsv              ?: null
params.grm_rdata                 = params.grm_rdata                 ?: null
params.resid_rda                 = params.resid_rda                 ?: null
params.herit_rdata               = params.herit_rdata               ?: null
params.taxonomic_file_taxonomy   = params.taxonomic_file_taxonomy   ?: null
params.taxonomic_file_phylotree  = params.taxonomic_file_phylotree  ?: null
params.taxonomic_parameter_ranks = params.taxonomic_parameter_ranks ?: "Phylum,Class,Order,Family,Genus,Species"

params.gcta_exec = params.gcta_exec ?: "/usr/local/bin/gcta64"
params.plink_exec = params.plink_exec ?: "plink"

params.trait               = params.trait               ?: ["Taxa"] // Taxa, Guild, Beta
params.genera              = params.genera              ?: ["Prevotella"]
params.genus_trait_prefix   = params.genus_trait_prefix   ?: "g__"
params.species_trait_prefix = params.species_trait_prefix ?: "s__"

params.maf_min     = params.maf_min     ?: 0.02
params.miss_max    = params.miss_max    ?: 0.10
params.mac_min     = params.mac_min     ?: 5
params.hwe_p_min   = params.hwe_p_min   ?: 1e-6
params.hwe_maf_min = params.hwe_maf_min ?: 0.05

params.seed                  = params.seed                  ?: 17
params.candidate_window_up   = params.candidate_window_up   ?: 0
params.candidate_window_down = params.candidate_window_down ?: 0

params.set_test  = params.set_test  ?: "SKATO"
params.var_block = params.var_block ?: 500

params.acat_weight  = params.acat_weight  ?: "h2"
params.acat_p_floor = params.acat_p_floor ?: 1e-15
params.q_gene       = params.q_gene       ?: 0.10

params.sentinel_ld_r2 = params.sentinel_ld_r2 ?: 0.10
params.sentinel_ld_kb = params.sentinel_ld_kb ?: 250
params.sentinel_p1    = params.sentinel_p1    ?: 1.0
params.sentinel_p2    = params.sentinel_p2    ?: 1.0
params.min_species    = params.min_species    ?: 6

params.kfolds        = params.kfolds        ?: 5
params.crossfit_seed = params.crossfit_seed ?: params.seed

// Optional inputs for synthesize_results.R; fallback to resid_rda if not provided
params.beta_rda     = params.beta_rda     ?: "/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Tax_Mat_Net/Diversity/residuals_qned_counts_beta.RData"
params.alpha_rda    = params.alpha_rda    ?: "/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Tax_Mat_Net/Diversity/residuals_qned_counts_alpha.RData"
params.clusters_rda = params.clusters_rda ?: "/users/abaud/fmorillo/NextFlow/Microbiome_Profiling/step2/Output_Tax_Mat_Net/Cluster_Analyses/residuals_qned_counts_clusters.RData"

// ===== Mediation (new) =====
params.mediation_enable         = params.mediation_enable         ?: true
params.med_bh_q                 = params.med_bh_q                 ?: 0.10
params.med_gate_from_pgls_q     = params.med_gate_from_pgls_q     ?: params.q_gene
params.med_method               = params.med_method               ?: 'genesis'
params.med_boot                 = params.med_boot                 ?: false
params.med_boot_sims            = params.med_boot_sims            ?: 0

// Outcomes & mediator names (must match rownames in residual matrices)
params.med_beta_trait_name      = params.med_beta_trait_name      ?: 'beta__PD_PC1'
params.med_alpha_trait_name     = params.med_alpha_trait_name     ?: 'alpha__PD_q2'
params.med_guild_trait_name     = params.med_guild_trait_name     ?: 'cluster__3'
params.med_glucose_rdata        = params.med_glucose_rdata        ?: '/users/abaud/fmorillo/paper_figures/pehotypes_residuals.RData'
params.med_glucose_column       = params.med_glucose_column       ?: 'Glucose'
params.med_bmi_column           = params.med_bmi_column           ?: 'BMI'
params.outcomes                 = params.outcomes                 ?: 'all'

params.container_bioconductor = params.container_bioconductor ?: "docker://morillofe/microbiome_geno_bioconductor:step4_v4"
params.container_plink        = params.container_plink        ?: "docker://quay.io/biocontainers/plink:1.90b6.21--h779adbc_1"
params.container_gcta         = params.container_gcta         ?: "docker://morillofe/microbiome_geno_gcta:step4_v1"
params.container_python       = params.container_python       ?: "docker://morillofe/python-ps:3.11-slim"
params.container_genesis      = params.container_genesis      ?: "docker://morillofe/genesis-nf:bioc3.18-ext.v2"

// Genetic correlations / BLUP
params.rg_enable         = params.rg_enable         ?: true
params.rg_partition      = params.rg_partition      ?: true
params.he_enable         = params.he_enable         ?: true
params.cvblup_enable     = params.cvblup_enable     ?: false
params.blup_enable       = params.blup_enable       ?: true
params.cvblup_align_enable = params.cvblup_align_enable ?: true
params.cvblup_outdir     = params.cvblup_outdir     ?: "${params.outdir}/phaseX_blup"

// traits we’ll build phenos for (columns in pheno file)
params.rg_traits = params.rg_traits ?: ["beta__PD_PC1","alpha__PD_q2","cluster__3",'Glucose','BMI']

// pairs to run
params.rg_pairs      = params.rg_pairs      ?: "genus~beta__PD_PC1,genus~alpha__PD_q2,genus~cluster__3,genus~other_genus,genus~Glucose,genus~BMI"
params.gcta_threads  = params.gcta_threads  ?: 4
params.rg_outdir     = params.rg_outdir     ?: "${params.outdir}/rg"
params.cvblup_outdir = params.cvblup_outdir ?: "${params.outdir}/cvblup" // (single definition)

/* ------------------------------ CONFIG ECHO ------------------------------- */
log.info """
================= CONFIG ==========================
h5_path                     : ${params.h5_path}
genes_xlsx                  : ${params.genes_xlsx}
metadata_tsv                : ${params.metadata_tsv}
grm_rdata                   : ${params.grm_rdata}
resid_rda                   : ${params.resid_rda}
herit_rdata                 : ${params.herit_rdata}
outdir                      : ${params.outdir}
trait                       : ${params.trait}
genera                      : ${params.genera}
genus_trait_prefix          : ${params.genus_trait_prefix}
species_trait_prefix        : ${params.species_trait_prefix}
maf_min/miss_max/mac_min    : ${params.maf_min} / ${params.miss_max} / ${params.mac_min}
hwe_p_min / hwe_maf_min     : ${params.hwe_p_min} / ${params.hwe_maf_min}
candidate_window_up/down    : ${params.candidate_window_up} / ${params.candidate_window_down}
set_test / var_block        : ${params.set_test} / ${params.var_block}
acat_weight / p_floor       : ${params.acat_weight} / ${params.acat_p_floor}
q_gene (BH within genus)    : ${params.q_gene}
min_species                 : ${params.min_species}
LD r2 / kb / p1 / p2        : ${params.sentinel_ld_r2} / ${params.sentinel_ld_kb} / ${params.sentinel_p1} / ${params.sentinel_p2}
bioconductor/plink/gcta/genesis : ${params.container_bioconductor} / ${params.container_plink} / ${params.container_gcta} / ${params.container_genesis}
====================================================
"""

/* ============================ PROCESSES ===================================== */

/* ---- Phase 0 ---- */
process H5_TO_GDS_CAND {
  container params.container_bioconductor
  publishDir "${params.outdir}/phase0", mode: 'copy'
  input:
    path h5_path
    path genes_xlsx
  output:
    path "genotypes_candidates.gds", emit: gds_candidates
    path "genotypes_candidates.seq.gds", emit: gds_candidates_seq
    path "genotypes_candidates_scan_lookup.tsv", emit: scan_lookup
    path "candidate_mapping_summary.tsv"
    path "candidate_mapping_total.tsv"
    path "candidate_mapping_per_gene.tsv"
  script:
  """
  echo "[phase0] Converting HDF5 -> GDS (candidate windows)"
  cp ${h5_path} tmp_cand.h5
  h5_to_gds_candidates.R \
    --h5          tmp_cand.h5 \
    --genes       "${genes_xlsx}" \
    --vers        "candidate" \
    --out         genotypes_candidates.gds \
    --window_up   ${params.candidate_window_up} \
    --window_down ${params.candidate_window_down}
  rm -f tmp_cand.h5
  """
}

process BUILD_SRM {
  container params.container_bioconductor
  publishDir "${params.outdir}/phase0", mode: 'copy'
  input:
    path metadata_tsv
  output:
    path "srm_dam.rds", emit: dam_rds
    path "srm_cage.rds", emit: cage_rds
    path "srm_ids.txt"
    tuple path("srm_dam.grm.gz"),  path("srm_dam.grm.id"),  emit: dam_tuple
    tuple path("srm_cage.grm.gz"), path("srm_cage.grm.id"),  emit: cage_tuple
  script:
  """
  echo "[phase0] Building SRM matrices (dam/cage)"
  build_srm_mats.R \
    --metadata "${metadata_tsv}" \
    --outdir   .
  """
}

process MAKE_GENESIS_COHORT_IDS {
  container params.container_bioconductor
  publishDir "${params.outdir}/cohort_ids", mode: 'copy'
  input:
    path gds_seq
    path scan_lookup
    path resid_rda
    path metadata_tsv
    path grm_rdata
    path dam_rds
    path cage_rds
  output:
    path "genesis_final_scan_ids.txt",  emit: scan_ids
    path "genesis_final_label_ids.txt", emit: label_ids
  script:
  """
  make_genesis_cohort_ids.R \
    --gds "${gds_seq}" --scan_lookup "${scan_lookup}" \
    --trait "${params.trait}" \
    --resid_rda "${resid_rda}" --metadata "${metadata_tsv}" \
    --grm_rdata "${grm_rdata}" --srm_dam "${dam_rds}" --srm_cage "${cage_rds}" \
    --outdir .
  """
}

/* ---- Phase 1 ---- */
process GDS_TO_PLINK {
  container params.container_bioconductor
  publishDir "${params.outdir}/phase1", mode: 'copy'
  input:
    path gds
    path cohort_scan_ids
  output:
    tuple path("candidates.bed"), path("candidates.bim"), path("candidates.fam"), emit: plink
  script:
  """
  echo "[phase1] Exporting cohort-filtered GDS -> PLINK (require 0/1/2 present)"
  gds_to_plink.R \
    --gds "${gds}" \
    --out candidates \
    --scan_ids "${cohort_scan_ids}" \
    --min_dosage_levels 3 \
    --block_size 5000
  """
}

process MAKE_GENEFILES {
  container params.container_bioconductor
  publishDir "${params.outdir}/phase1", mode: 'copy'
  input:
    tuple path(bed), path(bim), path(fam)
    path genes_xlsx
  output:
    path "magma.genes.loc", emit: gene_loc
  script:
  """
  echo "[phase1] Gene annotation (w/ windows)"
  make_gene_annotation.R \
    --bim        "${bim}" \
    --genes      "${genes_xlsx}" \
    --magma_up   ${params.candidate_window_up} \
    --magma_down ${params.candidate_window_down} \
    --outdir     .
  """
}

/* ---- Phase 2 ---- */
process PLINK_PREP_FOR_LOCO {
  container params.container_plink
  publishDir "${params.outdir}/phase1", mode: 'copy'
  input:
    tuple path(bed), path(bim), path(fam)
    path scan_ids
  output:
    tuple path("candidates_kept.bed"), path("candidates_kept.bim"), path("candidates_kept.fam"), emit: kept_plink
  script:
  """
  set -euo pipefail
  echo "[phase2] Subset PLINK to cohort"
  awk 'NR==FNR{ids[\$1]=1; next} (\$2 in ids){print \$1, \$2}' "${scan_ids}" candidates.fam > keep_2col.txt
  nmatch=\$(wc -l < keep_2col.txt || echo 0); echo "[phase2] matched to FAM: \${nmatch}"
  [ "\${nmatch}" -gt 0 ] || { echo "[phase2] ERROR: 0 matches"; exit 3; }
  plink --bfile candidates --keep keep_2col.txt --allow-no-sex --make-bed --out candidates_kept
  """
}

/* ---- Phase 3: Single-variant only (kept, optional) ---- */
process GENESIS_SPECIES_SINGLE {
  label "genesis"
  container params.container_bioconductor
  publishDir "${params.outdir}/phase3_species_single", mode: 'copy'
  input:
    path gds_seq; path scan_lookup; path cohort_scan_ids; path cohort_label_ids
    path resid_rda; path metadata_tsv; path grm_rdata; path dam_rds; path cage_rds
    path loco_rds; path loco_bins; path loco_Ns; path loco_ids
    val  species_prefix
    tuple val(species_one), val(genus_one)
  output:
    path "assoc_species_full__*", emit: assoc_dir_part
  tag { "species_single:${species_one}" }
  script:
  """
  set -euo pipefail
  sp="${species_one.trim()}"
  san="${species_one.trim().replaceAll('[^A-Za-z0-9_]', '_')}"
  outdir="assoc_species_full__\${san}"
  mkdir -p "\$outdir" "\$outdir/by_trait"

  run_genesis_species.R \
    --gds "${gds_seq}" \
    --scan_lookup "${scan_lookup}" \
    --scan_ids "${cohort_scan_ids}" \
    --label_ids "${cohort_label_ids}" \
    --trait "${params.trait}" \
    --resid_rda "${resid_rda}" \
    --metadata "${metadata_tsv}" \
    --grm_rdata "${grm_rdata}" \
    --srm_dam "${dam_rds}" \
    --srm_cage "${cage_rds}" \
    --genera "${genus_one.trim()}" \
    --trait_prefix "${species_prefix}" \
    --traits "\${sp}" \
    --maf_min ${params.maf_min} \
    --mac_min ${params.mac_min} \
    --miss_max ${params.miss_max} \
    --hwe_p_min ${params.hwe_p_min} \
    --hwe_maf_min ${params.hwe_maf_min} \
    --loco_grm_rds "${loco_rds}" \
    --outdir "\${san}"

  src="\${san}/by_trait/\${sp}_assoc.tsv.gz"
  if [[ -f "\$src" ]]; then
    cp -n "\$src" "\$outdir/by_trait/"
  fi
  """
}

process MERGE_ASSOC_SPECIES_SINGLE {
  container params.container_bioconductor
  publishDir "${params.outdir}/phase3_species_single", mode: 'copy'
  input:
    path assoc_parts, arity: '1..*'
  output:
    path "assoc_species_full", emit: merged
  script:
  """
  set -euo pipefail
  outdir="assoc_species_full"
  mkdir -p "\$outdir/by_trait"
  for d in ${assoc_parts}; do
    if [ -d "\$d/by_trait" ]; then
      cp -n "\$d"/by_trait/* "\$outdir/by_trait/" || true
    fi
  done
  ls -l "\$outdir/by_trait" >/dev/null
  """
}

/* ---- Phase 5: Genus sentinels (fold-tagged) ---- */
process GENESIS_GENUS_SINGLE_FOLD {
  label "genesis"
  container params.container_bioconductor
  publishDir "${params.outdir}/phase3_genus_single_folds", mode: 'copy'
  input:
    path gds_seq; path scan_lookup
    path resid_rda; path metadata_tsv; path grm_rdata; path dam_rds; path cage_rds
    path loco_rds; path loco_bins; path loco_Ns; path loco_ids
    val  genus_prefix
    tuple val(genus_one), path(train_ids_file), path(train_label_ids_file), val(fold_tag)
  output:
    path "fold*/assoc_genus_full__*", emit: assoc_dir_fold
  tag { "genus_single_fold:${fold_tag}:${genus_one}" }
  script:
  """
  set -euo pipefail
  fold="${fold_tag}"
  gn="${genus_one.trim()}"
  san=\$(echo "\$gn" | sed 's/[^A-Za-z0-9_]/_/g')
  outdir="fold_\${fold}/assoc_genus_full__\${san}"
  mkdir -p "\$outdir"
  run_genesis_genus.R \
    --gds          "${gds_seq}" \
    --scan_lookup  "${scan_lookup}" \
    --scan_ids     "${train_ids_file}" \
    --label_ids    "${train_label_ids_file}" \
    --trait        "${params.trait}" \
    --resid_rda    "${resid_rda}" \
    --metadata     "${metadata_tsv}" \
    --grm_rdata    "${grm_rdata}" \
    --srm_dam      "${dam_rds}" \
    --srm_cage     "${cage_rds}" \
    --genera       "\$gn" \
    --trait_prefix "${genus_prefix}" \
    --maf_min      ${params.maf_min} \
    --mac_min      ${params.mac_min} \
    --miss_max     ${params.miss_max} \
    --hwe_p_min    ${params.hwe_p_min} \
    --hwe_maf_min  ${params.hwe_maf_min} \
    --loco_grm_rds "${loco_rds}" \
    --outdir       "\$outdir"
  """
}

process SELECT_SENTINELS_GENUS_FOLD {
  container params.container_plink
  publishDir "${params.outdir}/phase5_sentinels_folds", mode: 'copy'
  input:
    tuple path(assoc_genus_fold_dir), val(fold_tag), val(genus_san)
    tuple path(bed), path(bim), path(fam)
    path gene_loc
    val  genus_prefix
  output:
    path "sentinels_genus__*__fold*.tsv", emit: sentinels_fold
    tuple val(genus_san), val(fold_tag), path("snps_for_med__*__fold*.tsv"), emit: snps_for_med
  tag { "sentinels_fold:${fold_tag}:${genus_san}" }
  script:
  """
  select_sentinels_genus.sh \
    "${assoc_genus_fold_dir}" \
    "${bim}" \
    "${gene_loc}" \
    "${genus_prefix}" \
    "${params.candidate_window_up}" \
    "${params.candidate_window_down}"
  gensan=\$(echo "${genus_san}" | sed 's/[^A-Za-z0-9_]/_/g')
  mv sentinels_genus.tsv "sentinels_genus__\${gensan}__fold${fold_tag}.tsv"
  INPUT_FILE="sentinels_genus__\${gensan}__fold${fold_tag}.tsv"
  OUTPUT_FILE="snps_for_med__\${gensan}__fold${fold_tag}.tsv"
  awk -F'\t' 'BEGIN {print "snp"} NR > 1 {print \$6 ":" \$7}' "\$INPUT_FILE" > "\$OUTPUT_FILE"
  """
}

/* ---- Phase 6: OOF species betas per fold ---- */
process MAKE_LABELS_FOR_IDS_TEST {
  container params.container_bioconductor
  publishDir "${params.outdir}/kfold_ids", mode:'copy'
  input:
    tuple val(fold_tag), path(scan_ids_subset)
    path cohort_scan_all
    path cohort_label_all
  output:
    tuple val(fold_tag), path("label_ids_test.txt"), emit: labels_test
  script:
  """
  set -euo pipefail
  paste "${cohort_scan_all}" "${cohort_label_all}" > cohort_scan_label.tsv
  awk 'NR==FNR{map[\$1]=\$2; next} { if(!(\$1 in map)){print "[fold ${fold_tag}] missing " \$1 > "/dev/stderr"; exit 2} print map[\$1] }' \
      cohort_scan_label.tsv "${scan_ids_subset}" > label_ids_test.txt
  nscan=\$(wc -l < "${scan_ids_subset}"   || echo 0)
  nlab=\$(wc -l < label_ids_test.txt      || echo 0)
  [ "\$nlab" -eq "\$nscan" ] || { echo "[fold ${fold_tag}] ERROR: label_ids length (\$nlab) != scan_ids length (\$nscan)"; exit 1; }
  """
}

process SPECIES_BETAS_FOR_SENTINELS_GENUS_FOLD {
  label "genesis"
  container params.container_bioconductor
  publishDir "${params.outdir}/phase6_beta_cf_folds", mode:'copy'
  input:
    path gds_seq; path scan_lookup
    path resid_rda; path metadata_tsv; path grm_rdata; path dam_rds; path cage_rds
    path loco_rds; path loco_bins; path loco_Ns; path loco_ids
    val  species_prefix
    tuple val(genus_san), val(fold_tag), path(sentinels_fold_tsv), path(test_ids), path(test_labels)
  output:
    path "beta_cf_by_species__from_genus__*__fold*.tsv", emit: beta_cf_fold
  tag { "species_oof:${fold_tag}:${genus_san}" }
  script:
  """
  set -euo pipefail
  gensan="${genus_san}"
  if [ "${params.trait}" = "Taxa" ]; then
    genus_name="\${gensan//_/ }"   # only Taxa wants spaces
  else
    genus_name="\$gensan"          # Beta/Guild keep underscores
  fi
  out="beta_cf_by_species__from_genus__\${gensan}__fold${fold_tag}.tsv"
  run_species_betas_for_sentinels_fold.R \
    --gds "${gds_seq}" \
    --scan_lookup "${scan_lookup}" \
    --scan_ids_test "${test_ids}" \
    --label_ids_test "${test_labels}" \
    --trait "${params.trait}" \
    --resid_rda "${resid_rda}" \
    --metadata "${metadata_tsv}" \
    --grm_rdata "${grm_rdata}" \
    --srm_dam "${dam_rds}" \
    --srm_cage "${cage_rds}" \
    --loco_grm_rds "${loco_rds}" \
    --genus "\${genus_name}" \
    --species_prefix "${species_prefix}" \
    --sentinels_tsv "${sentinels_fold_tsv}" \
    --out "\${out}"
  """
}

/* ---- Phase 6b: META (build consensus + meta) ---- */
process META_BETA_CF_FOLDS {
  container params.container_bioconductor
  publishDir "${params.outdir}/phase6_beta_cf", mode: 'copy'
  input:
    tuple val(genus_san), path(beta_cf_folds), path(sentinel_folds)
  output:
    path "beta_cf_by_species__from_genus_crossfit__*.tsv", emit: beta_cf_crossfit
    path "consensus_sentinels__*.tsv",                     emit: consensus_sentinels
  script:
  """
  set -euo pipefail
  gensan=\$(echo "$genus_san" | sed 's/[^A-Za-z0-9_]/_/g')
  out="beta_cf_by_species__from_genus_crossfit__\${gensan}.tsv"
  cons="consensus_sentinels__\${gensan}.tsv"

  printf "%s\n" ${beta_cf_folds}  > inputs.txt
  printf "%s\n" ${sentinel_folds} > sentinels.txt

  meta_beta_cf_across_folds.R \
    --inputs-list     inputs.txt \
    --sentinels-list  sentinels.txt \
    --out             "\${out}" \
    --consensus       "\${cons}"

  echo "[meta] WROTE: \$out and \$cons"
  """
}

/* ---- Phase 7: PGLS per genus on cross-fit β ---- */
process CORRELATE_H2_PGLS_GENUS_TAGGED {
  label 'big_mem'
  // You can swap container to params.container_genesis if desired
  container "/users/abaud/fmorillo/singularity_containers/genesis-nf_bioc3.18-ext.v3.sif"
  publishDir "${params.outdir}/phase7_pgls", mode: 'copy'
  input:
    tuple val(genus_san), path(beta_cf_by_species_genus)
    path herit_rdata
    path taxonomy_file
    path phylotree_file
  output:
    path "correlations_by_gene__*__EXCLUDED.tsv", emit: corr_excluded
    path "correlations_by_gene__*.tsv",           emit: corr_tsvs_genus
    path "correlation_plots/*"
  script:
  """
  set -euo pipefail
  out="correlations_by_gene__${genus_san}.tsv"
  correlate_h2_simple.R \
    --beta_cf ${beta_cf_by_species_genus} \
    --herit_rdata ${herit_rdata} \
    --taxonomy ${taxonomy_file} \
    --phylotree ${phylotree_file} \
    --tax_ranks "${params.taxonomic_parameter_ranks}" \
    --min_species ${params.min_species} \
    --branch genus_sentinels \
    --q ${params.q_gene} \
    --transform rankz \
    --report ML \
    --branchlen none \
    --outdir .
  mv correlations_by_gene__genus_sentinels.tsv "\${out}"
  """
}

process PSEUDO_CORRELATE {
  publishDir "${params.outdir}/phase7_pgls", mode: 'copy'
  input:
    tuple val(genus_san), path(beta_cf_by_species_genus)
  output:
    path "correlations_by_gene__*.tsv", emit: corr_tsvs_genus
  script:
  """
  out="correlations_by_gene__${genus_san}.tsv"
  cp ${beta_cf_by_species_genus} "\${out}"
  """
}


/* ---- Merge fold sentinels per genus ---- */
process MERGE_SENTINELS_GENUS_FOLDS {
  container params.container_bioconductor
  publishDir "${params.outdir}/phase5_sentinels_folds", mode:'copy'
  input:
    tuple val(genus_san), path(sentinels_fold_files)
  output:
    path "sentinels_genus__*__crossfit.tsv", emit: sentinels_cf_one
  script:
  """
  set -euo pipefail
  gensan=\$(echo "${genus_san}" | sed 's/[^A-Za-z0-9_]/_/g')
  out="sentinels_genus__\${gensan}__crossfit.tsv"
  files=( ${sentinels_fold_files} )
  if [ "\${#files[@]}" -lt 1 ]; then
    echo "[merge] ERROR: no per-fold inputs for genus ${genus_san}" >&2
    exit 2
  fi
  Rscript - "\${out}" "\${files[@]}" <<'RS'
    suppressPackageStartupMessages(library(data.table))
    args  <- commandArgs(trailingOnly = TRUE)
    out   <- args[1]
    files <- args[-1]
    if (length(files) == 0L) stop("no input files")
    DT <- rbindlist(lapply(files, fread), fill = TRUE)
    setDT(DT); setnames(DT, tolower(names(DT)))
    if (!"genus" %in% names(DT)) DT[, genus := NA_character_]
    if (!"gene"  %in% names(DT) && "gene_id" %in% names(DT)) setnames(DT, "gene_id", "gene")
    if (!"sentinel_snp" %in% names(DT) && "snp" %in% names(DT)) setnames(DT, "snp", "sentinel_snp")
    OUT <- unique(DT[, .(genus, gene, sentinel_snp)])
    fwrite(OUT, file = out, sep = "\t")
RS
  """
}

/* ---- Phase 8: Synthesis (GDS-based) ---- */
process SYNTHESIZE_RESULTS_SIMPLE {
  container params.container_bioconductor
  publishDir { "${params.outdir}/phase8_synthesis_simple/${genus_san}" }, mode: 'copy'
  input:
    path gds_seq
    path scan_lookup
    path cohort_scan_ids
    path cohort_label_ids
    tuple val(genus_name), val(genus_san), path(corr_by_gene_genus)
    path genes_xlsx
  output:
    tuple val(genus_san), path("selected_genes.RData"), emit: selected_rdata
  tag { "synthesize_simple:${genus_name}" }
  script:
  """
  set -euo pipefail
  synthesize_results.R \
    --gds "${gds_seq}" \
    --scan_lookup "${scan_lookup}" \
    --scan_ids "${cohort_scan_ids}" \
    --label_ids "${cohort_label_ids}" \
    --corr "${corr_by_gene_genus}" \
    --genes_xlsx "${genes_xlsx}"
  """
}

/* ---- Cross-fit helpers (RUN ONCE) ---- */
process MAKE_KFOLDS {
  container params.container_bioconductor
  publishDir "${params.outdir}/kfold_ids", mode: 'copy'
  input:
    path cohort_scan_ids
  output:
    path "fold*/scan_ids_train.txt", emit: train_ids
    path "fold*/scan_ids_test.txt",  emit: test_ids
  script:
  """
  make_kfolds.R \
    --scan_ids ${cohort_scan_ids} \
    --k ${params.kfolds} \
    --seed ${params.crossfit_seed} \
    --outdir .
  """
}

process MAKE_LABELS_FOR_IDS_TRAIN {
  container params.container_bioconductor
  publishDir "${params.outdir}/kfold_ids", mode:'copy'
  input:
    tuple val(fold_tag), path(scan_ids_subset)
    path cohort_scan_all
    path cohort_label_all
  output:
    tuple val(fold_tag), path("label_ids_train.txt"), emit: labels_train
  script:
  """
  set -euo pipefail
  paste "${cohort_scan_all}" "${cohort_label_all}" > cohort_scan_label.tsv
  awk 'NR==FNR{map[\$1]=\$2; next} { if(!(\$1 in map)){print "[fold \${fold_tag}] missing " \$1 > "/dev/stderr"; exit 2} print map[\$1] }' \
      cohort_scan_label.tsv "${scan_ids_subset}" > label_ids_train.txt
  nscan=\$(wc -l < "${scan_ids_subset}"   || echo 0)
  nlab=\$(wc -l < label_ids_train.txt     || echo 0)
  [ "\$nlab" -eq "\$nscan" ] || { echo "[fold \${fold_tag}] ERROR: label_ids length (\$nlab) != scan_ids length (\$nscan)"; exit 1; }
  """
}

process MAKE_LOCO_GRMS {
  label 'big_mem'
  container params.container_gcta
  publishDir "${params.outdir}/loco_grm", mode:'copy'
  input:
    tuple path(bed), path(bim), path(fam)
    path scan_ids
  output:
    path "grm_loco_index.tsv", emit:loco_index
    path "loco_chr*.grm.bin",  emit:loco_bins
    path "loco_chr*.grm.N.bin",emit:loco_Ns
    path "loco_chr*.grm.id",   emit:loco_ids
  script:
    "make_loco_grms.sh --bfile candidates_kept --outdir . --gcta_exec ${params.gcta_exec}"
}

process PACK_LOCO_RDS {
  container params.container_bioconductor
  publishDir "${params.outdir}/loco_grm", mode:'copy'
  input:
    path loco_index
  output:
    path "grm_loco.rds", emit:loco_rds
  script:
    "Rscript -e 'idx <- read.table(\"grm_loco_index.tsv\", header=TRUE, sep=\"\t\", stringsAsFactors=FALSE); vec <- setNames(idx\$prefix, idx\$chr); saveRDS(vec, file=\"grm_loco.rds\")'"
}

/* ---- Phase 9: Mediation gates & estimation ---- */

// 9.0 Gate by PGLS (per genus)
process PGLS_KEEP_VARIANTS {
  container params.container_bioconductor
  publishDir "${params.outdir}/phase9_mediation", mode: 'copy'
  input:
    tuple val(genus_name), val(genus_san), path(corr_tsv)
  output:
    tuple val(genus_san), path("pgls_keep__*.tsv"), emit: pgls_keep
  tag { "pgls-keep:${genus_san}" }
  script:
  """
  set -euo pipefail
  out="pgls_keep__${genus_san}.tsv"
  pgls_keep_variants.R \
    --corr "${corr_tsv}" \
    --q ${params.med_gate_from_pgls_q} \
    --out "\${out}"
  echo "[pgls-keep] WROTE \${out}"
  """
}

// 9.1 BH per fold on sentinels (Score.pval)
process FILTER_SENTINELS_BH_PER_FOLD {
  container "/users/abaud/fmorillo/singularity_containers/genesis-nf_bioc3.18-ext.v3.sif"
  publishDir "${params.outdir}/phase9_mediation", mode:'copy'
  input:
    tuple val(genus_san), val(fold_tag), path(sentinels_fold_tsv), path(pgls_keep_tsv)
  output:
    tuple val(genus_san), val(fold_tag), path("bh_keep__*__fold*.tsv"), emit: bh_kept
  tag { "bh:${fold_tag}:${genus_san}" }
  script:
  """
  set -euo pipefail
  gensan="${genus_san}"
  out="bh_keep__\${gensan}__fold${fold_tag}.tsv"
  Rscript - "${sentinels_fold_tsv}" "${pgls_keep_tsv}" "\${out}" "${params.med_bh_q}" <<'RS'
  suppressPackageStartupMessages(library(data.table))
  args <- commandArgs(trailingOnly=TRUE)
  sent_file <- args[1]; keep_file <- args[2]; outp <- args[3]; qth <- as.numeric(args[4])

  S <- fread(sent_file, sep="\t")
  K <- fread(keep_file,  sep="\t")
  S\$snp = paste0(S\$CHR, ":", S\$BP)
  print(paste0("Sentinel SNPs before PGLS filter: ",nrow(S)))

  # Filter by PGLS keep, THEN BH (independent filtering)
  S <- S[S\$snp %in% K\$snp,]
  print(paste0("Sentinel SNPs after PGLS filter: ",nrow(S)))
  
  S\$q = p.adjust(S\$P, method="BH")
  fwrite(S[,c("snp","q")], file=outp, sep="\t")
  RS
  echo "[bh] WROTE \${out}"
  """
}

// 9.2 Intersect BH-kept with PGLS-keep per genus, per fold
process INTERSECT_GATES_PER_FOLD {
  container params.container_bioconductor
  publishDir "${params.outdir}/phase9_mediation", mode:'copy'
  input:
    tuple val(genus_san), val(fold_tag), path(bh_keep_tsv), path(pgls_keep_tsv)
  output:
    tuple val(genus_san), val(fold_tag), path("snps_for_med__*__fold*.tsv"), emit: snps_for_med
  tag { "gate:${fold_tag}:${genus_san}" }
  script:
  """
  set -euo pipefail
  out="snps_for_med__${genus_san}__fold${fold_tag}.tsv"
  Rscript -e '
    suppressPackageStartupMessages(library(data.table))
    args <- commandArgs(trailingOnly=TRUE)
    bh <- fread(args[1]); pg <- fread(args[2])
    setnames(bh, tolower(names(bh))); setnames(pg, tolower(names(pg)))
    if (!"snp" %in% names(bh)) stop("bh file lacks snp")
    if (!"snp" %in% names(pg)) stop("pgls file lacks snp")
    M <- merge(unique(bh[, .(snp)]), unique(pg[, .(snp)]), by="snp")
    fwrite(M, file=args[3], sep="\t")
  ' "${bh_keep_tsv}" "${pgls_keep_tsv}" "\${out}"
  echo "[gate] WROTE \${out}"
  """
}

// 9.3 Run mediation per fold (GENESIS + LOCO)
process RUN_MEDIATION_FOLD {
  tag { "${genus_san}:${fold_tag}" }
  cpus 2
  memory '8 GB'
  errorStrategy 'finish'
  publishDir "${params.outdir}/phase9_mediation", mode: 'copy'
  container "/users/abaud/fmorillo/singularity_containers/microbiome_geno_bioconductor_step4_v7.sif"

  input:
  tuple val(genus_san), val(fold_tag), path(snps_tsv), path(selected_genes_rdata), val(other_genus)
  // shared/value-ish:
  path scan_lookup
  path resid_rda
  path beta_rda
  path alpha_rda
  path clusters_rda
  path phen_rda
  val  phen_glucose_col
  val  phen_bmi_col
  path dam_rds
  path cage_rds
  path loco_rds
  path loco_bins
  path loco_Ns
  path loco_ids
  val  trait_prefix
  val  beta_trait
  val  alpha_trait
  val  guild_trait
  val outcomes
  val  mr_enable
  val  mr_fstat_min

  output:
  tuple val(genus_san), val(fold_tag), path("med__${genus_san}__${fold_tag}.tsv"), emit: med_results_fold

  script:
  def out   = "med__${genus_san}__${fold_tag}.tsv"
  """
  gensan="${genus_san}"
  if [ "${params.trait}" = "Taxa" ]; then
    genus_name="\${gensan//_/ }"   # only Taxa wants spaces
  else
    genus_name="\$gensan"          # Beta/Guild keep underscores
  fi
  set -euo pipefail
  run_mediation_fold_genesis.R \
    --selected_genes_rdata "${selected_genes_rdata}" \
    --snps_tsv             "${snps_tsv}" \
    --scan_lookup          "${scan_lookup}" \
    --trait                "${params.trait}" \
    --resid_rda            "${resid_rda}" \
    --beta_rda             "${beta_rda}" \
    --alpha_rda            "${alpha_rda}" \
    --clusters_rda         "${clusters_rda}" \
    --phen_rda             "${phen_rda}" \
    --phen_glucose_col     "${phen_glucose_col}" \
    --phen_bmi_col         "${phen_bmi_col}" \
    --dam_rds              "${dam_rds}" \
    --cage_rds             "${cage_rds}" \
    --loco_rds             "${loco_rds}" \
    --genus                "\${genus_name}" \
    --other_genus          "${other_genus}" \
    --trait_prefix         "${trait_prefix}" \
    --beta_trait           "${beta_trait}" \
    --alpha_trait           "${alpha_trait}" \
    --guild_trait           "${guild_trait}" \
    --outcomes "${params.outcomes}" \
    --mr_enable        "${mr_enable}" \
    --mr_fstat_min     "${mr_fstat_min}" \
    --sig_metric q --fdr_level 0.10 \
    --out                  "${out}"
  echo "[med] wrote ${out}"
  """
}

// 9.4 Meta across folds
process META_MEDIATION_ACROSS_FOLDS {
  container "/users/abaud/fmorillo/singularity_containers/microbiome_geno_bioconductor_step4_v7.sif"
  publishDir "${params.outdir}/phase9_mediation", mode:'copy'
  input:
    tuple val(genus_san), path(med_results_for_genus)
  output:
    path "mediation_meta__*.tsv", emit: med_meta
  tag { "med:meta:${genus_san}" }
  script:
  """
  set -euo pipefail
  printf "%s\n" ${med_results_for_genus} > inputs.txt
  out="mediation_meta__${genus_san}.tsv"
  meta_mediation_across_folds.R \
  --inputs-list inputs.txt \
  --out "\${out}" \
  --sig_metric q \
  --fdr_level ${params.mr_fdr_level}
  echo "[med:meta] WROTE \${out}"
  """
}

process SIG_IV {
  container params.container_bioconductor
  publishDir "${params.outdir}/phase9_mediation", mode:'copy'
  input:
    path(med_meta)
  output:
    path "mediation_sig_iv.tsv"
  tag { "med:sig_iv" }
  script:
  """
  set -euo pipefail
  sig_iv.R \
  --meta ${med_meta} \
  --out  mediation_sig_iv.tsv \
  --q 0.10 --Fmin 10 --steiger_min_prop 0.60 --sign_coh_min 0.60 \
  --delta_frac_tol 0.05
  """
}

process SIG_MVMR {
  container "/users/abaud/fmorillo/singularity_containers/microbiome_geno_bioconductor_step4_v7.sif"
  publishDir "${params.outdir}/phase9_mediation", mode:'copy'
  input:
    path(med_meta)
  output:
    path "mvmr_sig.tsv"
  script:
  """
  Rscript - <<'RS' "${med_meta}" ${params.mr_fdr_level}
    suppressPackageStartupMessages(library(data.table))
    args <- commandArgs(trailingOnly=TRUE)
    meta <- args[1]; qthr <- as.numeric(args[2])
    DT <- fread(meta); setDT(DT); setnames(DT, tolower(names(DT)))
    X <- DT[analysis == "mvmr"]
    if (!nrow(X)) { fwrite(X, file="mvmr_sig.tsv", sep="\t"); quit(save="no") }
    X[, q1 := p.adjust(p1_meta, method="BH"), by=.(outcome)]
    X[, q2 := p.adjust(p2_meta, method="BH"), by=.(outcome)]
    keep <- X[(q1 < qthr) | (q2 < qthr)]
    fwrite(keep, file="mvmr_sig.tsv", sep="\t")
  RS
  """
}

/* ============================ RG / CVBLUP HELPERS ============================ */

/* Make all-SNP GRM from candidates_kept.* (single definition; reused by RG & BLUP) */
process MAKE_GRM_ALL {
  container params.container_gcta
  publishDir "${params.rg_outdir}/grm_all", mode:'copy'
  input:
    tuple path(bed), path(bim), path(fam)
  output:
    tuple path("grm_all.grm.bin"), path("grm_all.grm.N.bin"), path("grm_all.grm.id"), emit: grm_all
  script:
  """
  set -euo pipefail
  
  # Find a runnable gcta64 (works for plain binary or AppImage)
  GCTA_PATH="\$(command -v ${params.gcta_exec} || true)"
  [[ -n "\$GCTA_PATH" ]] || { echo "[mgrm] ERROR: cannot find ${params.gcta_exec}"; exit 127; }
  if "\$GCTA_PATH" --help >/dev/null 2>&1; then
    RUN_GCTA() { "\$GCTA_PATH" "\$@"; }
  elif "\$GCTA_PATH" --appimage-help >/dev/null 2>&1; then
    RUN_GCTA() { "\$GCTA_PATH" --appimage-extract-and-run "\$@"; }
  else
    echo "[mgrm] ERROR: gcta64 not runnable"; exit 127
  fi
  
  RUN_GCTA --bfile candidates_kept --autosome --make-grm-bin --thread-num ${params.gcta_threads} --out grm_all
  """
}

process BUILD_MGRM_WITH_SRM {
  container params.container_gcta
  publishDir "${params.rg_outdir}/grm_all", mode:'copy'
  input:
    tuple path(bed), path(bim), path(fam)
    tuple path(dam_grm_gz),  path(dam_id)
    tuple path(cage_grm_gz), path(cage_id)
    path  scan_lookup
    path  cohort_scan_ids
  output:
    path "grm_all_txt.grm.gz", emit: grm_all_txt_gz
    path "grm_all_txt.grm.id", emit: grm_all_txt_id
    path "grm_all_harmo.grm.gz", emit: grm_all_harmo_gz
    path "grm_all_harmo.grm.id", emit: grm_all_harmo_id
    path "srm_dam_harmo.grm.gz", emit: srm_dam_harmo_gz
    path "srm_dam_harmo.grm.id", emit: srm_dam_harmo_id
    path "srm_cage_harmo.grm.gz", emit: srm_cage_harmo_gz
    path "srm_cage_harmo.grm.id", emit: srm_cage_harmo_id
    path "mgrm_with_srm.txt",    emit: mgrm_with_srm
  script:
  """
  build_mgrm_with_srm.sh \
    --gcta_exec ${params.gcta_exec} \
    --threads   ${params.gcta_threads} \
    --bfile     candidates_kept \
    --dam_gz    ${dam_grm_gz} --dam_id  ${dam_id} \
    --cage_gz   ${cage_grm_gz} --cage_id ${cage_id} \
    --scan_lookup ${scan_lookup} \
    --scan_ids    ${cohort_scan_ids} \
    --mgrm_out  mgrm_with_srm.txt
  """
}

/* Build pheno with FID/IID + genus + requested rg_traits; and pairs.tsv */
process BUILD_GCTA_PHENOS {
  container params.container_bioconductor
  publishDir "${params.rg_outdir}/pheno", mode:'copy'
  input:
    tuple path(bed), path(bim), path(fam)  // <- accept PLINK tuple; no [2] indexing needed
    path scan_lookup
    path resid_rda
    path beta_rda
    path alpha_rda
    path clusters_rda
  output:
    path "gcta.pheno",    emit: pheno
    path "trait_map.tsv", emit: trait_map
    path "pairs.tsv",     emit: pairs
  script:
  """
  build_gcta_phenos.R \
    --fam ${fam} \
    --scan_lookup ${scan_lookup} \
    --resid_rda ${resid_rda} \
    --beta_rda ${beta_rda} \
    --alpha_rda ${alpha_rda} \
    --clusters_rda ${clusters_rda} \
    --phen_rda             "${phen_rda}" \
    --phen_glucose_col     "${phen_glucose_col}" \
    --phen_bmi_col         "${phen_bmi_col}" \
    --genus "\$(printf %s "${(params.genera as List)[0]}")" \
    --rg_traits "${params.rg_traits.join(",")}" \
    --rg_pairs "${params.rg_pairs}"
  """
}

/* Create SNP extract lists for the focus genus only */
process MAKE_GENESET_EXTRACTS_FOCUS {
  container params.container_bioconductor
  publishDir "${params.rg_outdir}/geneset", mode:'copy'
  input:
    tuple val(genus_san), path(selected_genes_rdata)
    tuple path(bed), path(bim), path(fam)
  output:
    path "snp_extract_focus.txt", emit: snp_focus
    path "snp_extract_bg.txt",    emit: snp_bg
  script:
  """
  make_geneset_extracts_focus.R \
    --selected_genes_rdata ${selected_genes_rdata} \
    --bim ${bim} \
    --out_prefix snp_extract
  """
}

/* Build focus vs background GRMs (for --mgrm) */
process MAKE_GRMS_GENESET_VS_BG {
  container params.container_gcta
  publishDir "${params.rg_outdir}/grm_sets", mode:'copy'
  input:
    tuple path(bed), path(bim), path(fam)
    path snp_focus
    path snp_bg
  output:
    tuple path("grm_focus.grm.bin"), path("grm_focus.grm.N.bin"), path("grm_focus.grm.id"), emit: grm_focus
    tuple path("grm_bg.grm.bin"),    path("grm_bg.grm.N.bin"),    path("grm_bg.grm.id"),    emit: grm_bg
    tuple path("grm_focus.grm.gz"),  path("grm_focus.grm.id"),    emit: grm_focus_gz
    tuple path("grm_bg.grm.gz"),     path("grm_bg.grm.id"),       emit: grm_bg_gz
  script:
  """
  set -euo pipefail
  
  # Find a runnable gcta64 (works for plain binary or AppImage)
  GCTA_PATH="\$(command -v ${params.gcta_exec} || true)"
  [[ -n "\$GCTA_PATH" ]] || { echo "[mgrm] ERROR: cannot find ${params.gcta_exec}"; exit 127; }
  if "\$GCTA_PATH" --help >/dev/null 2>&1; then
    RUN_GCTA() { "\$GCTA_PATH" "\$@"; }
  elif "\$GCTA_PATH" --appimage-help >/dev/null 2>&1; then
    RUN_GCTA() { "\$GCTA_PATH" --appimage-extract-and-run "\$@"; }
  else
    echo "[mgrm] ERROR: gcta64 not runnable"; exit 127
  fi
  
  RUN_GCTA --bfile candidates_kept --extract snp_extract_focus.txt --autosome --make-grm-bin --thread-num ${params.gcta_threads} --out grm_focus
  RUN_GCTA --bfile candidates_kept --extract snp_extract_bg.txt    --autosome --make-grm-bin --thread-num ${params.gcta_threads} --out grm_bg
  RUN_GCTA --bfile candidates_kept --extract snp_extract_focus.txt --autosome --make-grm-gz  --thread-num ${params.gcta_threads} --out grm_focus
  RUN_GCTA --bfile candidates_kept --extract snp_extract_bg.txt    --autosome --make-grm-gz  --thread-num ${params.gcta_threads} --out grm_bg
  """
}

/* Bivariate GREML on all-SNP GRM */
 process RUN_GCTA_BIVAR_RG {
  container params.container_gcta
  publishDir "${params.rg_outdir}/bivar_all", mode:'copy'
  input:
    path mgrm
    path pheno
    path trait_map
    path pairs
  output:
    path "rg_bivar_all__*.hsq", emit: hsq_all
    path "gcta_bivar_rg_summary.csv"
  script:
  """
  set -euo pipefail
  run_gcta_bivar.sh \
    --gcta_exec ${params.gcta_exec} \
    --grm_prefix ${mgrm} \
    --pheno ${pheno} \
    --pairs ${pairs} \
    --threads ${params.gcta_threads} \
    --out_prefix rg_bivar_all__
  retro_name_bivar.sh --pairs ${pairs} --dir . --out_prefix rg_bivar_all__
  summarize_gcta_bivar.R
  """
}



/* Bivariate GREML with two-component GRM */
process RUN_GCTA_BIVAR_RG_MGRM {
  when:
    params.rg_partition
  container params.container_gcta
  publishDir "${params.rg_outdir}/bivar_mgrm", mode:'copy'
  input:
    tuple path(f_gz), path(f_id)          // from grm_focus_gz
    tuple path(b_gz), path(b_id)          // from grm_bg_gz
    path pheno
    path pairs
    tuple path(dam_gz),  path(dam_id)     // from BUILD_SRM.out.dam_tuple
    tuple path(cage_gz), path(cage_id)    // from BUILD_SRM.out.cage_tuple
  output:
    path "rg_bivar_mgrm__*.hsq", emit: hsq_mgrm
  script:
  """
  set -euo pipefail

  # Drop the .grm.gz suffix using Bash parameter expansion
  prefix_f=\${f_gz%.grm.gz}
  prefix_b=\${b_gz%.grm.gz}
  prefix_dam=\${dam_gz%.grm.gz}
  prefix_cage=\${cage_gz%.grm.gz}

  {
    echo "\$prefix_f"
    echo "\$prefix_b"
    echo "\$prefix_dam"
    echo "\$prefix_cage"
  } > mgrm_sets_with_srm.txt

  run_gcta_bivar_mgrm.sh \
    --gcta_exec ${params.gcta_exec} \
    --mgrm mgrm_sets_with_srm.txt \
    --pheno gcta.pheno \
    --pairs pairs.tsv \
    --threads ${params.gcta_threads} \
    --out_prefix rg_bivar_mgrm__
  """
}


/* Bivariate HE regression */
 process RUN_GCTA_HE_BIVAR {
   when:
     params.he_enable
   container params.container_gcta
   publishDir "${params.rg_outdir}/he_bivar", mode:'copy'
   input:
     path mgrm
     path pheno
     path pairs
   output:
     path "he_bivar__*.HEreg" , emit: he_logs
   script:
   """
   run_gcta_he_bivar.sh \
     --gcta_exec ${params.gcta_exec} \
     --grm_prefix ${mgrm} \
     --pheno gcta.pheno \
     --pairs pairs.tsv \
     --threads ${params.gcta_threads} \
     --out_prefix he_bivar__
   """
 }


/* Optional cvBLUP for each phenotype column */
process RUN_GCTA_CVBLUP {
  when:
    params.cvblup_enable
  container params.container_gcta
  publishDir "${params.cvblup_outdir}", mode:'copy'

  input:
    path mgrm_with_srm       // <-- the mgrm_with_srm.txt you built
    path pheno
    path trait_map

  output:
    path "cvblup__*.indi.blp", emit: cv_ind_blp

  script:
  """
  run_gcta_cvblup.sh \
    --gcta_exec ${params.gcta_exec} \
    --grm_prefix ${mgrm_with_srm} \
    --pheno gcta.pheno \
    --trait_map trait_map.tsv \
    --threads ${params.gcta_threads} \
    --out_prefix cvblup__
  """
}


/* Tidy summaries */
process PARSE_RG_OUTPUTS {
  container params.container_bioconductor
  publishDir "${params.rg_outdir}/summary", mode:'copy'
  input:
    path hsq_all 
    path hsq_mgrm 
    path he_logs 
  output:
    path "rg_bivar_all_summary.tsv"  
    path "rg_bivar_mgrm_summary.tsv" 
    path "he_bivar_summary.tsv"      
  script:
  """
  parse_gcta_rg.R --glob "bivar_all/*.hsq"  --out rg_bivar_all_summary.tsv  || true
  parse_gcta_rg.R --glob "bivar_mgrm/*.hsq" --out rg_bivar_mgrm_summary.tsv || true
  parse_gcta_he_bivar.R --glob "he_bivar/*.HEreg" --out he_bivar_summary.tsv || true
  """
}

/* ===== BLUP + alignment ===== */
process MAKE_GCTA_PHENO_GENUS {
  container params.container_bioconductor
  publishDir "${params.cvblup_outdir}/pheno", mode:'copy'
  input:
    path scan_lookup
    path resid_rda
    val  genus_name
    val  trait_prefix
  output:
    path "pheno__g__*.txt", emit: pheno
  script:
  """
  set -euo pipefail
  make_gcta_pheno_genus.R \
    --scan_lookup "${scan_lookup}" \
    --resid_rda   "${resid_rda}" \
    --genus       "${genus_name}" \
    --trait_prefix "${trait_prefix}" \
    --out "pheno__g__${genus_name}.txt"
  """
}

process RUN_GCTA_BLUP_GENUS {
  container params.container_gcta
  publishDir "${params.cvblup_outdir}/blup", mode:'copy'
  input:
    path mgrm_with_srm
    path pheno_file
  output:
    path "blup__g__*.indi.blp", emit: indi_blp
    path "blup__g__*.hsq"
    path "blup__g__*.log"
  script:
  """
  set -euo pipefail
  GENUS=\$(basename "${pheno_file}" .txt | sed 's/^pheno__g__//')
  OUT="blup__g__\${GENUS}"

  # First GRM id file to harmonize IDs
  FIRST_PREFIX=\$(head -n1 "${mgrm_with_srm}")
  FIRST_ID="\${FIRST_PREFIX}.grm.id"
  [[ -s "\${FIRST_ID}" ]] || { echo "[blup] ERROR: missing \${FIRST_ID}"; exit 2; }

  harmonize_pheno_to_grm.sh --pheno "${pheno_file}" --grm_id "\${FIRST_ID}" --out pheno.harmo

  GCTA_PATH="\$(command -v ${params.gcta_exec} || true)"
  [[ -n "\$GCTA_PATH" ]] || { echo "[blup] ERROR: cannot find ${params.gcta_exec}"; exit 127; }
  if "\$GCTA_PATH" --help >/dev/null 2>&1; then
    RUN_GCTA() { "\$GCTA_PATH" "\$@"; }
  elif "\$GCTA_PATH" --appimage-help >/dev/null 2>&1; then
    RUN_GCTA() { "\$GCTA_PATH" --appimage-extract-and-run "\$@"; }
  else
    echo "[blup] ERROR: gcta64 not runnable"; exit 127
  fi

  RUN_GCTA \
    --mgrm-gz "${mgrm_with_srm}" \
    --pheno pheno.harmo \
    --reml \
    --reml-pred-rand \
    --reml-alg-inv 1 \
    --threads 2 \
    --out "\${OUT}"
  """
}


/* ============================ WORKFLOW ====================================== */

workflow {
  /* Inputs */
  def ch_h5    = Channel.fromPath(params.h5_path)
  def ch_genes = Channel.fromPath(params.genes_xlsx)
  def ch_meta  = Channel.fromPath(params.metadata_tsv)
  def ch_grm   = Channel.fromPath(params.grm_rdata)
  def ch_resid_file = switch(params.trait) {
    case 'Taxa'  -> Channel.fromPath(params.resid_rda)
    case 'Guild' -> Channel.fromPath(params.clusters_rda)
    case 'Beta'  -> Channel.fromPath(params.beta_rda)
    case 'Alpha' -> Channel.fromPath(params.alpha_rda)
    default      -> { error "Unknown --trait '${params.trait}'" }
  }
  def ch_resid = ch_resid_file.first()  // value channel
  def ch_herit = Channel.fromPath(params.herit_rdata)
  def ch_tax   = Channel.fromPath(params.taxonomic_file_taxonomy)
  def ch_tree  = Channel.fromPath(params.taxonomic_file_phylotree)

  /* Phase 0/1/2 */
  H5_TO_GDS_CAND(ch_h5, ch_genes)
  def gds_seq_ch = H5_TO_GDS_CAND.out.gds_candidates_seq
  def gds_ch     = H5_TO_GDS_CAND.out.gds_candidates
  def scan_lu_ch = H5_TO_GDS_CAND.out.scan_lookup

  BUILD_SRM(ch_meta)
  def dam_rds_ch  = BUILD_SRM.out.dam_rds
  def cage_rds_ch = BUILD_SRM.out.cage_rds

  MAKE_GENESIS_COHORT_IDS(gds_seq_ch, scan_lu_ch, ch_resid, ch_meta, ch_grm, dam_rds_ch, cage_rds_ch)
  def cohort_scan_ch  = MAKE_GENESIS_COHORT_IDS.out.scan_ids
  def cohort_label_ch = MAKE_GENESIS_COHORT_IDS.out.label_ids

  GDS_TO_PLINK(gds_ch, cohort_scan_ch)
  def plink_ch   = GDS_TO_PLINK.out.plink

  MAKE_GENEFILES(plink_ch, ch_genes)
  def gene_loc_ch = MAKE_GENEFILES.out.gene_loc

  PLINK_PREP_FOR_LOCO(plink_ch, cohort_scan_ch)
  def kept_plink_ch = PLINK_PREP_FOR_LOCO.out.kept_plink
  
  BUILD_MGRM_WITH_SRM(
    kept_plink_ch,
    BUILD_SRM.out.dam_tuple,
    BUILD_SRM.out.cage_tuple,
    H5_TO_GDS_CAND.out.scan_lookup.first(),
    MAKE_GENESIS_COHORT_IDS.out.scan_ids.first()
  )
  def mgrm_with_srm_ch = BUILD_MGRM_WITH_SRM.out.mgrm_with_srm

  MAKE_LOCO_GRMS(kept_plink_ch, cohort_scan_ch)
  def loco_index_ch = MAKE_LOCO_GRMS.out.loco_index
  def loco_bins_ch  = MAKE_LOCO_GRMS.out.loco_bins
  def loco_Ns_ch    = MAKE_LOCO_GRMS.out.loco_Ns
  def loco_ids_ch   = MAKE_LOCO_GRMS.out.loco_ids

  PACK_LOCO_RDS(loco_index_ch)
  def loco_rds_ch = PACK_LOCO_RDS.out.loco_rds

  MAKE_KFOLDS(cohort_scan_ch)
  def train_files_ch = MAKE_KFOLDS.out.train_ids.flatten()
  def test_files_ch  = MAKE_KFOLDS.out.test_ids.flatten()

  def foldTag = { p -> def m=(p.toString() =~ /fold_?(\d+)/); m ? m[0][1] : 'NA' }
  def train_keyed_ch = train_files_ch.map { f -> tuple(foldTag(f), f) }
  def test_keyed_ch  = test_files_ch .map { f -> tuple(foldTag(f), f) }

  MAKE_LABELS_FOR_IDS_TRAIN(train_keyed_ch, cohort_scan_ch.first(), cohort_label_ch.first())
  def labels_train_ch = MAKE_LABELS_FOR_IDS_TRAIN.out.labels_train
  MAKE_LABELS_FOR_IDS_TEST (test_keyed_ch,  cohort_scan_ch.first(), cohort_label_ch.first())
  def labels_test_ch  = MAKE_LABELS_FOR_IDS_TEST.out.labels_test

  def train_with_labels_ch = train_keyed_ch.join(labels_train_ch)
  def test_with_labels_ch  = test_keyed_ch .join(labels_test_ch)

  /* GENUS-level training per fold */
  def genera_ch = Channel.fromList(params.genera as List)
  def genus_fold_joined = genera_ch
    .combine(train_with_labels_ch)
    .map { genus, fold, tr, lab -> tuple(genus as String, tr, lab, fold) }

  GENESIS_GENUS_SINGLE_FOLD(
    gds_seq_ch.first(), scan_lu_ch.first(),
    ch_resid, ch_meta.first(), ch_grm.first(), dam_rds_ch.first(), cage_rds_ch.first(),
    loco_rds_ch.first(), loco_bins_ch.first(), loco_Ns_ch.first(), loco_ids_ch.first(),
    Channel.value(params.genus_trait_prefix),
    genus_fold_joined
  )
  
  def assoc_genus_fold_tuples = GENESIS_GENUS_SINGLE_FOLD.out.assoc_dir_fold.map { d ->
    def mFold = (d.toString() =~ /fold_(\d+)/)
    def mGen  = (d.toString() =~ /assoc_genus_full__([^\/]+)$/)
    tuple(d, mFold ? mFold[0][1] : 'NA', mGen ? mGen[0][1] : 'NA')
  }

  /* --------------- Sentinels per (genus × fold) -------------------- */
  SELECT_SENTINELS_GENUS_FOLD(
    assoc_genus_fold_tuples,
    kept_plink_ch.first(), gene_loc_ch.first(),
    Channel.value(params.genus_trait_prefix)
  )

  // Group sentinels by genus (reused later)
  def sent_grouped_for_merge = SELECT_SENTINELS_GENUS_FOLD.out.sentinels_fold
    .map { s ->
      def mG = (s.toString() =~ /sentinels_genus__(.+?)__/); def gs = mG ? mG[0][1] : 'NA'
      tuple(gs, s)
    }
    .groupTuple()
  
  MERGE_SENTINELS_GENUS_FOLDS(sent_grouped_for_merge)
  def sent_cf_one_ch = MERGE_SENTINELS_GENUS_FOLDS.out.sentinels_cf_one
  
  /* --------------- OOF species betas on TEST ----------------------- */
  def sent_by_fold = SELECT_SENTINELS_GENUS_FOLD.out.sentinels_fold.map { s ->
    def mG = (s.toString() =~ /sentinels_genus__(.+?)__/); def gs = mG ? mG[0][1] : 'NA'
    def mF = (s.toString() =~ /fold(\d+)/);                def fd = mF ? mF[0][1] : 'NA'
    tuple(fd, gs, s)
  }
  def sent_with_test = sent_by_fold.join(test_with_labels_ch)
  
    SPECIES_BETAS_FOR_SENTINELS_GENUS_FOLD(
    gds_seq_ch.first(), scan_lu_ch.first(),
    ch_resid.first(), ch_meta.first(), ch_grm.first(), dam_rds_ch.first(), cage_rds_ch.first(),
    loco_rds_ch.first(), loco_bins_ch.first(), loco_Ns_ch.first(), loco_ids_ch.first(),
    Channel.value(params.species_trait_prefix),
    sent_with_test.map { fd, gs, sfile, test_ids, test_labels -> tuple(gs, fd, sfile, test_ids, test_labels) }
  )
    
  /* --------------- META across folds (per genus) ------------------- */
  def beta_cf_per_genus = SPECIES_BETAS_FOR_SENTINELS_GENUS_FOLD.out.beta_cf_fold
    .map { f -> def m=(f.toString() =~ /from_genus__(.+?)__/); tuple(m ? m[0][1] : 'NA', f) }
    .groupTuple()
  
  def meta_inputs = beta_cf_per_genus.join(sent_grouped_for_merge)
  META_BETA_CF_FOLDS(meta_inputs)
  def beta_cf_crossfit_one = META_BETA_CF_FOLDS.out.beta_cf_crossfit
  
  def beta_cf_keyed = META_BETA_CF_FOLDS.out.beta_cf_crossfit.map { f ->
    def m = (f.toString() =~ /crossfit__(.+)\.tsv$/)
    def gensan = m ? m[0][1] : 'NA'
    tuple(gensan, f)
  }
  
  def corr_main_ch = null
  if (params.trait=="Taxa"){
    
    /* --------------- PGLS per genus on cross-fit β ------------------- */
    CORRELATE_H2_PGLS_GENUS_TAGGED(
      beta_cf_keyed,
      ch_herit.first(), ch_tax.first(), ch_tree.first()
    )
    
    // Pick the main correlation file per genus
    corr_main_ch = CORRELATE_H2_PGLS_GENUS_TAGGED.out.corr_tsvs_genus.map { x ->
      def xs = (x instanceof List) ? x : [x]
      xs.find { !it.toString().endsWith('__EXCLUDED.tsv') } ?: xs.first()
    }

  } else  {
    PSEUDO_CORRELATE(
      beta_cf_keyed.map { gensan, f -> tuple(gensan, f) }   // (val, path)
    )
    corr_main_ch = PSEUDO_CORRELATE.out.corr_tsvs_genus
  }
    
  // Attach genus to each correlation file
  def corr_keyed = corr_main_ch.map { f ->
    def m = (f.toString() =~ /correlations_by_gene__([^\/]+)\.tsv$/)
    def gensan = m ? m[0][1] : 'NA'
    def genus_nm = gensan.replaceAll('_', ' ')
    tuple(genus_nm, gensan, f)
  }
  
  /* ---------- Run the simple synthesize step ---------- */
  SYNTHESIZE_RESULTS_SIMPLE(
    gds_seq_ch.first(),
    scan_lu_ch.first(),
    MAKE_GENESIS_COHORT_IDS.out.scan_ids.first(),
    MAKE_GENESIS_COHORT_IDS.out.label_ids.first(),
    corr_keyed,
    Channel.fromPath(params.genes_xlsx).first()
  )
   // Already keyed by genus_san from the process output
  def selected_genes_keyed = SYNTHESIZE_RESULTS_SIMPLE.out.selected_rdata

   /* ---------- Phase 9: Mediation (parallel folds) ---------- */
  if (params.mediation_enable) {
  
    if (params.trait=="Taxa"){
       // 9.0 – PGLS keep (per genus)
      PGLS_KEEP_VARIANTS(corr_keyed)
      def pgls_keep_keyed = PGLS_KEEP_VARIANTS.out.pgls_keep     // (genus_san, keep.tsv)
      
      // Sentinels keyed by (genus_san, fold)
      def sent_by_fold_keyed = SELECT_SENTINELS_GENUS_FOLD.out.sentinels_fold.map { s ->
        def mG = (s.toString() =~ /sentinels_genus__(.+?)__/); def gensan = mG ? mG[0][1] : 'NA'
        def mF = (s.toString() =~ /fold(\d+)/);                def fold   = mF ? mF[0][1] : 'NA'
        tuple(gensan, fold, s)
      }
      
      /* 9.1 – BH per fold AFTER PGLS filtering
         CHANGE: .join(...)  ->  .combine(..., by: 0)   (broadcast per-genus keep to all folds) */
      def bh_inputs = sent_by_fold_keyed
                        .combine(pgls_keep_keyed, by: 0)    // (genus_san, fold, sfile, keep.tsv)
      FILTER_SENTINELS_BH_PER_FOLD( bh_inputs )
      def bh_kept = FILTER_SENTINELS_BH_PER_FOLD.out.bh_kept // (genus_san, fold, bh_keep.tsv)
      
      /* 9.2 – Intersect BH-kept with the SAME PGLS keep (per fold)
         CHANGE: .join(...)  ->  .combine(..., by: 0) */
      def gate_inputs = bh_kept
                          .combine(pgls_keep_keyed, by: 0)  // (genus_san, fold, bh_keep.tsv, keep.tsv)
      INTERSECT_GATES_PER_FOLD( gate_inputs )
      // (genus_san, fold, snps_for_med.tsv)
    }
    
    // 9.3 Build mediation inputs (directional Bacteroides<->Prevotella)
    def snps_for_med =null
    if (params.trait=="Taxa"){
      snps_for_med = INTERSECT_GATES_PER_FOLD.out.snps_for_med
    } else {
      snps_for_med = SELECT_SENTINELS_GENUS_FOLD.out.snps_for_med
    }
    
    def med_inputs =
      snps_for_med
        .map { gensan, fold, snps -> tuple(gensan, fold, snps) }
        .combine(selected_genes_keyed, by: 0)               // (genus_san, fold, snps, selected_genes.RData)
        .map { gensan, fold, snps, selRdata ->
          def genus_name = gensan.replaceAll('_',' ')
          def other = (genus_name == 'Bacteroides') ? 'Prevotella'
                     : (genus_name == 'Prevotella') ? 'Bacteroides' : 'Prevotella'
          tuple(gensan, fold, snps, selRdata, other)
        }
  
    // Shared/broadcast inputs (value channels are fine and reusable)  — see docs on value channels
    def resid_rda_file    = Channel.fromPath(params.resid_rda ?: params.resid_rda).first()
    def beta_rda_file     = Channel.fromPath(params.beta_rda ?: params.resid_rda).first()
    def alpha_rda_file     = Channel.fromPath(params.alpha_rda ?: params.resid_rda).first()
    def clusters_rda_file = Channel.fromPath(params.clusters_rda ?: params.resid_rda).first()
    def phen_rda_file     = Channel.fromPath(params.med_glucose_rdata).first()
    def dam_rds_file      = BUILD_SRM.out.dam_rds.first()
    def cage_rds_file     = BUILD_SRM.out.cage_rds.first()
    def loco_rds_file     = PACK_LOCO_RDS.out.loco_rds.first()
    def phen_glucose_col  = Channel.value(params.med_glucose_column)
    def phen_bmi_col      = Channel.value(params.med_bmi_column)
    def trait_prefix_val  = Channel.value(params.genus_trait_prefix)
    def beta_trait_val    = Channel.value(params.med_beta_trait_name)
    def alpha_trait_val    = Channel.value(params.med_alpha_trait_name)
    def guild_trait_val    = Channel.value(params.med_guild_trait_name)
    def outcomes          = Channel.value(params.outcomes)
  
    // 9.3 – Run mediation per (genus × fold) fully in parallel
    RUN_MEDIATION_FOLD(
      med_inputs,
      scan_lu_ch.first(),
      resid_rda_file, beta_rda_file, alpha_rda_file, clusters_rda_file,
      phen_rda_file, phen_glucose_col, phen_bmi_col,
      dam_rds_file, cage_rds_file,
      loco_rds_file, loco_bins_ch.first(), loco_Ns_ch.first(), loco_ids_ch.first(),
      trait_prefix_val, beta_trait_val, alpha_trait_val, guild_trait_val, outcomes,
      Channel.value(params.mr_enable), Channel.value(params.mr_fstat_min)
    )
  
    // 9.4 – Meta across folds per genus (wait for exactly k folds)
    def med_by_genus = RUN_MEDIATION_FOLD.out.med_results_fold
      .map { gensan, fold, f -> tuple(gensan, f) }
      .groupTuple(size: params.kfolds)  // emit once per genus after 5 folds arrive
    META_MEDIATION_ACROSS_FOLDS( med_by_genus )
    
    SIG_IV(META_MEDIATION_ACROSS_FOLDS.out.med_meta)
    SIG_MVMR(META_MEDIATION_ACROSS_FOLDS.out.med_meta)
  }

  /* ============================= RG / CVBLUP ============================== */
  rg_enable=false
  he_enable=false
  cvblup_enable=false
  if (rg_enable || he_enable || cvblup_enable) {

    def focus_genus = (params.genera as List)[0]
    def gensan      = focus_genus.replaceAll(' ','_')

    // One GRM for both RG and BLUP
    MAKE_GRM_ALL( kept_plink_ch )

    BUILD_GCTA_PHENOS(
      kept_plink_ch,
      H5_TO_GDS_CAND.out.scan_lookup.first(),
      Channel.fromPath(params.resid_rda).first(),
      Channel.fromPath(params.beta_rda).first(),
      Channel.fromPath(params.alpha_rda).first(),
      Channel.fromPath(params.clusters_rda).first()
    )

    MAKE_GENESET_EXTRACTS_FOCUS(
      SYNTHESIZE_RESULTS_SIMPLE.out.selected_rdata,
      kept_plink_ch
    )

    MAKE_GRMS_GENESET_VS_BG(
      kept_plink_ch,
      MAKE_GENESET_EXTRACTS_FOCUS.out.snp_focus,
      MAKE_GENESET_EXTRACTS_FOCUS.out.snp_bg
    )

    RUN_GCTA_BIVAR_RG(
      mgrm_with_srm_ch.first(),
      BUILD_GCTA_PHENOS.out.pheno.first(),
      BUILD_GCTA_PHENOS.out.trait_map.first(),
      BUILD_GCTA_PHENOS.out.pairs
    )

  //  RUN_GCTA_BIVAR_RG_MGRM(
  //    MAKE_GRMS_GENESET_VS_BG.out.grm_focus_gz,       // focus (gz + id) as before
  //    MAKE_GRMS_GENESET_VS_BG.out.grm_bg_gz,          // background (gz + id) as before
  //    BUILD_GCTA_PHENOS.out.pheno.first(),
  //    BUILD_GCTA_PHENOS.out.pairs,
  //    tuple( BUILD_MGRM_WITH_SRM.out.srm_dam_harmo_gz, 
  //           BUILD_MGRM_WITH_SRM.out.srm_dam_harmo_id),
  //    tuple( BUILD_MGRM_WITH_SRM.out.srm_cage_harmo_gz, 
  //           BUILD_MGRM_WITH_SRM.out.srm_cage_harmo_id)
  //  )
  //  
  //  
  //  RUN_GCTA_HE_BIVAR(
  //    mgrm_with_srm_ch.first(),
  //    BUILD_GCTA_PHENOS.out.pheno.first(),
  //    BUILD_GCTA_PHENOS.out.pairs
  //  )
  //  
  //  RUN_GCTA_CVBLUP(
  //    mgrm_with_srm_ch.first(),                 // mgrm_with_srm.txt
  //    BUILD_GCTA_PHENOS.out.pheno.first(),      // gcta.pheno
  //    BUILD_GCTA_PHENOS.out.trait_map.first()   // trait_map.tsv
  //  )
  //  
  //  
  //  PARSE_RG_OUTPUTS(
  //    RUN_GCTA_BIVAR_RG.out.hsq_all,
  //    RUN_GCTA_BIVAR_RG_MGRM.out.hsq_mgrm,
  //    RUN_GCTA_HE_BIVAR.out.he_logs
  //  )
  //}
    
  blup_enable=false
  /* ===== BLUP + alignment (optional) ===== */
  if (blup_enable) {
    // Reuse the same all-SNP GRM
    def genus_one = (params.genera as List)[0] as String
    MAKE_GCTA_PHENO_GENUS(
      H5_TO_GDS_CAND.out.scan_lookup.first(),
      ch_resid,
      Channel.value(genus_one),
      Channel.value(params.genus_trait_prefix)
    )
    
    RUN_GCTA_BLUP_GENUS(
      mgrm_with_srm_ch.first(),
      MAKE_GCTA_PHENO_GENUS.out.pheno.first()
    )
    
    //if (params.cvblup_align_enable) {
    //  CVBLUP_ALIGN_SUMMARY(
    //    RUN_GCTA_BLUP_GENUS.out.indi_blp.first(),
    //    H5_TO_GDS_CAND.out.scan_lookup.first(),
    //    Channel.fromPath(params.beta_rda).first(),
    //    Channel.fromPath(params.alpha_rda).first(),
    //    Channel.fromPath(params.clusters_rda).first()
    //  )
    }
  }
}

workflow.onComplete {
  println "Pipeline completed!"
  println "Started at  ${workflow.start}"
  println "Finished at ${workflow.complete}"
  println "Elapsed:    ${workflow.duration}"
  println "Status:     ${ workflow.success ? 'OK' : 'FAILED' }"
}
