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

/* -------------------- Params -------------------- */
params.help   = params.help   ?: false
params.resume = params.resume ?: false

params.outdir              = params.outdir              ?: './Results'
params.input_faa_glob      = params.input_faa_glob      ?: 'input/*.faa'
params.emapper_mode        = params.emapper_mode        ?: 'diamond'       // or 'mmseqs'
params.emapper_args_extra  = params.emapper_args_extra  ?: ''              // e.g. '--itype proteins'
params.threads             = params.threads             ?: 8

// Reference files
params.taxonomy_tsv        = params.taxonomy_tsv        ?: null            // GTDB mapping (Genome\tTaxonomy)
params.phylo_tree          = params.phylo_tree          ?: null            // GTDB R207 newick
params.herit_rdata         = params.herit_rdata         ?: null            // contains all_VCs_full
params.cog_def_tab         = params.cog_def_tab         ?: null            // e.g., cog-24.def.tab

// Which function spaces to build (comma-separated)
params.functions           = params.functions           ?: 'eggNOG_OGs,CAZy,EC,KEGG_ko,KEGG_Module,KEGG_Pathway,PFAMs,COG_category,COG_pathway'

// PGLS facets
params.subset_id           = params.subset_id           ?: 'ALL'
params.method              = params.method              ?: 'Shallow'
params.report_model        = params.report_model        ?: 'ML'            // 'ML' or 'AIC'
params.branchlen_transform = params.branchlen_transform ?: 'none'          // 'none' or 'grafen'
params.groups              = params.groups              ?: 'g__Prevotella,g__Bacteroides'

// Containers (optional): leave empty to use host env
params.container_emapper   = params.container_emapper   ?: ''
params.container_R         = params.container_R         ?: ''

/* -------------------- Helper printing -------------------- */
if (params.help) {
  log.info """
Usage:
  nextflow run main.nf -params-file params.yaml

Key params:
  input_faa_glob      : ${params.input_faa_glob}
  taxonomy_tsv        : ${params.taxonomy_tsv}
  phylo_tree          : ${params.phylo_tree}
  herit_rdata         : ${params.herit_rdata}
  cog_def_tab         : ${params.cog_def_tab}
  functions           : ${params.functions}
  groups              : ${params.groups}
  outdir              : ${params.outdir}
"""
  System.exit(0)
}

/* -------------------- Channels -------------------- */
Channel
  .fromPath( params.input_faa_glob, type: 'file', checkIfExists: true )
  .ifEmpty { error "No FASTA files matched: ${params.input_faa_glob}" }
  .set { ch_faa }

/* -------------------- Processes -------------------- */

process EMMAPPER_ANNOTATE {
  tag { "${fasta.simpleName}" }
  cpus params.threads
  memory '8 GB'
  errorStrategy 'terminate'
  publishDir "${params.outdir}/01_annotations", mode: 'copy'

  container params.container_emapper ?: null

  input:
    path fasta

  output:
    path "${fasta.simpleName}.annotations.tsv"

  script:
  """
  set -euo pipefail
  emapper_run.sh \\
    --fasta "${fasta}" \\
    --mode  "${params.emapper_mode}" \\
    --threads ${task.cpus} \\
    --extra "${params.emapper_args_extra}"
  """
}

process BUILD_MATRICES {
  tag "make_func_matrices"
  cpus 2
  memory '8 GB'
  errorStrategy 'terminate'
  publishDir "${params.outdir}/02_matrices", mode: 'copy'

  container params.container_R ?: null

  input:
    path(tsvs)

  output:
    path "func_matrices.RData"
    path "func_eggNOG_OGs.RData"
    path "qc_functions_summary.tsv"

  script:
  """
  set -euo pipefail
  mkdir -p ann_dir
  cp ${tsvs} ann_dir/

  make_func_matrices.R \\
    --ann_dir ./ann_dir \\
    --functions "${params.functions}" \\
    --cog_def "${params.cog_def_tab}" \\
    --out_prefix .
  """
}

process ENRICHMENT_PGLS {
  tag "enrichment_pgls"
  cpus 2
  memory '8 GB'
  errorStrategy 'terminate'
  publishDir "${params.outdir}/03_enrichment", mode: 'copy'

  container params.container_R ?: null

  input:
    tuple path(func_mats_rdata), path(func_enog_rdata)

  output:
    path "enrichment_*.RData"
    path "res_*.csv"
    path "*.pdf"

  when:
    params.taxonomy_tsv && params.phylo_tree && params.herit_rdata

  script:
  """
  set -euo pipefail
  enrichment_analysis_phylogenetic.R \\
    --taxonomy "${params.taxonomy_tsv}" \\
    --phylo    "${params.phylo_tree}" \\
    --herit    "${params.herit_rdata}" \\
    --func_mats "${func_mats_rdata}" \\
    --func_enog "${func_enog_rdata}" \\
    --subset   "${params.subset_id}" \\
    --method   "${params.method}" \\
    --report   "${params.report_model}" \\
    --branchlen "${params.branchlen_transform}" \\
    --groups   "${params.groups}" \\
    --functions "${params.functions}" \\
    --outdir   "."
  """
}


/* -------------------- Workflow -------------------- */
workflow {
  // 1) annotate one-by-one
  ch_annot = EMMAPPER_ANNOTATE(ch_faa).out
  // 2) collect all annotation files into a single list for one BUILD_MATRICES run
  ch_annot_all = ch_annot.collect()
  // 3) BUILD_MATRICES emits tuple (func_mats, func_enog)
  mats = BUILD_MATRICES(ch_annot_all)
  // 4) feed that tuple straight into ENRICHMENT_PGLS
  ENRICHMENT_PGLS(mats)
}

