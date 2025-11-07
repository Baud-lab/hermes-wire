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

module preprocessing          : ${params.module_preprocessing}
module taxonomic mapping      : ${params.module_taxonomic_mapping}
module functional mapping     : ${params.module_functional_mapping}

----------------------------------------------------------------------

2. TOOL OPTIONS

tool_opt                      : ${params.tool_opt}

----------------------------------------------------------------------

3. MODULE: PRE-PROCESSING

raw reads                     : ${params.preprocessing_file_raw_reads}
adapters                      : ${params.preprocessing_file_adapters}
host genome                   : ${params.preprocessing_file_host_genome}

----------------------------------------------------------------------

4. NON-HOST reads

non-host reads                : ${params.file_non_host_reads}

----------------------------------------------------------------------

5. MODULE: TAXONOMIC PROFILING

microbiome index              : ${params.taxonomic_file_microbiome_index}
genome sizes                  : ${params.taxonomic_file_genome_sizes}
normalization by genome size  : ${params.taxonomic_parameter_normalization_by_genome_size}

----------------------------------------------------------------------

6. MODULE: FUNCTIONAL PROFILING

functional catalogue          : ${params.functional_file_ref_catalogue}
functional index              : ${params.functional_file_index}
ec numbers                    : ${params.functional_file_ec}

----------------------------------------------------------------------

7. OUTPUT FOLDER
output folder                 : ${params.output_folder}

----------------------------------------------------------------------
"""

if (params.help) {
    log.info 'This is the MicroBiome Profiles\'s pipeline'
    log.info 'Enjoy it'
    log.info '\n'
    exit 1
}

/*
* Define output folders
*/
outputNonHost    = "${params.output_folder}/Non_Host"
outputCounts_Taxonomic     = "${params.output_folder}/Counts_Taxonomic"
outputCounts_Functional     = "${params.output_folder}/Counts_Functional"
outputReports_Taxonomic  = "${params.output_folder}/Reports_Taxonomic"
outputReports_Functional  = "${params.output_folder}/Reports_Functional"

/*
* Define modules
*/
def local_modules_file = file("${projectDir}/local_modules.nf")
def subworkflowsDir = "${projectDir}/../BioNextflow/subworkflows"

// Getting comand line options from a file
include { getParameters } from "${local_modules_file}" 
progPars = getParameters(params.tool_opt)

// Include processes from BioNextflow
include { FASTQCP as FASTQC_RAW; FASTQCP as FASTQC_TRIM} from "${subworkflowsDir}/qc/fastqc"
include { FILTERADAPTER as TRIMMOMATIC } from "${subworkflowsDir}/trimming/trimmomatic" addParams(EXTRAPARS: progPars["trimming--trimmomatic"], OUTPUT: "", LABEL: "big_cpus")
include { ALL as BOWTIE2_ALL } from "${subworkflowsDir}/alignment/bowtie2" addParams(LABEL: "big_cpus")
include { SORT; INDEX  } from "${subworkflowsDir}/misc/samtools"
include { FASTQ_PAIRS  } from "${subworkflowsDir}/misc/samtools" addParams(EXTRAPARS: progPars["remove_host--samtools"], OUTPUT: outputNonHost)
include { ALN as KAIJU_ALN  } from "${subworkflowsDir}/metagenomics/kaiju" addParams(LABEL: "big_mem_cpus", EXTRAPARS: progPars["alignment_microbiome--kaiju"])
include { REPORT as MULTIQC_TAX } from "${subworkflowsDir}/reporting/multiqc" addParams(EXTRAPARS: progPars["report--multiqc"], OUTPUT: outputReports_Taxonomic )
include { REPORT as MULTIQC_FUNC } from "${subworkflowsDir}/reporting/multiqc" addParams(EXTRAPARS: progPars["report--multiqc"], OUTPUT: outputReports_Functional )
// Include processes from local modules
//include { humann_profiling } from "${local_modules_file}" addParams(EXTRAPARS: progPars["alignment_microbiome--humann"])
include { prromenade_index; prromenade_classifier } from "${local_modules_file}"
include { makeKaijuTable as makeKaijuTable_tax; parse_kaiju_norm; parse_kaiju as parse_kaiju_tax; makeKaiTable as makeKaiTable_tax} from "${local_modules_file}" addParams(OUTPUT: outputCounts_Taxonomic )
include { get_stats as get_stats_tax} from "${local_modules_file}" addParams(OUTPUT: outputReports_Taxonomic )
include { makeKaijuTable as makeKaijuTable_func; parse_kaiju as parse_kaiju_func; makeKaiTable as makeKaiTable_func} from "${local_modules_file}" addParams(OUTPUT: "" )
include { makeECTables} from "${local_modules_file}" addParams(OUTPUT: outputCounts_Functional )
include { get_stats as get_stats_func} from "${local_modules_file}" addParams(OUTPUT: outputReports_Functional )

/*
* Getting data from parameters and storing into a channel
*/

//Modules
module_preprocessing = params.module_preprocessing
if (module_preprocessing!="NO" && module_preprocessing!="YES") {
  println "WARNING: The options for 'Module Preprocessing' should be 'YES' or 'NO'. By default, if none of these options is informed it is treated as 'NO'."
  module_preprocessing="NO"
}
module_taxonomic_mapping = params.module_taxonomic_mapping
if (module_taxonomic_mapping!="NO" && module_taxonomic_mapping!="YES") {
  println "WARNING: The options for 'Module Taxonomic Mapping' should be 'YES' or 'NO'. By default, if none of these options is informed it is treated as 'NO'."
  module_taxonomic_mapping="NO"
}
module_functional_mapping = params.module_functional_mapping
if (module_functional_mapping!="NO" && module_functional_mapping!="YES") {
  println "WARNING: The options for 'Module Functional Mapping' should be 'YES' or 'NO'. By default, if none of these options is informed it is treated as 'NO'."
  module_functional_mapping="NO"
}

// Pre-processing
fastq_files = params.preprocessing_file_raw_reads
adapters_file = params.preprocessing_file_adapters
host_genome = params.preprocessing_file_host_genome
multiqc_cfile = Channel.fromPath( "${projectDir}/config.yaml" )
// Non-host reads
non_host_reads = params.file_non_host_reads
// Taxonomic Profiling
microbiome_index = params.taxonomic_file_microbiome_index
genome_sizes = Channel.fromPath(params.taxonomic_file_genome_sizes)
norm_genome = params.taxonomic_parameter_normalization_by_genome_size
// Functional profiling
//functional_index = params.functional_file_index_folder
//functional_utility = params.functional_file_utility_folder
//functional_classification = params.functional_params_classification_types
functional_catalogue = params.functional_file_ref_catalogue
functional_index = params.functional_file_index
ec_numbers = params.functional_file_ec
//Empty variables
fq_raws=Channel.empty()
pairs_and_not=Channel.empty()
trimmed_fq=Channel.empty()
fq_trim=Channel.empty()
out=Channel.empty()
not_host_fq=Channel.empty()
ktable=Channel.empty()


workflow {
	
	if (module_preprocessing =="YES") {
	  	// 1. Pre-processing:
    	// 1.1 QC and Trimming
    	if (fastq_files=="") {
    	  error "ERROR: As module 'Pre-processing' was selected you need to inform the path to the raw reads.\n"
    	} else {
    	  if (adapters_file=="") {
    	    error "ERROR: As module 'Pre-processing' was selected you need to inform the file with the adapters sequences to be trimmed out.\n"
    	  } else {
        	fastq_files = Channel.fromFilePairs( params.preprocessing_file_raw_reads)
        	adapters_file = Channel.fromPath( params.preprocessing_file_adapters)
        	// 1.1.1 Quality check (QC) on raw reads
        	//fq_raws = FASTQC_RAW(fastq_files)
        	// 1.1.2 Trimming raw reads
        	trimmed_fq = TRIMMOMATIC(fastq_files, adapters_file)
        	// 1.1.3 Quality check (QC) on trimmed reads
        	//fq_trim = FASTQC_TRIM(trimmed_fq.trimmed_reads)
        	}
    	  }
    	// 1.2 Removal of host reads
    	// 1.2.1 Aligment with the host genome
    	pairs_and_not = trimmed_fq.trimmed_reads.join(trimmed_fq.trimmed_ureads).map{
    		[it[0], [ it[1][0], it[1][1], it[2][0], it[2][1] ] ] }
    	if (host_genome=="") {
    	  error "ERROR: As module 'Pre-processing' was selected you have to inform what fasta file should be used to build the index for the host genome."
    	} else {
    	  out = BOWTIE2_ALL(host_genome, pairs_and_not)
    	}

    	// 1.2.2 Sorting aligments
    	sorted_alns = SORT(out.aln)
    	// 1.2.3 New fastq files without host reads
    	not_host_fq = FASTQ_PAIRS(sorted_alns)
  }

  // 2. Taxonomic mapping
  if (module_taxonomic_mapping =="YES") {
    if (module_preprocessing =="NO" && non_host_reads=="") {
      error "ERROR: When module 'Pre-processing' is not selected and module 'Taxonomic Mapping' is selected you have to provide the path to the fastq files without host reads.\n"
    } else {
      if (module_preprocessing =="NO" && non_host_reads!="") {
        not_host_fq = Channel.fromFilePairs(params.file_non_host_reads).map{[
          it[0],
          it[1][0],
          it[1][1]
        ]}
      }
    }
    // 2.1 Alignment with the microbiome reference catalogue
    if (microbiome_index=="") {
      error "ERROR: As module 'Taxonomic Mapping' was selected you have to inform what index should be used.\n"
    } else {
      microbiome_index = Channel.fromPath(params.taxonomic_file_microbiome_index).first()
    }
    kaiju_out=KAIJU_ALN(not_host_fq,microbiome_index)
    // 2.2 Parsing read counts with adjusted counts (Unique + Multiple )
    if (norm_genome == "YES") {
	  	prof="Taxonomic"
	  	kaiju_res = parse_kaiju_norm(kaiju_out.combine(genome_sizes),prof)
	  } else {
	  	kaiju_res = parse_kaiju_tax(kaiju_out,prof)
	  }	
	  // 2.3 Make final counts table
	  kaitable = makeKaiTable_tax(kaiju_res.kai.map{it[1]}.collect(), norm_genome)	
	  if (module_preprocessing == "YES") {
	    // 2.4 Make table with mapping statistics
	    ktable = makeKaijuTable_tax(kaiju_res.kstats.map{it[1]}.collect())	
	    // 2.5 Make MULTIQC report
	    multiqc_files = multiqc_cfile.mix(ktable).mix(trimmed_fq.trim_log.map{it[1]}).mix(fq_raws).mix(fq_trim).mix(out.logs.map{it[1]}).collect()
	    multi_data = MULTIQC_TAX(multiqc_files).data
	    get_stats_tax(multi_data)	    
	  }
  }
	
	//// 3. Functional mapping
	//if (module_functional_mapping =="YES"){
	//  if (module_preprocessing =="NO" && non_host_reads=="") {
  //    error "ERROR: When module 'Pre-processing' is not selected and module 'Functional Mapping' is selected you have to provide the path to the fastq files without host reads.\n"
  //  } else {
  //    if (module_preprocessing =="NO" && non_host_reads!="") {
  //      not_host_fq = Channel.fromFilePairs(params.file_non_host_reads).map{[
  //        it[0],
  //        it[1][0],
  //        it[1][1]
  //      ]}
  //    }
  //  }
	//  // 3.1 Alignment with the microbiome reference catalogue
	//  if (functional_classification=="") {
  //  	    error "ERROR: Inform the functional classification types.\n"
  //  } else {
  //    functional_classification = Channel.fromPath(params.functional_params_classification_types)
  //  }
  //  if (functional_index=="") {
  //  	    error "ERROR: Inform the UniRef90 index folder.\n"
  //  } else {
  //    functional_index = Channel.fromPath(params.functional_file_index_folder)
  //  }
  //  if (functional_utility=="") {
  //  	    error "ERROR: Inform the utility table for functional regrouping of UniRef90 IDs into functional IDs.\n"
  //  } else {
  //    functional_utility = Channel.fromPath(params.functional_file_utility_folder)
  //  }
  //  humann_input=not_host_fq.combine(functional_index).combine(functional_utility).combine(functional_classification)
  //  kaiju_out=humann_profiling(humann_input)
	//}
	
	// 3. Functional mapping
	if (module_functional_mapping =="YES"){
	  if (module_preprocessing =="NO" && non_host_reads=="") {
      error "ERROR: When module 'Pre-processing' is not selected and module 'Functional Mapping' is selected you have to provide the path to the fastq files without host reads.\n"
    } else {
      if (module_preprocessing =="NO" && non_host_reads!="") {
        not_host_fq = Channel.fromFilePairs( params.file_non_host_reads)
      }
    }
	  // 3.1 Alignment with the reference catalogue for functional assignment of reads
	      if (functional_index=="") {
    	  if (functional_catalogue=="") {
    	    error "ERROR: As module 'Functional Mapping' was selected and no index was provided for the functional reference catalogue you have to inform what fasta file should be used to build the index.\n"
    	  } else {
    	    functional_catalogue = Channel.fromPath( params.functional_file_catalogue)
    	    functional_index = prromenade_index(functional_catalogue)
    	  }
    } else {
      functional_index = Channel.fromPath(params.functional_file_index)
    }
    prromenade_indexes = functional_index.collect()
	  if (ec_numbers=="") {
    	    error "ERROR: As module 'Functional Mapping' was selected a list with EC Numbers has to be informed.\n"
    } else {
      ec_numbers = Channel.fromPath(params.functional_file_ec)
    }
    prromenade_input=not_host_fq.combine(prromenade_indexes.toList()).combine(ec_numbers)
    kaiju_out=prromenade_classifier(prromenade_input)
	  // 3.2 Parsing read counts with adjusted counts (Unique + Multiple )
	  prof="Functional"
	  kaiju_res = parse_kaiju_func(kaiju_out,prof)
	  // 3.3 Make final counts table
	  norm_genome="NO"
	  kaitable = makeKaiTable_func(kaiju_res.kai.map{it[1]}.collect(), norm_genome)
	  makeECTables(kaitable[0])
	  if (module_preprocessing =="YES") {
	    // 3.4 Make table with mapping statistics
	    ktable = makeKaijuTable_func(kaiju_res.kstats.map{it[1]}.collect())	
	    // 3.5 Make MULTIQC report
	    multiqc_files = multiqc_cfile.mix(ktable).mix(trimmed_fq.trim_log.map{it[1]}).mix(fq_raws).mix(fq_trim).mix(out.logs.map{it[1]}).collect()
	    multi_data = MULTIQC_FUNC(multiqc_files).data
	    get_stats_tax(multi_data)
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
