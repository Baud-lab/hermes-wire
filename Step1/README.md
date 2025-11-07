🧬 HERMES-WIRE — Step 1: Data Pre-processing and Microbiome Profiling
Centre for Genomic Regulation (CRG), Barcelona — 2025
Authors: Felipe Morillo Sanz Dias, Luca Cozzuto, Hélène Tonnele, Amelie Baud
📘 Overview
HERMES-WIRE Step 1 is the first stage of the HERMES-WIRE (HERitable MicrobiomE Structure) workflow, designed to perform read-level pre-processing and microbiome profiling from shallow-shotgun or metagenomic sequencing data.
It provides an automated and reproducible Nextflow pipeline for:
Pre-processing of raw reads (QC, trimming, host read removal).
Taxonomic profiling using Kaiju and a GTDB-based protein index.
Functional profiling using PRROMenade and the IBM Functional Genomics Platform (IFGP) reference catalogue.
This step produces non-host filtered reads, taxonomic abundance tables, functional EC tables, and summary reports (MultiQC).
⚙️ Pipeline structure
Step1/
├── bin/                     # Auxiliary scripts and binaries
├── cmd                      # Command template
├── config.yaml              # MultiQC configuration file
├── dataset/                 # Example raw datasets
├── dockerfiles/             # Singularity/Docker recipes
├── input/                   # Input reference indices, adapters, host genome
├── local_modules.nf         # Locally defined Nextflow modules
├── main.nf                  # Main Nextflow pipeline
├── nextflow.config          # Default Nextflow configuration
├── params.yaml              # User-defined parameters
├── singularity_containers/  # Singularity images cache
├── tool_opt.tsv             # Tool-specific parameter definitions
├── submit_nf.sh             # Example SLURM submission script
└── work/                    # Nextflow working directory (auto-generated)
🚀 Quick start
1. Install dependencies
You’ll need:
Nextflow ≥ 23.10
Singularity or Apptainer (for container execution)
HPC or Unix environment with SLURM (optional)
# Load Nextflow and Singularity
module load Nextflow
module load Singularity
2. Download required reference indices
The pipeline depends on two large pre-computed indices hosted on Zenodo:
Profiling type	Description	Zenodo link
Taxonomic	GTDB R207 protein catalogue for Kaiju (gtdb_index.fmi)	https://zenodo.org/uploads/17545483
Functional	IFGP (IBM Functional Genomics Platform) index for PRROMenade (bactvirus2020*)	https://zenodo.org/uploads/17545032
Download and extract both archives into the input/indices/ directory:
cd input/indices/

# Taxonomic index
wget https://zenodo.org/uploads/17545483 -O GTDB207_kaiju_index.tar.gz
tar -xvzf GTDB207_kaiju_index.tar.gz -C ./kaiju/

# Functional index
wget https://zenodo.org/uploads/17545032 -O IFGP_prromenade_index.tar.gz
tar -xvzf IFGP_prromenade_index.tar.gz -C ./prromenade/
3. Configure parameters
Edit params.yaml to specify your dataset paths and options.
Example:
module_preprocessing: "YES"
module_taxonomic_mapping: "YES"
module_functional_mapping: "YES"

tool_opt: "./tool_opt.tsv"

preprocessing_file_raw_reads: "./dataset_raw/demo/*_R{1,2}_001.fastq.gz"
preprocessing_file_adapters: "./input/Adapters_Qiita.fa"
preprocessing_file_host_genome: "./input/mRatBN7.2.fna"

taxonomic_file_microbiome_index: "./input/indices/kaiju/gtdb_index.fmi"
taxonomic_file_genome_sizes: "./input/indices/kaiju/catalogue_stats.csv"
taxonomic_parameter_normalization_by_genome_size: "YES"

functional_file_index: "./input/indices/prromenade/bactvirus2020*"
functional_file_ec: "./input/indices/prromenade/taxid_name.txt"

output_folder: "./Results_Step1/"
4. Run the pipeline
Execute:
nextflow run main.nf -params-file params.yaml -profile singularity
or via SLURM using the provided submission script:
sbatch submit_nf.sh
🧩 Modules
Module	Description	Tools used
Pre-processing	Quality control, trimming, host read removal	FastQC, Trimmomatic, Bowtie2, Samtools
Taxonomic profiling	Protein-level classification against GTDB index	Kaiju, custom parsers
Functional profiling	Functional annotation via PRROMenade	PRROMenade, IFGP index
Reporting	QC summary and statistics	MultiQC, custom R scripts
📂 Output structure
Results are organized as follows:
Results_Step1/
├── Non_Host/                 # Reads without host contamination
├── Counts_Taxonomic/         # Kaiju count tables (taxonomic)
├── Counts_Functional/        # PRROMenade EC tables (functional)
├── Reports_Taxonomic/        # MultiQC and summary reports (taxonomic)
├── Reports_Functional/       # MultiQC and summary reports (functional)
└── pipeline_info/            # Logs and Nextflow reports
🧠 Citation
If you use this workflow or its pre-computed indices, please cite:
Morillo FMSD, Cozzuto L, Tonnele H, Baud A (2025).
HERMES-WIRE: HERitable MicrobiomE Structure — Workflow for Interpreting host–microbiome Relationships & Effects.
Centre for Genomic Regulation (CRG), Barcelona.
https://github.com/Baud-lab/hermes-wire
🧩 Acknowledgements
This pipeline was developed under the Baud Lab at the Centre for Genomic Regulation (CRG) and Universitat Pompeu Fabra (UPF).
We acknowledge the support of the HPC Core Facility at CRG.
🧾 License
© 2025 Centre for Genomic Regulation (CRG) and the authors.
Distributed under the MIT License.
