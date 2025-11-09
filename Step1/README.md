# 🧬 HERMES-WIRE — Step 1: Data Pre-processing and Microbiome Profiling

<b>Centre for Genomic Regulation (CRG), Barcelona — 2025</b>  
<b>Authors:</b> Felipe Morillo Sanz Dias, Luca Cozzuto, Hélène Tonnele, Amelie Baud  

---

## <b>📘 Overview</b>

<b>HERMES-WIRE Step 1</b> is the first stage of the <b>HERMES-WIRE</b> (<b>HER</b>itable <b>M</b>icrobiom<b>E</b> <b>S</b>tructure) workflow.  
It performs read-level pre-processing and microbiome profiling from shallow-shotgun or metagenomic sequencing data.

This automated and reproducible Nextflow pipeline executes:

1. <b>Pre-processing</b> of raw reads (QC, trimming, host-read removal).  
2. <b>Taxonomic profiling</b> using Kaiju and a GTDB-based protein index.  
3. <b>Functional profiling</b> using PRROMenade and the IBM Functional Genomics Platform (IFGP) reference catalogue.  

It produces non-host filtered reads, taxonomic abundance tables, functional EC tables, and MultiQC reports.

---

## <b>⚙️ Pipeline structure</b>

```

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

````

---

## <b>🚀 Quick start</b>

### <b>1. Install dependencies</b>

You’ll need:

- **Nextflow ≥ 23.10**  
- **Singularity or Apptainer** (for container execution)  
- **HPC or Unix environment** with SLURM (optional)

```bash
# Load Nextflow and Singularity
module load Nextflow
module load Singularity
```

---

### <b>2. Download required reference indices</b>

The pipeline depends on two large pre-computed indices hosted on <b>Zenodo</b>:

| Profiling type | Description                                                                     | Zenodo link                                                                |
| -------------- | ------------------------------------------------------------------------------- | -------------------------------------------------------------------------- |
| **Taxonomic**  | GTDB R207 protein catalogue for Kaiju (`gtdb_index.fmi`)                        | [https://zenodo.org/uploads/17545483](https://zenodo.org/uploads/17545483) |
| **Functional** | IFGP (IBM Functional Genomics Platform) index for PRROMenade (`bactvirus2020*`) | [https://zenodo.org/uploads/17545032](https://zenodo.org/uploads/17545032) |

Download and extract both archives into `input/indices/` using commands ```wget``` and ```tar -xvzf```, respectivelly.

---

### <b>3. Configure parameters</b>

Edit <b>`params.yaml`</b> to specify your dataset paths and options:

```yaml
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
```

---

### <b>4. Run the pipeline</b>

```bash
nextflow run main.nf -params-file params.yaml -profile singularity
```

Or submit through SLURM:

```bash
sbatch submit_nf.sh
```

---

## <b>🧩 Modules</b>

| Module                   | Description                                     | Tools used                             |
| ------------------------ | ----------------------------------------------- | -------------------------------------- |
| **Pre-processing**       | Quality control, trimming, host-read removal    | FastQC, Trimmomatic, Bowtie2, Samtools |
| **Taxonomic profiling**  | Protein-level classification against GTDB index | Kaiju, custom parsers                  |
| **Functional profiling** | Functional annotation via PRROMenade            | PRROMenade, IFGP index                 |
| **Reporting**            | QC summary and statistics                       | MultiQC, custom R scripts              |

---

## <b>📂 Output structure</b>

```
Results_Step1/
├── Non_Host/                 # Reads without host contamination
├── Counts_Taxonomic/         # Kaiju count tables (taxonomic)
├── Counts_Functional/        # PRROMenade EC tables (functional)
├── Reports_Taxonomic/        # MultiQC and summary reports (taxonomic)
├── Reports_Functional/       # MultiQC and summary reports (functional)
└── pipeline_info/            # Logs and Nextflow reports
```

---

## <b>🧠 Citation</b>

If you use this workflow or its pre-computed indices, please cite:

> Morillo FMSD, Cozzuto L, Tonnele H, Baud A (2025). <i>HERMES-WIRE: HERitable MicrobiomE Structure — Workflow for Interpreting host–microbiome Relationships & Effects.</i>
> Centre for Genomic Regulation (CRG), Barcelona.
> [https://github.com/Baud-lab/hermes-wire](https://github.com/Baud-lab/hermes-wire)

---

## 🧩 <b>Acknowledgements</b>

Developed under the <b>[Baud Lab](https://www.crg.eu/en/programmes-groups/baud-lab)</b> at the <b>[Centre for Genomic Regulation (CRG)](https://www.crg.eu)</b> and <b>[Universitat Pompeu Fabra (UPF)](https://www.upf.edu)</b>, Barcelona.
We acknowledge support from the <b>[Bioinformatics Core Facility](https://www.crg.eu/ca/programmes-groups/bioinformatics-unit)</b> and the <b>CRG [Scientific IT team](https://www.crg.eu/en/content/about-us-administration/scientific-information-technologies)</b>. Furthermore, the data used and tested for the development of this tool was generated in collaboration with <b>[NIDA](https://ratgenes.org)</b> and the <b>[Center for Microbiome Innovation](https://cmi.ucsd.edu)</b>, and the project received the support of <b>[La Caixa Foundation](https://lacaixafoundation.org/en/)</b>.

---

## 🧾 <b>License</b>

© 2025 Centre for Genomic Regulation (CRG) and the authors.
Distributed under the <b>[Apache License 2.0](https://github.com/Baud-lab/hermes-wire/blob/main/LICENSE)</b>.

---








