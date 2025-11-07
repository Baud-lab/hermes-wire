<h1 align="left" style="font-weight:400; font-size:28px; line-height:1.3;">
  <b style="font-size:30px;">HERMES-WIRE</b>
  (<b style="font-size:30px;">HER</b>itable
  <b style="font-size:30px;">M</b>icrobiom<b style="font-size:30px;">E</b>
  <b style="font-size:30px;">S</b>tructure —
  <b style="font-size:30px;">W</b>orkflow for
  <b style="font-size:30px;">I</b>nterpreting host–microbiome
  <b style="font-size:30px;">R</b>elationships &amp;
  <b style="font-size:30px;">E</b>ffects)
</h1>

<p align="center">
  <img src="assets/hermes-wire-logo.png" alt="HERMES-WIRE logo: Hermes with DNA and gut bacteria" width="320">
</p>

A modular suite of Nextflow pipelines that fuse cost-effective quantitative genetics with systems biology to map how host polygenic variation wires gut microbial community structure and model their causal impacts on metabolic health.

<b>Centre for Genomic Regulation (CRG), Barcelona — 2025</b>  
<b>Authors:</b> Felipe Morillo Sanz Dias, Luca Cozzuto, Hélène Tonnele, Amelie Baud  

---

## <b>STEPS:</b>

### <b>Step 1: Mapping bacterial taxa and/or functions</b>
#### Objective: To maximise microbiome profiling sensitivity from shallow-shotgun data, while using broad reference catalogues.
#### Main stages:
1. Data processing: Trimming and removal of host reads;
2. Taxonomic profiling: Alignment with a reference catalogue and read counting (Species-Phylum);
3. Functional profiling: Alignment with a reference catalogue and read counting (EC4-EC1);

### <b>Step 2: Microbiome characterisation</b>
#### Objectives: To minimise the risks of mapping false positives and characterise the microbiome at different levels.
#### Main stages:
1. Matrix processing: Sample and taxa filtering, aggregation to higher taxonomic/functional levels, data transformation, rank-normalisation and residualisation (removal of fixed-effects)
2. Analyses: Co-abundance guilds, alpha- and beta-diversity, enterotypes and networks
3. Heritability tests

### <b>Step 3: Bacterial functional enrichment analysis linked to heritability</b>
#### Objectives: To identify enriched functions (COGs/ENOGs) in the genomes of species with high- and low-heritability values within the same genus.
#### Main stages:
1. Functional annotation of genomes: Already using protein catalogues per genome;
2. Phylogenetic-aware enrichment analysis;

### <b> Step 4: A hypothesis-driven selection of host candidate genes for causal inference analyses</b>
#### Objectives: To delimitate, based on the results obtained in the previous steps, what host genes might be associated with polygenic effects in bacterial groups with high heritability and use such genes for causal inference analyses.
#### Main stages:
1. Filtering SNPs based on the list of candidate genes in the genotyping panel;
2. Select sentinel SNPs per gene on the genus level (Cross-validation: training set);
3.  Test whether the additive effects of such sentinel SNPs across species from the same genus were associated with the species heritability values (Cross-validation: testing set);
4.  Run mediation analysis and Mendelian randomisation (MR) crossing host genes, microbiome features and host phenotypes (Cross-validation: testing set);
5. Define a structured equation model (SEM) for the multi-level chain of events connecting keystone species with host health-associated traits after MR validation.

## <b>PREREQUISITES</b>
In order to run the pipelines, you will need Nextflow and Singularity for containerisation. 

### Installing Nextflow:

```
wget -qO- https://get.nextflow.io | bash
```

### Singularity: [here](https://docs.sylabs.io/guides/3.1/user-guide/quick_start.html#quick-installation-steps)


In CRG's cluster you can log to the NEXTFLOW node:

```
ssh -Y YOURUSER@nextflow.hpc.crg.es
```

Then, you just need to modify once your bashrc file for loading both Java and singularity. You need also to specify some folders for storing the cache for singularity. So you might want to create those folders first.

```
mkdir $HOME/tmp
mkdir $HOME/singularity_containers/
mkdir $HOME/tmp_singu
```
Finally,

```
vi $HOME/.bashrc

# ADD THIS FOR LOADING JAVA AND SINGULARITY
module use /software/as/el7.2/EasyBuild/CRG/modules/all
module load Java/11.0.2
module load Singularity/3.7.0
# Add this for using more space for temporary images
export SINGULARITY_CACHEDIR=$HOME/tmp
export NXF_SINGULARITY_CACHEDIR="$HOME/singularity_containers/"
export SINGULARITY_TMPDIR=$HOME/tmp_singu

```

Log-out and log-in to reset:

```
exit 

# The system will print this: 
logout
Connection to nextflow.hpc.crg.es closed.

# Log in again
ssh -Y YOURUSER@nextflow.hpc.crg.es
```

Then you can install NextFlow as indicated previously.

## <b>INSTALLING THE SUITE</b>

```
git clone --recurse-submodules git@github.com:Baud-lab/hermes-wire.git
```

## <b>NEXTFLOW CONFIG FILE</b>

The file <b>`nextflow.config`</b> indicates how the pipeline’s processes will run in parallel according to their channels, in accordance with the cluster specifications being used, such as the different queues available for job scheduling. At CRG, the workload manager and job scheduler used for the institution’s HPC cluster is SLURM (Simple Linux Utility for Resource Management), and for this reason, the environments/profiles defined in this file are currently set to this system, but this can be changed in case of other HCP managing systems being used.

## <b>RUNNING A PIPELINE (Example: Step1)</b>

```
sbatch submit_nf.sh main.nf -profile singularity,slurm_genoa -params-file params.yaml -resume -w work
```

This command line can be seen in the <b>`cmd`</b> files present in all step directories. The <b>`main.nf`</b> file, also known as the pipeline “manifest”, is the primary code that defines how processes will be organised in different workflows and how initial parameters will determine the various directions the pipeline will follow. Each <b>`main.nf`</b> file is supported by a <b>`local_modules.nf`</b>, loading the different processes and workflows these modules comprehend, along with the specific requirements in terms of container, running time, memory allocation, inputs and scripts used. The latter are packed in the folder <b>`bin`</b> and can be made of many programming languages, such as R, Python or bash.

To stop a run, you need to kill the Nextflow main process. The **PID** is indicated by the hidden file <b>`.nextflow.pid`</b> 

```
cat .nextflow.pid | xargs kill
```

the <b>`params.yaml`</b> file contains the information for running the pipeline, such as input files and arguments to be incorporated by the scripts and the corresponding tools/libraries/packages of each process. Here is an example, considering Step1:

```
### Parameters of the HERMES-WIRE Step 1 Pipeline

### STEP 1: Data preprocessing and profiling

### Modules
module_preprocessing:         "YES"
module_taxonomic_mapping:     "YES"
module_functional_mapping:    "YES"

### Tool options
tool_opt: "./tool_opt.tsv"

### Pre-processing
preprocessing_file_raw_reads: "./dataset_raw/demo/*_R{1,2}_001.fastq.gz"
preprocessing_file_adapters: "./input/Adapters_Qiita.fa"
preprocessing_file_host_genome: "./input/mRatBN7.2.fna" # Dowload https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/227/675/GCF_015227675.2_mRatBN7.2/GCF_015227675.2_mRatBN7.2_rna.fna.gz and unzip

### Non-host reads (If you don't select Module Pre-processing)
file_non_host_reads: "./dataset_non_host/*_{1,2}.fq.gz"  # Keep empty if no files are available

### Taxonomic
taxonomic_file_microbiome_index: "./input/indices/kaiju/gtdb_index.fmi" # Download from https://zenodo.org/uploads/17545483 and unzip
taxonomic_file_genome_sizes: "./input/indices/kaiju/catalogue_stats.csv"
taxonomic_parameter_normalization_by_genome_size: 'YES'

### Functional
functional_file_ref_catalogue: ""
functional_file_index: "./input/indices/prromenade/bactvirus2020*" # Download from https://zenodo.org/uploads/17545032 and unzip
functional_file_ec: "./input/indices/prromenade/taxid_name.txt"

### Output folder
output_folder: "./Results_Step1/"
```

## <b>USING THE NEXTFLOW TOWER</b>

Log in to the Tower website: https://tower.nf/. Go on `Sign-in`. Then authenticate with either **GitHub** or **Google**.  
Then go on the icon on the top right with your face or a symbol and **Your tokens**. Add a token, choosing a name, and copy the value. 

Add to your `$HOME/.bashrc` as previosly indicated. 

```
export TOWER_ACCESS_TOKEN=<YOUR ACCESS TOKEN>
```

Log-off and Log-in.

Then launch the pipeline using `-with-tower` as an argument for the command line above. You just need to set the token once, and then you can use Tower every time you run a pipeline.

You can then follow the results going to link indicated when executing the pipeline in your log file (in the case of SLURM users the submit_nf.sh file will generate files with names starting with “slurm” on the corresponding step directory).

---

## <b>🧠 Citation</b>

If you use this workflow or its pre-computed indices, please cite:

> Morillo FMSD, Cozzuto L, Tonnele H, Baud A (2025). <i>HERMES-WIRE: HERitable MicrobiomE Structure — Workflow for Interpreting host–microbiome Relationships & Effects.</i>
> Centre for Genomic Regulation (CRG), Barcelona.
> [https://github.com/Baud-lab/hermes-wire](https://github.com/Baud-lab/hermes-wire)

---

## <b>🧩 Acknowledgements</b>

Developed under the <b>Baud Lab</b> at the <b>Centre for Genomic Regulation (CRG)</b> and <b>Universitat Pompeu Fabra (UPF)</b>.
We acknowledge the support of the <b>Bioinformatics Core Facility at CRG</b>.

---

## <b>🧾 License</b>

© 2025 Centre for Genomic Regulation (CRG) and the authors.
Distributed under the <b>MIT License</b>.

---
