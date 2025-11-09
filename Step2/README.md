# 🧬 HERMES-WIRE — Step 2: Profile Matrix Processing and Analyses

<b>Centre for Genomic Regulation (CRG), Barcelona — 2025</b>  
<b>Authors:</b> Felipe Morillo Sanz Dias, Luca Cozzuto, Hélène Tonnele, Amelie Baud  

---

## <b>📘 Overview</b>

<b>HERMES-WIRE Step 2</b> is the second stage of the <b>HERMES-WIRE</b> suite (<b>HER</b>itable <b>M</b>icrobiom<b>E</b> <b>S</b>tructure — <b>W</b>orkflow for <b>I</b>nterpreting host–microbiome <b>R</b>elationships &amp; <b>E</b>ffects).

It integrates processed microbiome profiles (taxonomic or functional) and performs downstream analyses, including:

1. **Matrix processing** — Filtering, depth normalization, transformation, and residualization.  
2. **Community structure** — Clustering, network construction, and enterotype detection.  
3. **Diversity metrics** — Alpha and beta diversity estimation.  
4. **Host–microbiome interaction modelling** — Heritability estimation and GWAS.  
5. **Cross-study harmonization** — Comparing Shallow vs 16S metagenomes when applicable.

---

## <b>⚙️ Directory structure</b>

```

Step2/
├── bin/                     # Executable scripts (R, Python, Bash)
├── cmd                      # Command template
├── dockerfiles/             # Singularity/Docker recipes
├── input/                   # Metadata, taxonomy files, and input matrices
├── local_modules.nf         # Step-specific modules
├── main.nf                  # Main Nextflow pipeline
├── nextflow.config          # Runtime profiles (SLURM, Singularity, etc.)
├── parameters/              # Example parameter YAMLs (Shallow, Functional, Harmonized, etc.)
├── singularity_containers/  # Local Singularity image cache
├── submit_nf.sh             # Example SLURM submission script
└── work/                    # Nextflow working directory (auto-generated)

```

---

## <b>🚀 Quick Start</b>

### <b>1 – Dependencies</b>

Requirements:

- **Nextflow ≥ 23.10**
- **Singularity / Apptainer**
- **R ≥ 4.2**, Python ≥ 3.8 (via containers)
- Optional: HPC cluster with SLURM

```bash
module load Nextflow
module load Singularity
````

---

### <b>2. Download required HDF5 file (for analyses with the selected HS rats)</b>

The pipeline depends on a big HDF5 file hosted on <b>Zenodo</b>:

| File type | Contents                                                                     | Zenodo link                                                                |
| -------------- | ------------------------------------------------------------------------------- | -------------------------------------------------------------------------- |
| **HDF5**  | GRM and Gwnotypes (`P50_rats_Rn7.h5`)                        | [[https://zenodo.org/uploads/17545483](https://zenodo.org/records/17546688)] |

Download and extract the archive into `input/` using commands ```wget``` and ```tar -xvzf```, respectivelly.

---

### <b>3 – Select a configuration</b>

Templates are available under the <code>parameters/</code> folder:

| Configuration                                  | Purpose                                               | YAML file                                      |
| ---------------------------------------------- | ----------------------------------------------------- | ---------------------------------------------- |
| **Shallow Full Taxonomic**                     | Standard run using taxonomic profiles from Step 1     | `parameters/full_shallow/params_shallow_tax.yaml`     |
| **Shallow Full Functional**                    | Standard run using functional profiles (EC tables)    | `parameters/full_shallow/params_shallow_func.yaml`    |
| **Shallow ↔ 16S Harmonization (full)**         | Comparisons between Shallow and 16S with samples and taxa harmonised | `parameters/comparative/params_<shallow or 16S>_full_harm.yaml`    |
| **Shallow ↔ 16S Harmonization (samples only)** | Comparisons between Shallow and 16S with only samples harmonised       | `parameters/comparative/params_<shallow or 16S>_sample_harm.yaml` |

Each YAML file defines all modules (Matrix Processing, Clusters, Networks, Enterotypes, Heritability, GWAS) and their parameters.

---

### <b>4 – Run the pipeline</b>

```bash
nextflow run main.nf -params-file parameters/full_shallow/params_shallow_tax.yaml -profile singularity
```

or via SLURM:

```bash
sbatch submit_nf.sh
```

You can monitor execution on **Nextflow Tower** by adding `-with-tower`.

---

## <b>🧩 Modules and workflows</b>

| Module                          | Description                                                            | Key tools/processes        |
| :------------------------------ | :--------------------------------------------------------------------- | :------------------------- |
| **Comparisons (Harmonization)** | Align and harmonize taxa across datasets (e.g. 16S vs Shallow)         | R (custom scripts)         |
| **Matrix Processing**           | Filtering, prevalence cutoff, CLR transformation, covariate regression | R, HDF5 output             |
| **Clusters (UMAP + HDBSCAN)**   | Identify guilds / clusters based on co-abundance structure             | umap-learn, hdbscan, R     |
| **Enterotypes**                 | Define enterotypes (k-means clustering on residualized data)           | R (kmeans)                 |
| **Networks**                    | Build and analyze co-abundance networks (modules, hubs)                | igraph, WGCNA              |
| **Alpha/Beta Diversity**        | Shannon, Simpson, Hill numbers and PCoA based metrics                  | vegan, phyloseq            |
| **Heritability**                | Host genetic variance partition using GCTA-style LMM                   | GENESIS LMM, custom R      |
| **GWAS (optional)**             | Genome-wide mapping of microbial traits                                | GPKronSum, LMM via GENESIS |

---

## <b>📂 Outputs</b>

Results are written to the directory defined in <code>output_folder</code>:

```
Results_*/ 
├── Matrix_Processing/         # Filtered, transformed, residualized matrices (H5/RData)
├── Diversity/                 # Alpha and Beta diversity results
├── Clusters/                  # UMAP coordinates, cluster membership, DBCV/ARI stats
├── Networks/                  # Co-abundance network files and hub metrics
├── Enterotypes/               # Enterotype assignments and stability metrics
├── Heritability/              # h2 estimates, FDR-significant traits
├── GWAS/                      # SNP–trait association results (if enabled)
└── pipeline_info/             # Nextflow logs, trace, and reports
```

---

## <b>🧠 Citation</b>

If you use this module or its outputs, please cite:

> **Morillo FMSD, Cozzuto L, Tonnele H, Baud A (2025)**
> *HERMES-WIRE: HERitable MicrobiomE Structure — Workflow for Interpreting host–microbiome Relationships & Effects*
> Centre for Genomic Regulation (CRG), Barcelona.
> [https://github.com/Baud-lab/hermes-wire](https://github.com/Baud-lab/hermes-wire)

---

## 🧩 <b>Acknowledgements</b>

Developed under the <b>[Baud Lab](https://www.crg.eu/en/programmes-groups/baud-lab)</b> at the <b>[Centre for Genomic Regulation (CRG)](https://www.crg.eu)</b> and <b>[Universitat Pompeu Fabra (UPF)](https://www.upf.edu)</b>, Barcelona.
We acknowledge support from the <b>[Bioinformatics Core Facility](https://www.crg.eu/ca/programmes-groups/bioinformatics-unit)</b> and the <b>CRG [Scientific IT team](https://www.crg.eu/en/content/about-us-administration/scientific-information-technologies)</b>. Furthermore, the data used and tested for the development of this tool was generated in collaboration with <b>[NIDA](https://ratgenes.org)</b> and the <b>[Center for Microbiome Innovation](https://cmi.ucsd.edu)</b>. The project received the support of <b>[La Caixa Foundation](https://lacaixafoundation.org/en/)</b>.

---

## 🧾 <b>License</b>

© 2025 Centre for Genomic Regulation (CRG) and the authors.
Distributed under the <b>[Apache License 2.0](https://github.com/Baud-lab/hermes-wire/blob/main/LICENSE)</b>.

---












