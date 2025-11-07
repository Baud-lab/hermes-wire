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

### <b>2 – Select a configuration</b>

Templates are available under the <code>parameters/</code> folder:

| Configuration                                  | Purpose                                               | YAML file                                      |
| ---------------------------------------------- | ----------------------------------------------------- | ---------------------------------------------- |
| **Shallow Full Taxonomic**                     | Standard run using taxonomic profiles from Step 1     | `parameters/full_shallow/params_shallow_taxonomic.yaml`     |
| **Shallow Full Functional**                    | Standard run using functional profiles (EC tables)    | `parameters/full_shallow/params_shallow_functional.yaml`    |
| **Shallow ↔ 16S Harmonization (full)**         | Joint taxonomic harmonization between Shallow and 16S | `parameters/comparative/params_<shallow or 16S>_full_harm.yaml`    |
| **Shallow ↔ 16S Harmonization (samples only)** | Only harmonizes shared samples between datasets       | `parameters/comparative/params_<shallow or 16S>_sample_harm.yaml` |

Each YAML file defines all modules (Matrix Processing, Clusters, Networks, Enterotypes, Heritability, GWAS) and their parameters.

---

### <b>3 – Run the pipeline</b>

```bash
nextflow run main.nf -params-file parameters/params_shallow_taxonomic.yaml -profile singularity
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

## <b>📊 Workflow logic (Mermaid)</b>

```mermaid
flowchart TD
    A([Start]) --> B{Comparisons?}
    B -- YES --> B1[Samples & taxa harmonization<br>(Shallow ↔ 16S)]
    B -- NO --> C{Matrix Processing?}
    B1 --> C

    C -- YES --> D[Filtering & Transformation]
    C -- NO --> Z[Skip to existing H5 / residuals]

    D --> E{Alpha / Beta Diversity?}
    E -- YES --> F[Compute α, β diversity]
    E -- NO --> G

    F --> G{Clusters / Enterotypes / Networks?}
    G -- YES --> H[UMAP + HDBSCAN → Guilds / Networks]
    G -- NO --> I

    H --> I{Heritability or GWAS?}
    I -- YES --> J[H5 creation → LMM heritability + GWAS]
    I -- NO --> K([End])

    Z --> I
    J --> K([End])
```

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

## <b>🧩 Acknowledgements</b>

Developed under the <b>Baud Lab</b> (CRG & UPF).
We acknowledge the support of the <b>CRG HPC Core Facility</b> and the Genetics of Complex Traits group.

---

## <b>🧾 License</b>

© 2025 Centre for Genomic Regulation (CRG) and the authors.
Distributed under the <b>MIT License</b>.

---

