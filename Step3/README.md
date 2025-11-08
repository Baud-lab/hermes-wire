Great — thanks for clarifying. Here is the **README.md for Step 3** of the HERMES‑WIRE pipeline, matching the style of Step 4 and incorporating your directives.

---

# 🧬 HERMES-WIRE — Step 3: Genome functional annotation & gene-family profiling (COG/KEGG/CAZy)

**Centre for Genomic Regulation (CRG), Barcelona — 2025**
**Authors:** Felipe Morillo Sanz Dias, Luca Cozzuto, Hélène Tonnele, Amelie Baud

---

## 📘 Overview

**Step 3** of the HERMES-WIRE pipeline performs functional annotation of microbial genomes and builds species-level gene-family matrices. Specifically:

* We begin from protein FASTA files (in our case, the protein FASTAs from the GTDB R207 reference set).
* We annotate each protein with functional categories (COGs/ENOGs, KEGG orthology/EC) via eggNOG‑mapper.
* We annotate carbohydrate‐active enzyme (CAZyme) families via dbCAN2 (HMM + DIAMOND approach). ([OUP Academic][1])
* We summarise results to produce matrices of *function vs species* (for the set of species represented in the input FASTAs).
* These matrices serve as input for subsequent quantitative and comparative analyses (used in Step 4 and downstream).

---

## ⚙️ Pipeline Summary

| Phase | Description                                                                           | Output folder                |
| ----- | ------------------------------------------------------------------------------------- | ---------------------------- |
| **0** | Gather input protein FASTA files (GTDB R207) and prepare directory structure.         | `phase0`                     |
| **1** | Run eggNOG-mapper functional annotation on the protein FASTAs; generate .tsv outputs. | `phase1_eggnog`              |
| **2** | Run CAZyme annotation (dbCAN2) on the protein FASTAs; parse and summarise hits.       | `phase2_cazyme`              |
| **3** | Consolidate annotation tables: link species ↔ COG/KEGG/ENOG/CAZy hits.                | `phase3_species_annotations` |
| **4** | Build gene‐family abundance matrices (species × function) for COG/KEGG and CAZy.      | `phase4_matrices`            |
| **5** | Quality control & filtering (e.g., remove low‐coverage species, normalise counts).    | `phase5_filtered_matrices`   |

---

## 📁 Directory Layout

```
Step3/
├── main.nf                         # Main Nextflow pipeline (DSL2)
├── nextflow.config                 # Execution profiles, defaults, resources
├── submit_nf.sh                    # Example SLURM submission wrapper
├── bin/                            # R and bash helper scripts (parsing, summarising)
├── parameters/
│   └── params_step3.yaml           # Configure run (input paths, thresholds etc.)
├── singularity_containers/         # (Optional) Local .sif containers if not using Docker
├── input/
│   ├── proteins/                    # Input protein FASTA files (GTDB R207)
│   └── metadata/                    # Species metadata (IDs, names, GTDB taxonomy)
└── work/                           # Nextflow work directory (auto-generated)
```

---

## 🚀 Running the Pipeline

### 1. Load Nextflow environment

```bash
module load nextflow
```

### 2. Submit a job (example)

```bash
bash submit_nf.sh \
  -p parameters/params_step3.yaml \
  -w /scratch/fmorillo/Step3_work \
  -o results_step3
```

### 3. Alternatively, run interactively

```bash
nextflow run main.nf -params-file parameters/params_step3.yaml -resume
```

---

## 🔧 Input Requirements

| Parameter                  | Description                                                             |
| -------------------------- | ----------------------------------------------------------------------- |
| `protein_fastas_dir`       | Directory containing all species’ protein FASTA files (e.g., GTDB R207) |
| `species_metadata_tsv`     | Metadata table: species ID, GTDB taxonomy, file path etc.               |
| `eggnog_db`                | Path to eggNOG-mapper database reference                                |
| `cazyme_db`                | Path to CAZyme (dbCAN2) database files                                  |
| `min_proteins_per_species` | Minimum number of proteins to retain a species                          |
| `function_output_prefix`   | Prefix for functional annotation tables                                 |

---

## 🧠 Analysis Modules

* **eggNOG annotation:** Uses eggNOG-mapper on the input protein FASTAs to assign COGs/ENOGs, KEGG orthologs, and EC numbers; summarises per‐species counts.
* **CAZyme annotation:** Utilises the dbCAN2 pipeline (HMMER vs dbCAN HMMdb + DIAMOND vs CAZy database) for carbohydrate‐active enzyme family identification. ([OUP Academic][1])
* **Matrix construction:** Aggregates per‐species annotation hits into abundance matrices (species × function families), suitable for comparative analysis (taxon‐function associations, diversity metrics, guild profiling).
* **Quality filtering & normalisation:** Optionally removes species with insufficient proteins or annotation coverage, and normalises counts (e.g., per-million, log-transform) to ensure downstream compatibility.

---

## 🧬 Example Runs

| Run Name           | YAML file                     | Description                                            |
| ------------------ | ----------------------------- | ------------------------------------------------------ |
| **Bacteroides**    | `parameters/params_bac.yaml`  | Full annotation matrix for *Bacteroides* species group |
| **Prevotella**     | `parameters/params_prev.yaml` | Annotation for *Prevotella* species group              |
| **Paraprevotella** | `parameters/params_para.yaml` | Annotation for *Paraprevotella* group                  |

---

## 🧾 Outputs (key files)

| Folder                      | Key Output File                          | Description                                   |
| --------------------------- | ---------------------------------------- | --------------------------------------------- |
| `phase1_eggnog/`            | `species_annotation_eggnog.tsv`          | Species × eggNOG annotation summary           |
| `phase2_cazyme/`            | `species_annotation_cazy.tsv`            | Species × CAZy family summary                 |
| `phase4_matrices/`          | `species_x_function_matrix.tsv`          | Combined annotation matrix (COG/KEGG or CAZy) |
| `phase5_filtered_matrices/` | `filtered_matrix_species_x_function.tsv` | QC-filtered and normalised matrix             |

---

## 🧱 Optional: Fetching Genomes from NCBI and Generating Protein FASTAs

Although Step 3 typically uses the GTDB R207 protein FASTAs, the following workflow is provided for users who wish to retrieve genomes directly from NCBI (for example to extend beyond GTDB or include custom accessions) and convert them into protein FASTAs suitable for annotation.

### Step A – Select accession list

1. Choose your target taxa (e.g., *Bacteroides*, *Prevotella*, *Paraprevotella*).
2. Create a list of all genomes (clusters) you want to download (accession numbers).

### Step B – Download genome FASTA files

```bash
# To retrieve GCA_ accessions (from GenBank)
singularity exec ncbi-genome-download.sif \
  ncbi-genome-download --assembly-accessions genomes_v3.txt \
    -s genbank -F fasta -o /path/to/output/genomes_comparison_v3 bacteria

# To retrieve GCF_ accessions (from RefSeq)
singularity exec ncbi-genome-download.sif \
  ncbi-genome-download --assembly-accessions genomes_v2.txt \
    -s refseq -F fasta -o /path/to/output/genomes_comparison_v2 bacteria
```

> **Note:** GCF_ accession numbers come from RefSeq assemblies; GCA_ accessions come from GenBank. Both databases may contain bacterial genomes. In the command you must specify `bacteria` as the taxonomic group.

### Step C – Move extracted files

```bash
find /path/to/output/genomes_comparison_v2 -mindepth 2 -type f -name "*.gz" \
     -exec mv -i {} /path/to/output/genomes_comparison_v2 \;
```

### Step D – Predict proteins with Prodigal

```bash
INPUT_P=$(ls /path/to/output/genomes_comparison_v3/*_genomic.fna | sort -u | sed -n ${SGE_TASK_ID}p)
O=$(echo $INPUT_P | sed -r 's/.{69}//' | sed 's/_genomic.fna//')
OUT="/path/to/output/genomes_comparison_v3"

echo "Started"
prodigal -i ${INPUT_P} -a ${OUT}/proteins/${O}_protein.fasta -t $NSLOTS
echo "Prodigal ended"
```

Once you have the protein FASTA files (`*_protein.fasta`), you can feed them into Step 3 of the pipeline (instead of the GTDB R207 set).

---

## 🧩 Acknowledgements

Developed under the [Baud Lab](https://www.crg.eu/en/programmes-groups/baud-lab) at the [Centre for Genomic Regulation (CRG)](https://www.crg.eu) and [Universitat Pompeu Fabra (UPF)](https://www.upf.edu), Barcelona. We thank the CRG Bioinformatics Core and Scientific IT Team for compute infrastructure. The GTDB data, eggNOG, and dbCAN resources were essential for this workflow.

---

## 🧠 Citation

> Morillo FMSD, Cozzuto L, Tonnele H, Baud A (2025). *HERMES-WIRE: HERitable MicrobiomE Structure — Workflow for Interpreting host–microbiome Relationships & Effects.* Centre for Genomic Regulation (CRG), Barcelona.
> [https://github.com/Baud-lab/hermes-wire](https://github.com/Baud-lab/hermes-wire)

---

## 🧾 License

© 2025 Centre for Genomic Regulation (CRG) and the authors. Distributed under the [Apache License 2.0](https://github.com/Baud-lab/hermes-wire/blob/main/LICENSE).

---

If you like, I can prepare a **markdown version with badges** (CI status, license, version) and upload ready for your repo. Would you like me to do that?

[1]: https://academic.oup.com/nar/article/46/W1/W95/4996582?utm_source=chatgpt.com "dbCAN2: a meta server for automated carbohydrate-active enzyme ..."
