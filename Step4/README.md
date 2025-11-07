<h1 align="center" style="font-weight:400; font-size:28px; line-height:1.3;">
  <b style="font-size:30px;">🧬 HERMES-WIRE</b> — Step 4:<br>
  <b style="font-size:30px;">H</b>eritability-based <b style="font-size:30px;">E</b>xploration of host–microbiome <b style="font-size:30px;">R</b>elationships via <b style="font-size:30px;">M</b>ediation and <b style="font-size:30px;">E</b>ffect-modelling <b style="font-size:30px;">S</b>tructures
</h1>

<p align="center" style="font-weight:400; font-size:16px;">
  <i>Integrative quantitative-genetic analysis of microbial taxa, guilds, and diversity traits using host genome data, 
  cross-fit sentinel mapping, and structural mediation inference.</i>
</p>

---

## 🧩 Overview

**Step 4** of the **HERMES-WIRE** pipeline performs the final *genome-to-microbiome causal dissection* by integrating:
- **Heritable microbiome traits** (taxa, guilds, α/β-diversity)  
- **Candidate host genomic windows**
- **Cross-validated sentinel selection**
- **Mediation and Mendelian Randomization (MR) analyses**
- **Host–microbiome genetic correlation and BLUP estimation**

It bridges host quantitative-genetic signals with microbial community traits to identify **functional genetic mediators**, **crypt-colonisation effects**, and **metabolic outcome pathways**.

---

## ⚙️ Pipeline Summary

| Phase | Description | Output folder |
|:------|:-------------|:--------------|
| **0** | Convert HDF5 → GDS (candidate windows) and build SRM matrices (dam/cage effects). | `phase0` |
| **1–2** | Prepare PLINK data, annotate genes, subset individuals, and create LOCO GRMs. | `phase1`, `loco_grm` |
| **3–5** | Run GENESIS association tests per *species/genus*, perform k-fold crossfit sentinel selection. | `phase3_*`, `phase5_*` |
| **6–7** | Compute out-of-fold (OOF) species betas, meta-analyze across folds, and correlate β with heritability (PGLS). | `phase6_*`, `phase7_pgls` |
| **8** | Integrate PGLS results and synthesize per-genus gene selection for mediation. | `phase8_synthesis_simple` |
| **9** | Perform mediation / MR meta-analysis across folds, producing causal chains and significant IV sets. | `phase9_mediation` |
| **Optional** | Estimate host–microbiome genetic correlations and polygenic scores (RG / BLUP / cvBLUP). | `rg_*`, `cvblup_*` |

---

## 📁 Directory Layout

```

Step4/
├── main.nf                 # Main Nextflow pipeline (DSL2)
├── nextflow.config          # Execution profiles, defaults, and resources
├── submit_nf.sh             # Example SLURM submission wrapper
├── bin/                     # R and bash helper scripts
├── parameters/
│   ├── params_bac.yaml      # Run configuration – Bacteroides
│   ├── params_prev.yaml     # Run configuration – Prevotella
│   ├── params_alpha.yaml    # Run configuration – Alpha diversity (PD q2)
│   ├── params_beta.yaml     # Run configuration – Beta diversity (PCoA1)
│   └── params_guild.yaml    # Run configuration – Guild 3 (co-abundance)
├── dockerfiles/             # Optional custom build contexts
├── singularity_containers/  # Local .sif containers (if not using Docker)
├── input/                   # Input data files (HDF5, metadata, etc.)
└── work/                    # Nextflow work directory (auto-generated)

````

---

## 🚀 Running the Pipeline

### 1. Activate Nextflow environment
```bash
module load nextflow
````

### 2. Submit a job (example)

```bash
bash submit_nf.sh \
  -p parameters/params_bac.yaml \
  -w /scratch/fmorillo/Step4_work \
  -o results_bac
```

### 3. Alternatively, run interactively

```bash
nextflow run main.nf -params-file parameters/params_prev.yaml -resume
```

---

## 🔧 Input Requirements

| Parameter                                            | Description                                              |
| :--------------------------------------------------- | :------------------------------------------------------- |
| `h5_path`                                            | HDF5 genotype matrix (from Step 1)                       |
| `genes_xlsx`                                         | Curated list of candidate host genes (strict intragenic) |
| `metadata_tsv`                                       | Cohort metadata (ID, sex, batch, etc.)                   |
| `grm_rdata`                                          | Host genomic relationship matrix (RData)                 |
| `resid_rda`, `beta_rda`, `alpha_rda`, `clusters_rda` | Residual microbiome trait matrices from Step 2           |
| `herit_rdata`                                        | Species-level heritability estimates                     |
| `taxonomic_file_*`                                   | GTDB taxonomy and phylogenetic tree                      |
| `med_glucose_rdata`                                  | Phenotype residuals (e.g. glucose, BMI)                  |

---

## 🧠 Analysis Modules

* **Crossfit GENESIS association tests:**
  Identify species/genus-level associations within folds.

* **Sentinel detection & ACAT aggregation:**
  Summarize multi-species evidence into *per-gene sentinel SNPs*.

* **PGLS correlation of β vs h²:**
  Test whether species with high heritability have stronger host-gene associations.

* **Mediation / MR:**
  Estimate direct, indirect, and total effects from host genes → microbiome → metabolic traits.

* **Genetic correlation / BLUP (optional):**
  Evaluate polygenic covariance between microbial and metabolic phenotypes.

---

## 🧬 Example Runs

| Trait           | YAML                           | Description                                                |
| :-------------- | :----------------------------- | :--------------------------------------------------------- |
| **Prevotella**  | `parameters/params_prev.yaml`  | Full causal analysis for *Prevotella* species              |
| **Bacteroides** | `parameters/params_bac.yaml`   | Includes genetic correlation + BLUP estimation             |
| **α-diversity** | `parameters/params_alpha.yaml` | Mediation using PD q2 as microbial trait                   |
| **β-diversity** | `parameters/params_beta.yaml`  | Mediation using PCoA1                                      |
| **Guild 3**     | `parameters/params_guild.yaml` | Host–guild interactions (Prevotella + Bacteroides cluster) |

---

## 🧾 Outputs (key files)

| Folder                     | Key Output                                       | Description                                          |
| :------------------------- | :----------------------------------------------- | :--------------------------------------------------- |
| `phase5_sentinels_folds/`  | `sentinels_genus__*__fold*.tsv`                  | Sentinel SNPs per genus/fold                         |
| `phase6_beta_cf/`          | `beta_cf_by_species__from_genus_crossfit__*.tsv` | Cross-fit β values                                   |
| `phase7_pgls/`             | `correlations_by_gene__*.tsv`                    | PGLS β–h² correlations                               |
| `phase8_synthesis_simple/` | `selected_genes.RData`                           | Significant host genes for mediation                 |
| `phase9_mediation/`        | `mediation_meta__*.tsv`, `mediation_sig_iv.tsv`  | Meta-Mediation & IV results                          |
| `rg_* / cvblup_*`          | GREML / BLUP summaries                           | (optional) genetic correlation and effect prediction |

---

## 🧠 Citation and Credits

If you use this workflow, please cite:

> **Morillo F. M. S. D.**, *et al.* (in prep).
> *Host genetics shapes crypt niche colonisation by keystone gut bacteria, influencing metabolic health.*
> Centre for Genomic Regulation (CRG) & Universitat Pompeu Fabra (UPF), Barcelona.

For code and containers:

* **GitHub:** [Baud-lab/hermes-wire](https://github.com/Baud-lab/hermes-wire)
* **DockerHub:** `morillofe/microbiome_geno_bioconductor`

---

Would you like me to extend this README with the **SLURM `submit_nf.sh` usage block and parameter explanations** (like those in Step 2, e.g., `--trait`, `--genera`, `--outdir`, etc.)?
I can also format it for **Zenodo or GitHub rendering** (with collapsible `<details>` sections).
