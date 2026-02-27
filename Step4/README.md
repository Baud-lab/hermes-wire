# 🧬 HERMES-WIRE — Step 4: Integrative quantitative-genetic analysis of microbial taxa, guilds, and diversity traits using host genome data, cross-fit sentinel mapping, and structural mediation inference

<b>Centre for Genomic Regulation (CRG), Barcelona — 2025</b>  
<b>Authors:</b> Dr Felipe Morillo Sanz Dias, Dr Luca Cozzuto, Dr Hélène Tonnelé, Dr Amelie Baud 

---

## 📘 Overview

**Step 4** of the **HERMES-WIRE** pipeline performs the final *genome-to-microbiome causal dissection* by integrating:
- **Heritable microbiome traits** (taxa, guilds, α/β-diversity)  
- **Candidate host genomic windows**
- **Cross-validated sentinel selection**
- **Mediation and Mendelian Randomization (MR) analyses**

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

---

## 🧬 Example Runs

| Trait           | YAML                           | Description                                                |
| :-------------- | :----------------------------- | :--------------------------------------------------------- |
| **Prevotella**  | `parameters/params_prev.yaml`  | Causal analysis for the *Prevotella* genus                 |
| **Bacteroides** | `parameters/params_bac.yaml`   | Causal analysis for the *Bacteroides* genus                |
| **α-diversity** | `parameters/params_alpha.yaml` | Causal analysis for phylogenetic diversity                 |
| **β-diversity** | `parameters/params_beta.yaml`  | Causal analysis for beta diversity PCoA1                   |
| **Guild 3**     | `parameters/params_guild.yaml` | Causal analysis for a bacterial species guild              |

---

## 🧾 Outputs (key files)

| Folder                     | Key Output                                       | Description                                          |
| :------------------------- | :----------------------------------------------- | :--------------------------------------------------- |
| `phase5_sentinels_folds/`  | `sentinels_genus__*__fold*.tsv`                  | Sentinel SNPs per genus/fold                         |
| `phase6_beta_cf/`          | `beta_cf_by_species__from_genus_crossfit__*.tsv` | Cross-fit β values                                   |
| `phase7_pgls/`             | `correlations_by_gene__*.tsv`                    | PGLS β–h² correlations                               |
| `phase8_synthesis_simple/` | `selected_genes.RData`                           | Significant host genes for mediation                 |
| `phase9_mediation/`        | `mediation_meta__*.tsv`, `mediation_sig_iv.tsv`  | Meta-Mediation & IV results                          |

---

## ✅ **Validation:**

The results reported on the thesis: <i>Host genetics shapes mucin niche colonisation by keystone gut bacteria, influencing metabolic health</i> (Dias, FMS; 2026) were already obtained by using the full Nextflow pipeline. Please, let us know if you found any problems trying to run it or if you miss something that you cannot find here on in the extra-codes that you consider necessary to replicate the results of the thesis.

<b>Citation:</b> Morillo Sanz Dias, F. (2026). Host genetics shapes crypt niche colonisation by keystone gut bacteria, influencing metabolic health. Zenodo. https://doi.org/10.5281/zenodo.18663497

---

## 🧩 <b>Acknowledgements</b>

Developed under the <b>[Baud Lab](https://www.crg.eu/en/programmes-groups/baud-lab)</b> at the <b>[Centre for Genomic Regulation (CRG)](https://www.crg.eu)</b> and <b>[Universitat Pompeu Fabra (UPF)](https://www.upf.edu)</b>, Barcelona.
We acknowledge support from the <b>[Bioinformatics Core Facility](https://www.crg.eu/ca/programmes-groups/bioinformatics-unit)</b> and the <b>CRG [Scientific IT team](https://www.crg.eu/en/content/about-us-administration/scientific-information-technologies)</b>. Furthermore, the data used and tested for the development of this tool was generated in collaboration with <b>[NIDA](https://ratgenes.org)</b> and the <b>[Center for Microbiome Innovation](https://cmi.ucsd.edu)</b>, being available on <b>[Qiita - Study 11479](https://qiita.ucsd.edu/study/description/11479)</b>. The project received the support of <b>[La Caixa Foundation](https://lacaixafoundation.org/en/)</b>.

---

## 🧠 <b>Citation</b>

> Morillo FMSD, Cozzuto L, Tonnelé H, Baud A (2026). <i>HERMES-WIRE: HERitable MicrobiomE Structure — Workflow for Interpreting host–microbiome Relationships & Effects.</i>
> Centre for Genomic Regulation (CRG), Barcelona.
> [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18801938.svg)](https://doi.org/10.5281/zenodo.18801938)

---

## 🧾 <b>License</b>

© 2025 Centre for Genomic Regulation (CRG) and the authors.
Distributed under the <b>[Apache License 2.0](https://github.com/Baud-lab/hermes-wire/blob/main/LICENSE)</b>.







