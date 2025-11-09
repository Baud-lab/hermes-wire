<h1 style="font-weight:400; font-size:28px; line-height:1.3;">
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

<i>Named after Hermes, the Greek messenger of the gods, HERMES-WIRE tracks how information flows from multiple host genes to gut microbiome assembly, facilitating the interpretation of how they affect each other. It consists of a modular Nextflow suite combining quantitative genetics with systems biology approaches, revealing how host polygenic variation functionally wires gut microbial communities with impacts in metabolic health.</i>

<b>Centre for Genomic Regulation (CRG), Barcelona — 2025</b><br>
<b>Authors:</b> Felipe Morillo Sanz Dias, Dr Luca Cozzuto, Dr Hélène Tonnelé, Dr Amelie Baud  

---

## 📘 <b>Pipeline Overview</b>

<b>HERMES-WIRE</b> consists of **four modular steps**, each of which can be run independently or sequentially.  
Together they provide a full *host–microbiome systems genetics* framework — from raw reads to host–microbiome causal inference.

This work forms part of the PhD thesis <i>Host genetics shapes mucin niche colonisation by keystone gut bacteria, influencing metabolic health</i> by <b>Felipe Morillo Sanz Dias</b> (submitted November 2025), which was conducted at the Centre for Genomic Regulation (CRG), Department of Systems and Synthetic Biology / Universitat Pompeu Fabra (UPF), Department of Medicine and Life Sciences, with inputs of Dr Amelie Baud. The thesis is currently under evaluation and yet to be published.
Users are kindly invited to report any issues encountered when running the pipeline and/or report missing parts for the total reproduction of the results reported on the thesis.

---

### 🔍 <b>Step 1 — Mapping bacterial taxa and/or functions</b>
#### 🎯 Objective
Maximise microbiome profiling sensitivity from shallow-shotgun data using comprehensive reference catalogues.  
#### 🔬 Main stages
1. **Data pre-processing** — quality trimming and host-read removal.  
2. **Taxonomic profiling** — alignment to a GTDB-based protein catalogue (Species → Phylum).  
3. **Functional profiling** — alignment to IFGP/PRROMenade catalogue (EC4 → EC1 levels).  

---

### 🔍 <b>Step 2 — Microbiome characterisation</b>
#### 🎯 Objective
Minimise mapping artefacts and characterise microbiome structure at multiple levels.  
#### 🔬 Main stages
1. **Matrix processing** — sample/taxon filtering, aggregation, CLR transformation, residualisation.  
2. **Community analyses** — co-abundance guilds, α/β-diversity, enterotypes, and network inference.  
3. **Heritability analysis** — estimate host genetic effects on microbiome features.  

---

### 🔍 <b>Step 3 — Functional enrichment linked to heritability</b>
#### 🎯 Objective
Identify functional gene enrichments (COGs/ENOGs) distinguishing high- vs low-heritability species within the same genus.  
#### 🔬 Main stages
1. **Functional annotation** of genomes using protein catalogues.  
2. **Phylogenetic-aware enrichment** analyses controlling for lineage structure.  

---

### 🔍 <b>Step 4 — Host gene prioritisation and causal inference</b>
#### 🎯 Objective
Select candidate host genes associated with bacterial heritability and model causal host → microbiome → phenotype links.  
#### 🔬 Main stages
1. **Candidate SNP selection** — filter genotyping panel by gene list.  
2. **Sentinel SNP assignment** per gene/genus via cross-validation.  
3. **Heritability correlation tests** — SNP additive effects vs species h².  
4. **Mediation + MR analyses** — host ↔ microbiome ↔ phenotype.  
5. **Structural Equation Modelling (SEM)** — integrate validated paths into a causal multi-layer model.  

---

## ⚙️ <b>Prerequisites</b>

### 🧰 Dependencies
You will need **Nextflow** and **Singularity/Apptainer**.

#### 📦 Install Nextflow
```bash
wget -qO- https://get.nextflow.io | bash
````

#### 📦 Install Singularity

Follow [these instructions](https://docs.sylabs.io/guides/3.1/user-guide/quick_start.html#quick-installation-steps)

---

### 💻 <b>CRG HPC Configuration</b>

Connect to the dedicated Nextflow node:

```bash
ssh -Y YOURUSER@nextflow.hpc.crg.es
```

Then prepare your environment:

```bash
mkdir -p $HOME/tmp $HOME/singularity_containers $HOME/tmp_singu
```

Add to `~/.bashrc`:

```bash
module use /software/as/el7.2/EasyBuild/CRG/modules/all
module load Java/11.0.2
module load Singularity/3.7.0
export SINGULARITY_CACHEDIR=$HOME/tmp
export NXF_SINGULARITY_CACHEDIR="$HOME/singularity_containers/"
export SINGULARITY_TMPDIR=$HOME/tmp_singu
```

Log out and back in, then install Nextflow.

---

## 📦 <b>Installing the Suite</b>

```bash
git clone --recurse-submodules git@github.com:Baud-lab/hermes-wire.git
```

---

## ⚙️ <b>Nextflow Configuration</b>

The <code>nextflow.config</code> file defines how processes run in parallel, setting memory/time requirements and cluster queues.
HERMES-WIRE ships with **SLURM-compatible profiles** for the CRG HPC cluster, but these can be adapted for other schedulers.

---

## ▶️ <b>Running a Pipeline (Example: Step 1)</b>

```bash
sbatch submit_nf.sh main.nf -profile singularity,slurm_genoa -params-file params.yaml -resume -w work
```

To stop a run:

```bash
cat .nextflow.pid | xargs kill
```

---

## 🧾 <b>Parameter File Example (Step 1)</b>

```yaml
module_preprocessing: "YES"
module_taxonomic_mapping: "YES"
module_functional_mapping: "YES"
tool_opt: "./tool_opt.tsv"
preprocessing_file_raw_reads: "./dataset_raw/demo/*_R{1,2}_001.fastq.gz"
preprocessing_file_host_genome: "./input/mRatBN7.2.fna"
taxonomic_file_microbiome_index: "./input/indices/kaiju/gtdb_index.fmi"
functional_file_index: "./input/indices/prromenade/bactvirus2020*"
output_folder: "./Results_Step1/"
```

---

## 🌐 <b>Monitoring with Nextflow Tower</b>

1. Log in at [https://tower.nf](https://tower.nf)
2. Generate a token → copy it
3. Add to your `~/.bashrc`:

   ```bash
   export TOWER_ACCESS_TOKEN=<YOUR_ACCESS_TOKEN>
   ```
4. Run with:

   ```bash
   nextflow run main.nf -with-tower ...
   ```

Tower provides real-time monitoring, job tracing, and report visualisation.

---

## 🧩 <b>Acknowledgements</b>

Developed in the <b>[Baud Lab](https://www.crg.eu/en/programmes-groups/baud-lab)</b> at the <b>[Centre for Genomic Regulation (CRG)](https://www.crg.eu)</b> and <b>[Universitat Pompeu Fabra (UPF)](https://www.upf.edu)</b>, Barcelona.
We acknowledge support from the <b>[Bioinformatics Core Facility](https://www.crg.eu/ca/programmes-groups/bioinformatics-unit)</b> and the <b>CRG [Scientific IT team](https://www.crg.eu/en/content/about-us-administration/scientific-information-technologies)</b>. Furthermore, the data used and tested for the development of this tool was generated in collaboration with <b>[NIDA](https://ratgenes.org)</b> and the <b>[Center for Microbiome Innovation](https://cmi.ucsd.edu)</b>. The project received the support of <b>[La Caixa Foundation](https://lacaixafoundation.org/en/)</b>.

---

## 🧠 <b>Citation</b>

> Morillo FMSD, Cozzuto L, Tonnelé H, Baud A (2025). <i>HERMES-WIRE: HERitable MicrobiomE Structure — Workflow for Interpreting host–microbiome Relationships & Effects.</i>
> Centre for Genomic Regulation (CRG), Barcelona.
> [https://github.com/Baud-lab/hermes-wire](https://github.com/Baud-lab/hermes-wire)

---

## 🧾 <b>License</b>

© 2025 Centre for Genomic Regulation (CRG) and the authors.
Distributed under the <b>[Apache License 2.0](https://github.com/Baud-lab/hermes-wire/blob/main/LICENSE)</b>.

---
