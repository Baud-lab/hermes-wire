# 🧬 HERMES-WIRE — Step 3: Genome functional annotation → function×species matrices → phylogeny-aware enrichment

**Centre for Genomic Regulation (CRG), Barcelona — 2025**  
**Authors:** Felipe Morillo Sanz Dias, Luca Cozzuto, Hélène Tonnele, Amelie Baud

---

## 📘 Overview

**Step 3** takes **protein FASTA** files (we used **GTDB r207 proteins**) and:

1) **Annotates** proteins with **eggNOG-mapper** (COG/ENOG, KEGG KO, EC, PFAM).  
2) **Builds matrices** of **functions × species** (from all species represented in the FASTAs).  
3) **Tests phylogeny-aware enrichment** (PGLS) of functions in **high- vs low-heritability** species using the GTDB taxonomy, phylogenetic tree and Step-2 heritability results.

> ℹ️ CAZy: The core pipeline here **does not run dbCAN2** itself. If you want CAZy, annotate externally (e.g. with `run_dbcan`) and ensure the parsed tables are available for `make_func_matrices.R` to ingest. See “Optional: CAZy” below.

---

## ⚙️ What this pipeline actually runs

- **Process 1 — EMMAPPER_ANNOTATE**  
  Runs eggNOG-mapper on each `*.faa` and writes `*.annotations.tsv`.

- **Process 2 — BUILD_MATRICES**  
  Reads all annotation TSVs and creates:
  - `func_matrices.RData` (all selected spaces)
  - `func_eggNOG_OGs.RData` (eggNOG OG-focused)
  - `qc_functions_summary.tsv`

- **Process 3 — ENRICHMENT_PGLS** *(enabled only if `taxonomy_tsv`, `phylo_tree`, `herit_rdata` are set)*  
  Runs phylogeny-aware enrichment (PGLS) across species and writes:
  - `res_*.csv` tables
  - `enrichment_*.RData`
  - figures `*.pdf`

---

## 📁 Directory layout (repo)

```

Step3/
├── main.nf
├── nextflow.config
├── params.yaml
├── submit_nf.sh
├── bin/                     # helper scripts called by processes
├── input/                   # e.g., cog-24.def.tab, small demos
├── dockerfiles/             # optional
├── singularity_containers/  # optional
└── work/                    # nf work dir (auto)

````

---

## 🔧 Inputs (params)

These mirror `main.nf` and `params.yaml`.

**Core:**
- `input_faa_glob` : Glob for protein FASTAs (e.g. `input/*.faa`)
- `threads` : CPU threads per eggNOG-mapper job
- `emapper_mode` : `"diamond"` (default) or `"mmseqs"`
- `emapper_args_extra` : Extra args to pass to eggNOG-mapper (optional)

**Optional (for enrichment):**
- `taxonomy_tsv` : GTDB mapping (`Genome\tTaxonomy`)
- `phylo_tree` : GTDB r207 Newick tree
- `herit_rdata` : RData with heritability objects (from Step 2)

**Function space building:**
- `functions` : Comma-sep list (default:
  `eggNOG_OGs,CAZy,EC,KEGG_ko,KEGG_Module,KEGG_Pathway,PFAMs,COG_category,COG_pathway`)
- `cog_def_tab` : `cog-24.def.tab` (functional names for COGs)

**Enrichment options:**
- `subset_id` : Label for subset (e.g., `ALL`)
- `method` : e.g., `Shallow`
- `report_model` : `ML` or `AIC`
- `branchlen_transform` : `none` or `grafen`
- `groups` : Taxa of interest (e.g., `g__Prevotella,g__Bacteroides`)

**Containers (optional):**
- `container_emapper` : e.g., `docker://eggnogmapper/emapper:2.1.12`
- `container_R` : e.g., `docker://.../genesis-nf:bioc3.18-ext.v2`

---

## 🚀 How to run

### A) With the provided `params.yaml`

```bash
nextflow run main.nf -params-file params.yaml -with-report -resume
````

### B) On a cluster (example wrapper)

```bash
bash submit_nf.sh
```

> Adjust `nextflow.config` (profiles, queue, resources) as needed.

---

## 🧾 Outputs

```
Results/
├── 01_annotations/
│   └── <GENOME>.annotations.tsv
├── 02_matrices/
│   ├── func_matrices.RData
│   ├── func_eggNOG_OGs.RData
│   └── qc_functions_summary.tsv
└── 03_enrichment/   # only if taxonomy + tree + herit provided
    ├── res_*.csv
    ├── enrichment_*.RData
    └── *.pdf
```

**Final matrices are** **functions × species** (matching the set of species present in the protein FASTAs).

---

## 🧪 Notes on function spaces

* **eggNOG-mapper** supplies COG/ENOG, KEGG KO, EC, PFAM.
* **CAZy** requires external annotation (e.g., `run_dbcan`) or a curated mapping layer. If you include CAZy in `functions` but don’t supply CAZy hits, those slots will be empty.

---

## 🌳 Phylogeny-aware enrichment (PGLS)

If you provide:

* `taxonomy_tsv` (GTDB mapping),
* `phylo_tree` (GTDB r207 Newick),
* `herit_rdata` (from Step 2; includes heritability metrics),

then **ENRICHMENT_PGLS** runs PGLS across species and writes CSV tables and PDF figures in `03_enrichment/`.

---

## 📦 Optional: fetch genomes from NCBI and produce proteins

> We used **GTDB r207 protein FASTAs** for our analyses. If you want to start from NCBI accessions instead:

### 1) Download assemblies (GenBank or RefSeq)

```bash
# accessions.txt contains one GCA_/GCF_ per line
ncbi-genome-download --assembly-accessions accessions.txt \
  -s genbank -F fasta -o ./ncbi_dl bacteria   # or -s refseq
```

### 2) Collect FASTA files into one place

```bash
mkdir -p genomes
find ./ncbi_dl -type f -name "*.fna.gz" -exec cp {} genomes/ \;
gunzip -f genomes/*.fna.gz
```

### 3) Call proteins with Prodigal

```bash
mkdir -p proteins
for f in genomes/*.fna; do
  base=$(basename "$f" .fna)
  prodigal -i "$f" -a "proteins/${base}.faa" -q
done
```

Then set `input_faa_glob: ./proteins/*.faa` in `params.yaml` and run Step 3.

---

## (Optional) CAZy with dbCAN2

If you want **CAZy**:

```bash
# Example per-species protein file
run_dbcan proteins/<species>.faa protein --out_dir cazy_out/<species> --cpu 8
```

Parse/aggregate dbCAN outputs and ensure `make_func_matrices.R` sees them (e.g., place compatible TSVs into the annotation folder the script ingests). See dbCAN2 docs for usage details.

---

## 🧩 <b>Acknowledgements</b>

Developed under the <b>[Baud Lab](https://www.crg.eu/en/programmes-groups/baud-lab)</b> at the <b>[Centre for Genomic Regulation (CRG)](https://www.crg.eu)</b> and <b>[Universitat Pompeu Fabra (UPF)](https://www.upf.edu)</b>, Barcelona.
We acknowledge support from the <b>[Bioinformatics Core Facility](https://www.crg.eu/ca/programmes-groups/bioinformatics-unit)</b> and the <b>CRG [Scientific IT team](https://www.crg.eu/en/content/about-us-administration/scientific-information-technologies)</b>. Furthermore, the data used and tested for the development of this tool was generated in collaboration with <b>[NIDA](https://ratgenes.org)</b> and the <b>[Center for Microbiome Innovation](https://cmi.ucsd.edu)</b>, with the support of <b>[La Caixa Foundation](https://lacaixafoundation.org/en/)</b>.

---

## 🧠 <b>Citation</b>

> Morillo FMSD, Cozzuto L, Tonnele H, Baud A (2025). <i>HERMES-WIRE: HERitable MicrobiomE Structure — Workflow for Interpreting host–microbiome Relationships & Effects.</i>
> Centre for Genomic Regulation (CRG), Barcelona.
> [https://github.com/Baud-lab/hermes-wire](https://github.com/Baud-lab/hermes-wire)

---

## 🧾 <b>License</b>

© 2025 Centre for Genomic Regulation (CRG) and the authors.
Distributed under the <b>[Apache License 2.0](https://github.com/Baud-lab/hermes-wire/blob/main/LICENSE)</b>.
