#!/bin/bash
#SBATCH --job-name=sem_ff
#SBATCH --open-mode=append
#SBATCH --time=00:30:00
#SBATCH --mem=32G
#SBATCH --qos=vshort

set -euo pipefail
echo "[$(date)] starting node=$(hostname)"

# A container with all required packages in sem_feedfrward.R must run from this point on with ```singularity run```

CODE="./sem_feedforward.R"

# --- inputs (same as before) ---
COHORT="../Step4/Results_bacteriodes/cohort_ids/genesis_final_label_ids.txt" # Output obtained in Step 4
MED_BAC="../Step4/Results_bacteriodes/phase9_mediation/mediation_sig_iv.tsv" # Output obtained in Step 4
MED_PREV="../Step4/Results_prevotella/phase9_mediation/mediation_sig_iv.tsv" # Output obtained in Step 4
PUTATIVE="../Step4/input/putative_genes_final.xlsx" # Same input as in Step 4
GDS="../Step4/Results_bacteriodes/phase0/genotypes_candidates.gds" # Output obtained in Step 4

RESID_RDA="../Step2/Results_Taxonomic/Matrix_Processing/residuals_qned_counts.RData" # Output obtained in Step 2
BETA_RDA="../Step2/Results_Taxonomic/Diversity/residuals_qned_counts_beta.RData" # Output obtained in Step 2
ALPHA_RDA="../Step2/Results_Taxonomic/Diversity/residuals_qned_counts_alpha.RData" # Output obtained in Step 2
CLUSTERS_RDA="../Step2/Results_Taxonomic/Cluster_Analyses/residuals_qned_counts_clusters.RData" # Output obtained in Step 2
PHEN_RDA="../Step4/input/pehotypes_residuals.RData" # Same input as in Step 4
GRM_RDA="../Step2/input/input/grm.RData" # Same input as in Step 2
META="../Step2/input/input/metadata_Shallow.txt" # Same input as in Step 2

OUTDIR="../Results"

# ---- optional custom model control ----
# Provide EITHER MODEL_FILE or MODEL_TEXT (leave both empty to use the feed-forward builder)
MODEL_FILE="../input/my_model.lav"
cat ${MODEL_FILE}
MODEL_TEXT=""   # e.g. $'Bact ~~ Prev\nBeta ~~ Alpha\nGuild ~ Bact + Prev\n...'

# Only knob for inclusion/exclusion:
#   Empty  -> use all variables available in data
#   Non-empty -> intersect with available ones
COMPONENTS="Bact, Prev, Guild, Beta, Alpha, Glucose" # e.g. "Bact, Prev, Guild, Beta, Alpha, Glucose, BMI"

# ---- run ----
srun -u stdbuf -oL -eL Rscript --vanilla "${CODE}" \
  --cohort "${COHORT}" \
  --med_bac_tsv "${MED_BAC}" \
  --med_prev_tsv "${MED_PREV}" \
  --putative_genes_xlsx "${PUTATIVE}" \
  --selected_genes_rdata "${GDS}" \
  --resid_rda "${RESID_RDA}" \
  --beta_rda "${BETA_RDA}" \
  --alpha_rda "${ALPHA_RDA}" \
  --clusters_rda "${CLUSTERS_RDA}" \
  --phen_rda "${PHEN_RDA}" \
  --grm_rda "${GRM_RDA}" \
  --metadata "${META}" \
  --meta_rfid_col RFID \
  --meta_dam_col dam \
  --meta_cage_col cage \
  --bac_name g__Bacteroides \
  --prev_name g__Prevotella \
  --guild_name cluster__3 \
  --beta_trait beta__PD_PC1 \
  --alpha_trait alpha__PD_q2 \
  --glucose_col Glucose \
  --bmi_col BMI \
  --drop_genes FALSE \
  --iv_gold_tier FALSE \
  --missing_mode fiml \
  --drop_corr_threshold 0.99 \
  --qr_drop TRUE \
  --jitter_sd 0.0 \
  --print_mmer_summary TRUE \
  --fdr_level 0.10 \
  --mi_threshold 10.0 \
  --max_iter 10 \
  --model_file "${MODEL_FILE}" \
  --model_text "${MODEL_TEXT}" \
  --components "${COMPONENTS}" \
  --outdir "${OUTDIR}"

echo "[$(date)] finished"
