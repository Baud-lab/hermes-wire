#!/usr/bin/env bash
#SBATCH --no-requeue
#SBATCH --mem 6G
#SBATCH -p genoa64
#SBATCH --qos=pipelines
#SBATCH --time=03:00:00
#SBATCH -J submit_n

set -euo pipefail

# ---- live logging ----
#LAUNCH_LOG="log"
#exec > >(tee -a "$LAUNCH_LOG") 2>&1

_term() {
  echo "[submit_nf] Caught SIGTERM, forwarding to Nextflow ($pid)"
  kill -s SIGTERM $pid || true
  wait $pid || true
}
trap _term TERM

module load Java || true

# Limit JVM for the Nextflow *manager* process
export NXF_JVM_ARGS="-Xms2g -Xmx5g"

# Singularity cache on shared storage (fast + reused across nodes)
export NXF_SINGULARITY_CACHEDIR="/users/abaud/fmorillo/singularity_containers"
export SINGULARITY_CACHEDIR="$NXF_SINGULARITY_CACHEDIR"

echo "[submit_nf] Host: $(hostname)"
echo "[submit_nf] PWD : $(pwd)"
echo "[submit_nf] SLURM_JOBID: ${SLURM_JOB_ID:-N/A}"
echo "[submit_nf] Singularity cache: $SINGULARITY_CACHEDIR"

# Which Nextflow binary to use
NF_BIN="${NF_BIN:-$HOME/nextflow}"
echo "[submit_nf] Nextflow binary: ${NF_BIN}"
"$NF_BIN" -v || true

# Launch Nextflow — pass all args exactly as you provide them
echo "[submit_nf] Launching: $NF_BIN run -ansi-log false $*"
"$NF_BIN" run -ansi-log false "$@" & pid=$!

echo "[submit_nf] Waiting for Nextflow PID ${pid}"
wait $pid
rc=$?
echo "[submit_nf] Nextflow exit code: $rc"
exit $rc
