#!/bin/bash
#SBATCH --job-name=pdrep_r
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --output=logs/run_r_%x_%j.out
#SBATCH --error=logs/run_r_%x_%j.err
# Generic R runner. Usage: sbatch -J <name> [--export=RENV=<env path>] slurm/run_r.sbatch <script.R> [args...]
# Submit from the repository root: log paths above are relative to it.
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
eval "$(python3 -I src/pdrep/config.py --env)"   # paths from config/project.toml
RENV=${RENV:-$PDREP_R_ENV}
# node /tmp can be full ("cannot create 'R_TempDir'"): keep R temp files on lustre
export TMPDIR=$PDREP_SCRATCH/tmp_$SLURM_JOB_ID
mkdir -p "$TMPDIR"; trap 'rm -rf "$TMPDIR"' EXIT
"$RENV/bin/Rscript" --vanilla "$@"
