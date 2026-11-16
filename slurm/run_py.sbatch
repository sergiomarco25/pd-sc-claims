#!/bin/bash
#SBATCH --job-name=pdrep_py
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/run_py_%x_%j.out
#SBATCH --error=logs/run_py_%x_%j.err
# Generic runner (adapted from revision/slurm/run_py.sbatch). Usage: sbatch -J <name> slurm/run_py.sbatch <script.py> [args...]
# Submit from the repository root: log paths above are relative to it.
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
eval "$(python3 -I src/pdrep/config.py --env)"   # paths from config/project.toml
export TMPDIR=$PDREP_SCRATCH/tmp_$SLURM_JOB_ID
mkdir -p "$TMPDIR"; trap 'rm -rf "$TMPDIR"' EXIT
"$PDREP_PYTHON" -I "$@"
