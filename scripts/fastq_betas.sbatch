#!/bin/bash -l
#SBATCH --job-name=spire_betas
#SBATCH --partition=genoa64
#SBATCH --qos=vlong
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --no-requeue
set -euo pipefail
export PROJECT_DIR=${1:?Usage: fastq_betas.sbatch PROJECT_DIR}
cd "$PROJECT_DIR"
source config/fastq_cluster.env
export SPIRE_FASTQ_INCLUSION_TABLE BETAS_DATASET=spire_fastq
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export TMPDIR=$SCRATCH_DIR/maria_MH_Adel_spire12ko_2026/betas/${SLURM_JOB_ID:?Run in a Slurm allocation}
mkdir -p "$TMPDIR"
run_dir=results/spire_fastq/betas_N10_nsim1000
if [[ ! -f $run_dir/run.rds ]]; then
    bash scripts/run_R.sh scripts/01_filter.R
    bash scripts/run_R.sh scripts/00_init_run.R "$run_dir"
fi
bash scripts/run_R.sh scripts/02_simulate.R "$run_dir"
exec bash scripts/run_R.sh scripts/03_collect.R "$run_dir"
