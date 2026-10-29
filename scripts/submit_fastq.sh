#!/usr/bin/env bash
set -euo pipefail
stage=${1:-all}
case "$stage" in
    all|rnaseq|vast|dge|betas) ;;
    *) echo 'Usage: bash scripts/submit_fastq.sh [all|rnaseq|vast|dge|betas]' >&2; exit 2 ;;
esac
cd "$(dirname "${BASH_SOURCE[0]}")/.."
export PROJECT_DIR=$PWD
source config/fastq_cluster.env
export RAW_DIR VAST_PIPELINE VASTDB SCRATCH_DIR GENOME_FASTA GENES_GTF SPIRE_FASTQ_INCLUSION_TABLE
command -v sbatch >/dev/null
if [[ $stage == all || $stage == rnaseq || $stage == vast ]]; then
    python3 -I scripts/prepare_fastq_samples.py "$RAW_DIR" data/processed/spire_fastq/samples
fi
if [[ $stage == all || $stage == rnaseq || $stage == dge ]]; then
    for path in "$GENOME_FASTA" "$GENES_GTF"; do
        if [[ ! -s $path ]]; then echo "Missing reference: $path" >&2; exit 1; fi
        printf 'Using RNA-seq reference: %s\n' "$path"
    done
fi
if [[ $stage == betas && ! -s $SPIRE_FASTQ_INCLUSION_TABLE ]]; then
    echo "Missing VAST inclusion table: $SPIRE_FASTQ_INCLUSION_TABLE" >&2
    exit 1
fi
if [[ $stage == all || $stage == vast ]]; then
    for path in "$VAST_PIPELINE/main.nf" "$VAST_PIPELINE/nextflow.config"; do
        if [[ ! -s $path ]]; then echo "Missing VAST workflow: $path" >&2; exit 1; fi
    done
    if [[ ! -s $VASTDB/Mm2/FILES/gDNA.1.ebwt ]]; then
        echo "Missing or incomplete mm10 VAST library: $VASTDB/Mm2" >&2
        exit 1
    fi
    printf 'Using VASTDB: %s\n' "$VASTDB/Mm2"
fi
mkdir -p logs runs/spire_fastq
dependency=()
if [[ $stage == all || $stage == rnaseq ]]; then
    job=$(sbatch --parsable --export=ALL --job-name=spire_rnaseq \
        --output="$PROJECT_DIR/logs/spire_rnaseq.%j.log" \
        scripts/fastq_pipeline.sbatch rnaseq "$PROJECT_DIR")
    printf 'rnaseq\t%s\n' "$job" | tee -a runs/spire_fastq/slurm_jobs.tsv
    dependency=(--dependency="afterok:${job%%;*}")
fi
for downstream in vast dge; do
    if [[ $stage == all || $stage == "$downstream" ]]; then
        job=$(sbatch --parsable --export=ALL --job-name="spire_$downstream" "${dependency[@]}" \
            --output="$PROJECT_DIR/logs/spire_$downstream.%j.log" \
            scripts/fastq_pipeline.sbatch "$downstream" "$PROJECT_DIR")
        printf '%s\t%s\n' "$downstream" "$job" | tee -a runs/spire_fastq/slurm_jobs.tsv
        if [[ $downstream == vast ]]; then vast_job=${job%%;*}; fi
    fi
done
if [[ $stage == all || $stage == betas ]]; then
    dependency=()
    if [[ $stage == all ]]; then dependency=(--dependency="afterok:$vast_job"); fi
    job=$(sbatch --parsable --export=ALL --job-name=spire_betas "${dependency[@]}" \
        --output="$PROJECT_DIR/logs/spire_betas.%j.log" \
        scripts/fastq_betas.sbatch "$PROJECT_DIR")
    printf 'betas\t%s\n' "$job" | tee -a runs/spire_fastq/slurm_jobs.tsv
fi
