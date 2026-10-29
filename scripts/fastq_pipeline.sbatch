#!/bin/bash -l
#SBATCH --partition=genoa64
#SBATCH --qos=pipelines
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --no-requeue
set -euo pipefail
stage=${1:?Usage: fastq_pipeline.sbatch rnaseq|vast|dge PROJECT_DIR}
export PROJECT_DIR=${2:?Missing project directory}
cd "$PROJECT_DIR"
source config/fastq_cluster.env
export RAW_DIR VAST_PIPELINE VASTDB SCRATCH_DIR GENOME_FASTA GENES_GTF
module load Nextflow
module load Java
module load Python/3.12.3-GCCcore-13.2.0
command -v singularity >/dev/null
export NXF_VER=25.10.4
export NXF_JVM_ARGS='-Xms1g -Xmx5g'
export NXF_OPTS='-XX:ActiveProcessorCount=2'
export TERM=dumb SLURM_EXPORT_ENV=ALL
export SLURM_CONF=${SLURM_CONF:-/etc/slurm/slurm.conf}
export NXF_SINGULARITY_CACHEDIR=$SCRATCH_DIR/singularity
export READONLY_BINDS="$PROJECT_DIR:$PROJECT_DIR:ro,$RAW_DIR:$RAW_DIR:ro"
if [[ $stage != vast ]]; then
    export READONLY_BINDS="$READONLY_BINDS,$GENOME_FASTA:$GENOME_FASTA:ro,$GENES_GTF:$GENES_GTF:ro"
fi
samples=$PROJECT_DIR/data/processed/spire_fastq/samples
outdir=$PROJECT_DIR/results/spire_fastq
mkdir -p "$NXF_SINGULARITY_CACHEDIR" "$PROJECT_DIR/runs/spire_fastq/$stage"
seqera_args=()
if [[ -n ${TOWER_ACCESS_TOKEN:-} ]]; then seqera_args=(-with-tower); fi
cd "$PROJECT_DIR/runs/spire_fastq/$stage"
common=(-c "$PROJECT_DIR/config/fastq_crg.config" \
    -work-dir "$SCRATCH_DIR/maria_MH_Adel_spire12ko_2026/work_$stage" \
    -resume -ansi-log false "${seqera_args[@]}")
case "$stage" in
    rnaseq)
        python3 -I "$PROJECT_DIR/scripts/prepare_fastq_samples.py" "$RAW_DIR" "$samples"
        exec nextflow run nf-core/rnaseq -r 3.27.0 -profile singularity \
            -params-file "$PROJECT_DIR/config/spire_rnaseq.json" \
            --input "$samples/rnaseq.csv" --fasta "$GENOME_FASTA" --gtf "$GENES_GTF" \
            --outdir "$outdir/rnaseq" "${common[@]}"
        ;;
    vast)
        python3 -I "$PROJECT_DIR/scripts/prepare_fastq_samples.py" "$RAW_DIR" "$samples" \
            --trimmed "$outdir/rnaseq/fastp"
        exec nextflow run "$VAST_PIPELINE/main.nf" \
            -c "$PROJECT_DIR/config/spire_vast.config" \
            --sample_csv "$samples/vast.csv" --data_dir "$outdir/rnaseq/fastp" \
            --vastdb_path "$VASTDB" --species mm10 --project_name spire_fastq \
            --outdir "$outdir/vast" --skip_trimming --skip_fastqc --skip_compare --skip_matt --skip_betas \
            "${common[@]}"
        ;;
    dge)
        counts=$outdir/rnaseq/star_salmon/salmon.merged.gene_counts.tsv
        lengths=$outdir/rnaseq/star_salmon/salmon.merged.gene_lengths.tsv
        for path in "$counts" "$lengths"; do
            if [[ ! -s $path ]]; then echo "Missing RNA-seq output: $path" >&2; exit 1; fi
        done
        exec nextflow run nf-core/differentialabundance -r 2.0.0 -profile rnaseq,singularity \
            --input "$PROJECT_DIR/config/spire_fastq.tsv" \
            --contrasts "$PROJECT_DIR/config/spire_contrasts.csv" \
            --matrix "$counts" --feature_length_matrix "$lengths" --gtf "$GENES_GTF" \
            --exploratory_main_variable condition \
            --study_name SPIRE_DKO --outdir "$outdir/dge" "${common[@]}"
        ;;
    *) echo "Unknown stage: $stage" >&2; exit 2 ;;
esac
