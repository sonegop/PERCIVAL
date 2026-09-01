#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --output=/mnt/nfs3/sonegop/projects/PERCIVAL/logs/fastqc_fastp_%x_%A_%a.out
#SBATCH --error=/mnt/nfs3/sonegop/projects/PERCIVAL/logs/fastqc_fastp_%x_%A_%a.err
#SBATCH --time=1:00:00          # Maximum execution time (4 hours)
#SBATCH --ntasks=1              # Number of tasks (processes)
#SBATCH --cpus-per-task=12       # CPUs per task
#SBATCH --mem=4G                # Memory allocation
#SBATCH --partition=cpu
#SBATCH --array=0-119           # Adjust to number of FASTQ files minus one

# ==============================
# USER PARAMETERS
# ==============================
#FASTQ_DIR=/mnt/nfs3/sonegop/projects/PERCIVAL/data/fastq/raw  # directory containing fastqs
FASTQ_DIR=/mnt/nfs3/sonegop/projects/PERCIVAL/data/fastq/fastq  # directory containing fastqs
#OUT_DIR=/mnt/nfs3/sonegop/projects/PERCIVAL/data/fastqc_results/raw       # output directory
OUT_DIR=/mnt/nfs3/sonegop/projects/PERCIVAL/data/fastqc_results/fastq      # output directory
LOG_DIR=/mnt/nfs3/sonegop/projects/PERCIVAL/logs
CONTAINER=/mnt/nfs3/sonegop/images/apptainer/fastqc_0.12.1--hdfd78af_0.sif  # FastQC container (adjust path if needed)

# ==============================
# PREPARE ENV
# ==============================
# module load apptainer   # only if required by your cluster

mkdir -p $OUT_DIR
mkdir -p $LOG_DIR

# ==============================
# CREATE SAMPLE LIST
# ==============================
# Generate a list of all FASTQ files (both R1 and R2)
SAMPLE_LIST=$(ls ${FASTQ_DIR}/*.fq.gz | sort)
# Convert to array
FASTQS=($SAMPLE_LIST)

# Select file for this array task
fastq=${FASTQS[$SLURM_ARRAY_TASK_ID]}
sample=$(basename $fastq .fq.gz)

echo "Running FastQC on $sample (task $SLURM_ARRAY_TASK_ID)"

# ==============================
# RUN FASTQC
# ==============================
echo "Starting FastQC for sample $sample..."
apptainer exec \
    --bind $FASTQ_DIR:/data \
    --bind $OUT_DIR:/out \
    $CONTAINER \
    fastqc \
        --threads $SLURM_CPUS_PER_TASK \
        --outdir /out \
        /data/$(basename $fastq)

echo "FastQC completed for sample $sample"
