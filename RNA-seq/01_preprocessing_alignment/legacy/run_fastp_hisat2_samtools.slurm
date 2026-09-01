#!/bin/bash
#SBATCH --job-name=fastp_hisat2
#SBATCH --output=/mnt/nfs3/sonegop/projects/PERCIVAL/logs/fastp_hisat2_%x_%A_%a.out
#SBATCH --error=/mnt/nfs3/sonegop/projects/PERCIVAL/logs/fastp_hisat2_%x_%A_%a.err
#SBATCH --time=12:00:00 # Maximum execution time (12 hours)
#SBATCH --ntasks=1               # Number of tasks (processes)
#SBATCH --cpus-per-task=16
#SBATCH --mem=28G                 # Increased memory for fastp + hisat2
#SBATCH --partition=cpu
#SBATCH --array=0-119   # adjust to number of samples minus one: here 100/2-1

# ==============================
# USER PARAMETERS
# ==============================
GENOME_INDEX=/mnt/nfs3/sonegop/projects/PERCIVAL/references/selected_references/indexes/Slycopersicum_SL5.0_ITAG5.0_hisat2/Slycopersicum_SL5.0_ITAG5.0_hisat2
ANNOTATION=/mnt/nfs3/sonegop/projects/PERCIVAL/references/selected_references/Slycopersicum_796_ITAG5.0.gene_exons.gff3
FASTQ_DIR=/mnt/nfs3/sonegop/projects/PERCIVAL/data/fastq/raw
FASTP_OUT_DIR=/mnt/nfs3/sonegop/projects/PERCIVAL/data/fastq/fastp # fastp output directory
HISAT2_OUT_DIR=/mnt/nfs3/sonegop/projects/PERCIVAL/data/hisat2_results/fastp
LOG_DIR=/mnt/nfs3/sonegop/projects/PERCIVAL/logs
CONTAINER=/mnt/nfs3/sonegop/images/apptainer/fastp_hisat2_samtools_subread.sif 


# FastP parameters
MIN_READ_LENGTH=50      # Minimum read length after trimming
QUALITY_THRESHOLD=20    # Quality threshold for base trimming
WINDOW_SIZE=4          # Sliding window size for quality filtering

# HISAT2 parameters
MIN_INTRON=20          # Minimum intron length
MAX_INTRON=15000      # Maximum intron length

# ==============================
# PREPARE ENV
# ==============================
# module load apptainer   # only if required by your cluster

mkdir -p $FASTP_OUT_DIR
mkdir -p $HISAT2_OUT_DIR
mkdir -p $LOG_DIR

# ==============================
# CREATE SAMPLE LIST
# ==============================
# We generate a list of R1 FASTQ files, one per line
SAMPLE_LIST=$(ls ${FASTQ_DIR}/*_1.fq.gz | sort)
# Convert to array
FASTQS=($SAMPLE_LIST)

# Select file for this array task
fq1=${FASTQS[$SLURM_ARRAY_TASK_ID]}
fq2=${fq1/_1.fq.gz/_2.fq.gz}
sample=$(basename $fq1 _1.fq.gz)

echo "Processing sample $sample (task $SLURM_ARRAY_TASK_ID)"
echo "Raw R1: $fq1"
echo "Raw R2: $fq2"

# Check if R2 file exists
if [ ! -f "$fq2" ]; then
    echo "ERROR: R2 file not found: $fq2"
    exit 1
fi

# Define cleaned fastq file names
clean_fq1=${FASTP_OUT_DIR}/${sample}_clean_1.fq.gz
clean_fq2=${FASTP_OUT_DIR}/${sample}_clean_2.fq.gz
fastp_html=${FASTP_OUT_DIR}/${sample}_fastp_report.html
fastp_json=${FASTP_OUT_DIR}/${sample}_fastp_report.json

echo "Clean R1: $clean_fq1"
echo "Clean R2: $clean_fq2"

# ==============================
# RUN FASTP PREPROCESSING
# ==============================
echo "Starting FastP preprocessing for sample $sample..."

apptainer exec \
    --bind $FASTQ_DIR:/raw_data \
    --bind $FASTP_OUT_DIR:/clean_data \
    --bind $HISAT2_OUT_DIR:/alignment_out \
    --bind $(dirname $GENOME_INDEX):/ref \
    $CONTAINER \
    bash -c "
        echo 'Running FastP...'
        fastp \
            --in1 /raw_data/$(basename $fq1) \
            --in2 /raw_data/$(basename $fq2) \
            --out1 /clean_data/$(basename $clean_fq1) \
            --out2 /clean_data/$(basename $clean_fq2) \
            --html /clean_data/$(basename $fastp_html) \
            --json /clean_data/$(basename $fastp_json) \
            --thread $SLURM_CPUS_PER_TASK \
            --length_required $MIN_READ_LENGTH \
            --qualified_quality_phred $QUALITY_THRESHOLD \
            --cut_window_size $WINDOW_SIZE \
            --cut_mean_quality $QUALITY_THRESHOLD \
            --trim_poly_g \
            --trim_poly_x \
            --correction \
            --detect_adapter_for_pe
    "

# Check if fastp completed successfully
if [ $? -ne 0 ]; then
    echo "ERROR: FastP failed for sample $sample"
    exit 1
fi

echo "FastP preprocessing completed for sample $sample"

# ==============================
# RUN HISAT2 ALIGNMENT
# ==============================
echo "Starting HISAT2 alignment for sample $sample..."

apptainer exec \
    --bind $FASTQ_DIR:/raw_data \
    --bind $FASTP_OUT_DIR:/clean_data \
    --bind $HISAT2_OUT_DIR:/alignment_out \
    --bind $(dirname $GENOME_INDEX):/ref \
    $CONTAINER \
    bash -c "
        echo 'Running HISAT2 alignment...'
        hisat2 -x /ref/$(basename $GENOME_INDEX) \
            --very-sensitive \
            --seed 42 \
            --min-intronlen $MIN_INTRON \
            --max-intronlen $MAX_INTRON \
            --novel-splicesite-outfile /alignment_out/${sample}_novel_splicesites.txt \
            --summary-file /alignment_out/${sample}_alignment_summary.txt \
            -1 /clean_data/$(basename $clean_fq1) \
            -2 /clean_data/$(basename $clean_fq2) \
            -p $SLURM_CPUS_PER_TASK \
            -S /alignment_out/${sample}.sam
    "

# Check if hisat2 completed successfully
if [ $? -ne 0 ]; then
    echo "ERROR: HISAT2 alignment failed for sample $sample"
    exit 1
fi

echo "HISAT2 alignment completed for sample $sample"

# ==============================
# CONVERT SAM TO BAM AND SORT
# ==============================
echo "Converting SAM to sorted BAM for sample $sample..."

apptainer exec \
    --bind $HISAT2_OUT_DIR:/alignment_out \
    $CONTAINER \
    bash -c "
        echo 'Converting SAM to sorted BAM...'
        samtools view -@ $SLURM_CPUS_PER_TASK -bS /alignment_out/${sample}.sam | \
        samtools sort -@ $SLURM_CPUS_PER_TASK -o /alignment_out/${sample}_sorted.bam -
        
        echo 'Indexing BAM file...'
        samtools index /alignment_out/${sample}_sorted.bam
        
        echo 'Generating alignment statistics...'
        samtools flagstat /alignment_out/${sample}_sorted.bam > /alignment_out/${sample}_flagstat.txt
        samtools stats /alignment_out/${sample}_sorted.bam > /alignment_out/${sample}_stats.txt
    "

# Check if SAM to BAM conversion completed successfully
if [ $? -ne 0 ]; then
    echo "ERROR: SAM to BAM conversion failed for sample $sample"
    exit 1
fi

# ==============================
# CLEANUP
# ==============================
echo "Cleaning up intermediate files for sample $sample..."

# Remove SAM file to save space (optional - comment out if you want to keep it)
rm ${HISAT2_OUT_DIR}/${sample}.sam

echo "Pipeline completed successfully for sample $sample"
echo "Output files:"
echo "  - Cleaned reads: $clean_fq1, $clean_fq2"
echo "  - FastP report: $fastp_html"
echo "  - Aligned BAM: ${HISAT2_OUT_DIR}/${sample}_sorted.bam"
echo "  - Alignment stats: ${HISAT2_OUT_DIR}/${sample}_flagstat.txt"
