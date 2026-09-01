#!/bin/bash
#SBATCH --job-name=align_star
#SBATCH --output=/mnt/nfs3/sonegop/projects/PERCIVAL/logs/align_star_%x_%A_%a.out
#SBATCH --error=/mnt/nfs3/sonegop/projects/PERCIVAL/logs/align_star_%x_%A_%a.err
#SBATCH --time=6:00:00           # Adjusted for STAR (typically faster than HISAT2)
#SBATCH --ntasks=1                # STAR uses threads, not MPI tasks
#SBATCH --cpus-per-task=18        # STAR benefits from more threads
#SBATCH --mem=42G                 # STAR needs more RAM for genome loading
#SBATCH --partition=cpu
#SBATCH --array=0-119             # Adjust to number of samples minus one

# ==============================
# USER PARAMETERS
# ==============================
GENOME_INDEX=/mnt/nfs3/sonegop/projects/PERCIVAL/references/selected_references/indexes/Slycopersicum_SL5.0_ITAG5.0_STAR
ANNOTATION=/mnt/nfs3/sonegop/projects/PERCIVAL/references/selected_references/Slycopersicum_796_ITAG5.0.gene_exons.gff3
FASTQ_DIR=/mnt/nfs3/sonegop/projects/PERCIVAL/data/fastq/fastp
OUT_DIR=/mnt/nfs3/sonegop/projects/PERCIVAL/data/star_results/fastp
LOG_DIR=/mnt/nfs3/sonegop/projects/PERCIVAL/logs
CONTAINER=/mnt/nfs3/sonegop/images/apptainer/fastp_star_samtools_subread.sif

# STAR-specific parameters
READ_LENGTH=150                   # Read length for your data
MAX_INTRON=20000                  # Maximum intron size for Vitis
MIN_INTRON=20                     # Minimum intron size
OVERHANG=$((READ_LENGTH - 1))     # Splice junction overhang (typically readLength - 1)

# Multimapping parameters
MULTIMAP_MAX=20                   # Maximum number of loci a read is allowed to map to
OUT_FILTER_TYPE="BySJout"         # Reduces spurious junctions
OUT_SAM_MULTIMAP_MAX=1            # Max number of alignments reported for multimappers (-1 for all)

# ==============================
# PREPARE ENVIRONMENT
# ==============================
# module load apptainer   # Only if required by your cluster

mkdir -p $OUT_DIR
mkdir -p $LOG_DIR
mkdir -p ${OUT_DIR}/bam
mkdir -p ${OUT_DIR}/logs
mkdir -p ${OUT_DIR}/junctions

# ==============================
# CREATE SAMPLE LIST
# ==============================
# Generate a list of R1 FASTQ files, one per line
SAMPLE_LIST=$(ls ${FASTQ_DIR}/*clean_1.fq.gz | sort)
# Convert to array
FASTQS=($SAMPLE_LIST)

# Select files for this array task
fq1=${FASTQS[$SLURM_ARRAY_TASK_ID]}
fq2=${fq1/clean_1.fq.gz/clean_2.fq.gz}
sample=$(basename $fq1 clean_1.fq.gz)

echo "================================================"
echo "Running STAR alignment for sample: $sample"
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Read 1: $(basename $fq1)"
echo "Read 2: $(basename $fq2)"
echo "================================================"

# Create sample-specific output directory
SAMPLE_OUT=${OUT_DIR}/${sample}
mkdir -p ${SAMPLE_OUT}

# ==============================
# RUN STAR ALIGNMENT
# ==============================
apptainer exec \
    --bind $FASTQ_DIR:/data \
    --bind $OUT_DIR:/out \
    --bind $GENOME_INDEX:/genome \
    --bind $(dirname $ANNOTATION):/annotation \
    $CONTAINER \
    STAR \
        --runThreadN $SLURM_CPUS_PER_TASK \
        --genomeDir /genome \
        --readFilesIn /data/$(basename $fq1) /data/$(basename $fq2) \
        --readFilesCommand zcat \
        --outFileNamePrefix /out/${sample}/${sample}_ \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes NH HI AS nM MD \
        --outFilterType $OUT_FILTER_TYPE \
        --outFilterMultimapNmax $MULTIMAP_MAX \
        --outMultimapperOrder Random \
        --outSAMmultNmax $OUT_SAM_MULTIMAP_MAX \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --alignIntronMin $MIN_INTRON \
        --alignIntronMax $MAX_INTRON \
        --alignMatesGapMax 1000000 \
        --sjdbGTFfile /annotation/$(basename $ANNOTATION) \
        --sjdbGTFtagExonParentTranscript Parent \
        --sjdbGTFtagExonParentGene ID \
        --sjdbOverhang $OVERHANG \
        --quantMode GeneCounts \
        --twopassMode Basic \
        --limitBAMsortRAM 40000000000

# ==============================
# POST-PROCESSING
# ==============================
echo "================================================"
echo "Post-processing for sample: $sample"
echo "================================================"

# Index the sorted BAM file
apptainer exec \
    --bind $OUT_DIR:/out \
    $CONTAINER \
    samtools index /out/${sample}/${sample}_Aligned.sortedByCoord.out.bam

# Generate alignment statistics
apptainer exec \
    --bind $OUT_DIR:/out \
    $CONTAINER \
    bash -c "
        samtools flagstat /out/${sample}/${sample}_Aligned.sortedByCoord.out.bam > /out/${sample}/${sample}_flagstat.txt
        samtools stats /out/${sample}/${sample}_Aligned.sortedByCoord.out.bam > /out/${sample}/${sample}_stats.txt
    "

# Move key files to organized locations
mv ${SAMPLE_OUT}/${sample}_Aligned.sortedByCoord.out.bam ${OUT_DIR}/bam/${sample}.bam
mv ${SAMPLE_OUT}/${sample}_Aligned.sortedByCoord.out.bam.bai ${OUT_DIR}/bam/${sample}.bam.bai
mv ${SAMPLE_OUT}/${sample}_Log.final.out ${OUT_DIR}/logs/${sample}_final.log
mv ${SAMPLE_OUT}/${sample}_SJ.out.tab ${OUT_DIR}/junctions/${sample}_SJ.tab

# Clean up temporary files if needed (optional)
# rm -f ${SAMPLE_OUT}/${sample}_Log.progress.out
# rm -f ${SAMPLE_OUT}/${sample}_Log.out

echo "================================================"
echo "Finished processing sample: $sample"
echo "Output BAM: ${OUT_DIR}/bam/${sample}.bam"
echo "Final log: ${OUT_DIR}/logs/${sample}_final.log"
echo "================================================"
