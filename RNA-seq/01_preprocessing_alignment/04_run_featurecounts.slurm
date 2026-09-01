#!/bin/bash
#SBATCH --job-name=subread
#SBATCH --output=/mnt/nfs3/sonegop/projects/PERCIVAL/logs/subread_%x_%A_%a.out
#SBATCH --error=/mnt//nfs3/sonegop/projects/PERCIVAL/logs/subread_%x_%A_%a.err
#SBATCH --time=6:00:00 # Tempo massimo di esecuzione (12 ore)
#SBATCH --ntasks=1               # Numero di task (processi)
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --partition=cpu
#SBATCH --array=0-119   # adjust to number of samples minus one: here 100/2-1

# ==============================
# USER PARAMETERS
# ==============================
SUBREAD_INDEX=/mnt/nfs3/sonegop/projects/PERCIVAL/references/PN40024/subread_index/subread_T2T_idx   # built Subread index prefix
ANNOTATION=/mnt/nfs3/sonegop/projects/PERCIVAL/references/selected_references/Slycopersicum_796_ITAG5.0.gene_exons.gff3  # GFF3 annotation file
FASTQ_DIR=/mnt/nfs3/sonegop/projects/PERCIVAL/data/fastq/fastp  # directory containing fastqs
OUT_DIR=/mnt/nfs3/sonegop/projects/PERCIVAL/data/subread_results  # output directory
LOG_DIR=/mnt/nfs3/sonegop/projects/PERCIVAL/logs
CONTAINER=/mnt/nfs3/sonegop/images/apptainer/hisat2.2.1_samtools_1.22.1_subread_2.1.1.sif
# ==============================
# PREPARE ENV
# ==============================
# module load apptainer   # only if required by your cluster

mkdir -p $OUT_DIR
mkdir -p $LOG_DIR

# ==============================
# CREATE SAMPLE LIST
# ==============================
# We generate a list of R1 FASTQ files, one per line (here *_1.fq.gz)
SAMPLE_LIST=$(ls ${FASTQ_DIR}/*_1.fq.gz | sort)
# Convert to array
FASTQS=($SAMPLE_LIST)

# Select file for this array task
fq1=${FASTQS[$SLURM_ARRAY_TASK_ID]}
fq2=${fq1/_1.fq.gz/_2.fq.gz}
sample=$(basename $fq1 _1.fq.gz)

echo "Running sample $sample (task $SLURM_ARRAY_TASK_ID)"

# ==============================
# RUN ALIGNMENT WITH SUBREAD
# ==============================
echo "Starting alignment for sample $sample..."
apptainer exec \
    --bind $FASTQ_DIR:/data \
    --bind $OUT_DIR:/out \
    --bind $(dirname $SUBREAD_INDEX):/ref \
    --bind $(dirname $ANNOTATION):/annotation \
    $CONTAINER \
    subread-align \
        -T $SLURM_CPUS_PER_TASK \
        -t 0 \
        -i /ref/$(basename $SUBREAD_INDEX) \
        -r /data/$(basename $fq1) \
        -R /data/$(basename $fq2) \
        -a /annotation/$(basename $ANNOTATION) \
        --gtfAttr Parent \
        -o /out/${sample}.bam

echo "Alignment completed for sample $sample"

# Index the BAM file
echo "Indexing BAM file for sample $sample..."
apptainer exec \
    --bind $OUT_DIR:/out \
    $CONTAINER \
    samtools index /out/${sample}.bam

echo "Indexing completed for sample $sample"

# ==============================
# RUN FEATURE COUNTING
# ==============================
echo "Starting feature counting for sample $sample..."
apptainer exec \
    --bind $OUT_DIR:/out \
    --bind $(dirname $ANNOTATION):/annotation \
    $CONTAINER \
    featureCounts \
        -a /annotation/$(basename $ANNOTATION) \
        -o /out/${sample}_counts.txt \
        -t exon \
        --countReadPairs \
        -g Parent \
	--countReadPairs \
        -T $SLURM_CPUS_PER_TASK \
        -p \
        -s 2 \
        /out/${sample}.bam

echo "Feature counting completed for sample $sample"

# ==============================
# WAIT FOR ALL JOBS TO COMPLETE AND CREATE COUNTS MATRIX
# ==============================
# This section will only run for the last array task (task 5 in this case)
if [ $SLURM_ARRAY_TASK_ID -eq 5 ]; then
    echo "Waiting for all alignment and counting jobs to complete..."
    
    # Wait for all BAM files to be created and indexed
    expected_samples=6  # Adjust this to your actual number of samples
    
    while true; do
        bam_count=$(ls ${OUT_DIR}/*.bam 2>/dev/null | wc -l)
        counts_count=$(ls ${OUT_DIR}/*_counts.txt 2>/dev/null | wc -l)
        
        if [ $bam_count -eq $expected_samples ] && [ $counts_count -eq $expected_samples ]; then
            echo "All files ready. Creating combined counts matrix..."
            break
        else
            echo "Waiting... Found $bam_count BAM files and $counts_count count files out of $expected_samples expected"
            sleep 30
        fi
    done
    
    # Create combined counts matrix
    echo "Creating combined counts matrix..."
    apptainer exec \
        --bind $OUT_DIR:/out \
        --bind $(dirname $ANNOTATION):/annotation \
        $CONTAINER \
        featureCounts \
            -a /annotation/$(basename $ANNOTATION) \
            -o /out/combined_counts_matrix.txt \
            -t exon \
            --countReadPairs \
            -g Parent \
            -T $SLURM_CPUS_PER_TASK \
            --countReadPairs \
            -p \
            -s 2 \
            /out/*.bam
    

    echo "Combined counts matrix created: ${OUT_DIR}/combined_counts_matrix.txt"
    
    # Create summary statistics
    echo "Creating alignment summary..."
    apptainer exec \
        --bind $OUT_DIR:/out \
        $CONTAINER \
        bash -c "
            echo 'Sample,Total_Reads,Mapped_Reads,Mapping_Rate' > /out/alignment_summary.csv
            for bam in /out/*.bam; do
                sample=\$(basename \$bam .bam)
                total=\$(samtools view -c \$bam)
                mapped=\$(samtools view -c -F 4 \$bam)
                rate=\$(echo \"scale=2; \$mapped * 100 / \$total\" | bc -l)
                echo \"\$sample,\$total,\$mapped,\$rate%\" >> /out/alignment_summary.csv
            done
        "
    
    echo "Analysis complete!"
    echo "Output files:"
    echo "  - Individual BAM files: ${OUT_DIR}/*.bam"
    echo "  - Individual count files: ${OUT_DIR}/*_counts.txt"
    echo "  - Combined counts matrix: ${OUT_DIR}/combined_counts_matrix.txt"
    echo "  - Alignment summary: ${OUT_DIR}/alignment_summary.csv"
fi

echo "Finished sample $sample"
