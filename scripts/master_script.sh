#!/bin/bash
#SBATCH --array=1-86
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --time=48:48:48
#SBATCH --mem-per-cpu=6G

# This script performs read trimmining, alignment and outputs filtered variants

# Create directories if they don't exist
make_dirs() {
    for dir in trimmed_seq alignments duplicate_metrics marked_duplicates raw_var raw_snp filtered_snp bqsr_snps recal_data recal_reads raw_variants_recal raw_snps_recal; do
        [ -e "$dir" ] || mkdir "$dir"
    done
}

# Trim reads
trim_reads() {
    accession=$1
    cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG \
                -q 25 \
                -o trimmed_seq/"$accession"_1.fastq.gz \
                -p trimmed_seq/"$accession"_2.fastq.gz \
                "$accession"_1.fastq.gz \
                "$accession"_2.fastq.gz
}

# Align reads
align_reads() {
    accession=$1
    bwa mem -k 19 -r 2.5 \
                    -R "@RG\tID:$accession\tSM:$accession\tPL:illumina\tLB:$accession\tPU:1" \
                    -o alignments/"$accession".sam \
                    reference/reference.fasta \
                    trimmed_seq/"$accession"_1.fastq.gz trimmed_seq/"$accession"_2.fastq.gz
}

# Sort BAM file
sort_bam() {
    accession=$1
    samtools sort alignments/"$accession".sam -o alignments/"$accession".bam
}

# Mark duplicates
mark_duplicates() {
    accession=$1
    picard MarkDuplicates \
                    -I alignments/"$accession".bam \
                    -O marked_duplicates/"$accession".bam \
                    -M duplicate_metrics/marked_dup_metrics_"$accession".txt
}

# Call variants
call_variants() {
    accession=$1
    gatk HaplotypeCaller \
        --input marked_duplicates/"$accession".bam \
        --output raw_var/"$accession".vcf \
        --reference reference/reference.fasta \
        -ploidy 1 \
        --standard-min-confidence-threshold-for-calling 30
}

# Filter so only SNPs are selected for further analysis
filter_snps() {
    accession=$1
    gatk SelectVariants \
        -R reference/reference.fasta \
        -V raw_var/"$accession".vcf \
        -select-type SNP \
        -O raw_snp/"$accession".vcf
}

# Perform variant filtration
variant_filtration() {
    accession=$1
    gatk VariantFiltration \
        -R reference/reference.fasta \
        -V raw_snp/"$accession".vcf \
        -O filtered_snp/"$accession".vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 55.0" \
        -filter-name "MQ_filter" -filter "MQ < 50.0" \
        -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -filter-name "MQRankSum_filter" -filter "MQRankSum < -8.0" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"
}

# Exclude filtered variants
exclude_filtered_variants() {
    accession=$1
    gatk SelectVariants \
        --exclude-filtered \
        -V filtered_snp/"$accession".vcf \
        -O bqsr_snps/"$accession".vcf
}

# Perform base quality score recalibration
bqsr() {
    accession=$1
    gatk BaseRecalibrator \
        -R reference/reference.fasta \
        -I marked_duplicates/"$accession".bam \
        --known-sites bqsr_snps/"$accession".vcf \
        -O recal_data/"$accession".table
}

# Apply base quality score recalibration
apply_bqsr() {
    accession=$1
    gatk ApplyBQSR \
        -R reference/reference.fasta \
        -I marked_duplicates/"$accession".bam \
        -bqsr recal_data/"$accession".table \
        -O recal_reads/"$accession".bam
}

# Call variants after BQSR
call_variants_after_bqsr() {
    accession=$1
    gatk HaplotypeCaller \
        -R reference/reference.fasta \
        -I recal_reads/"$accession".bam \
        -O raw_variants_recal/"$accession".vcf
}

# Select only SNPs after BQSR
filter_variants_after_bqsr() {
    accession=$1
    gatk SelectVariants \
        -R reference/reference.fasta \
        -V raw_variants_recal/"$accession".vcf \
        -select-type SNP \
        -O raw_snps_recal/"$accession".vcf
}

# Read the file names from the text file into an array
mapfile -t files < accessions.txt

# Calculate the start and end index for the current job. In this example
# 3 sequences are processed per job.
start=$((SLURM_ARRAY_TASK_ID * 3 - 2))
end=$((SLURM_ARRAY_TASK_ID * 3))

# Process files in the current partition
for ((i = start - 1; i < end; i++)); do
    file="${files[i]}"
    accession=$(basename "${file%.*}")
    
    # Create necessary directories if they don't exist
    create_directories
    
    # Trim reads
    trim_reads "$accession"
    if [ $? -eq 0 ]; then
        # Delete input FASTQ files
        rm "$accession"_1.fastq.gz "$accession"_2.fastq.gz
        # Align reads
        align_reads "$accession"
        if [ $? -eq 0 ]; then
            # Sort BAM file
            sort_bam "$accession"
            if [ $? -eq 0 ]; then
                # Remove SAM file
                rm alignments/"$accession".sam
                # Mark duplicates
                mark_duplicates "$accession"
                if [ $? -eq 0 ]; then
                    # Index BAM file
                    samtools index marked_duplicates/"${accession}.bam"
                    if [ $? -eq 0 ]; then
                        # Call variants
                        call_variants "$accession"
                        if [ $? -eq 0 ]; then
                            # Filter SNPs
                            filter_snps "$accession"
                            if [ $? -eq 0 ]; then
                                # Perform variant filtration
                                variant_filtration "$accession"
                                if [ $? -eq 0 ]; then
                                    # Exclude filtered variants
                                    exclude_filtered_variants "$accession"
                                    if [ $? -eq 0 ]; then
                                        # Perform BQSR
                                        bqsr "$accession"
                                        if [ $? -eq 0 ]; then
                                            # Apply BQSR
                                            apply_bqsr "$accession"
                                            if [ $? -eq 0 ]; then
                                                # Call variants after BQSR
                                                call_variants_after_bqsr "$accession"
                                                if [ $? -eq 0 ]; then
                                                    # Filter variants after BQSR
                                                    filter_variants_after_bqsr "$accession"
                                                fi
                                            fi
                                        fi
                                    fi
                                fi
                            fi
                        fi
                    fi
                fi
            fi
        fi
    fi
done
