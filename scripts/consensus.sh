#!/bin/bash
#SBATCH --array=1-86
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:48:48
#SBATCH --mem-per-cpu=2G

#This script outputs consensus sequence from the variants 

# Read the file names from the text file into an array
mapfile -t files < accessions.txt

# Create directories if they don't exist
make_dirs() {
    for dir in fastas zipped; do
        [ -e "$dir" ] || mkdir "$dir"
    done
}

# Calculate the start and end index for the current job. In this example
# 3 sequences are processed per job.
start=$((SLURM_ARRAY_TASK_ID * 3 - 2))
end=$((SLURM_ARRAY_TASK_ID * 3))

# Process files in the current partition
for ((i = start - 1; i < end; i++)); do
    file="${files[i]}"
    accession=$(basename "${file%.*}")
     # Index the VCF file using bcftools
    bcftools index -f raw_snp_recal/"$accession".vcf.gz

    # Generate a consensus FASTA sequence based on the reference genome
    # and the variants found in the VCF file using bcftools
    cat reference/reference.fasta | bcftools consensus raw_snp_recal/"$accession".vcf.gz > fastas/"$accession".fasta
done