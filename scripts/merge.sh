#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:48:48
#SBATCH --mem-per-cpu=2G

# This script merges the consensus sequences into one file 

mapfile -t files < accessions.txt
for accession in "${files[@]}"; do
        cat fastas/"${accession}.fasta" >> "merged.fasta"
done