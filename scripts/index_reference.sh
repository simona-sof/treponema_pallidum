#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=5:5:5
#SBATCH --mem-per-cpu=2G

# Function to create sequence dictionary
create_sequence_dictionary() {
    picard CreateSequenceDictionary R=reference.fasta O=reference.dict
}

# Function to create FASTA index
create_fasta_index() {
    samtools faidx reference.fasta
}

# Function to create BWA index
create_bwa_index() {
    bwa index reference.fasta
}

# Call the functions
create_sequence_dictionary
create_fasta_index
create_bwa_index
