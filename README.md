# Treponema pallidum Nextstrain build

[comment]: <> (This is a Nextstrain build for Treponema pallidum accompanied by the data of two subspieces of Treponema Pllidum: TPE and TPA. TPA subspieces is then devided into two lineages: SS14 and Nichols.)

[comment]: <> (Scripts folder contains scripts used to process reads. Read processing in detail is described in the project report. In short, the reference was stored in a separate folder, and indexed. The sequence files were downladed with down.sh script. The master_script.sh contains the commands used to process the reads and get a recalibrated vcf file containing recalibrated SNPs. The consensus sequences then can be found with consensus.sh, the sequences then can be merged merge.sh and renamed for Nextstrain analysis rename.py.)

[comment]: <> (The Snakefile contains Augur commands that can be run on Nextstrain container to generate Auspice files for phylogeny viewing. It is a general command file therefore, the paths should be adjusted depending on which lineage is being analysed. The metadata for each phylogeny, the combined fasta file of all sequences, the raw variants and SNPs after recalibration are each in separate folder for the respective subspieces/lineage. The configurations folder contains the Auspice configurations. Additional bash script, augur.sh, was created to run the Augur commands without the Nextstrain container.)

[comment]: <> (The metadata was cleaned up, the python scripts in pre-processing folder contain the meta pre-processing, clean_meta.py contains script to find longitutes and latitudes, adjust the date format, so it is usable on Nextstrain.)

[comment]: <> (Comparison folder contains the SS14-lineage phylogeny JSON files for 13 sequences processed in this project and phylogeny from Majander et al 2023. The sequences used for the comparison are not masked for recombinant sites.)

[comment]: <> (This repository contains scripts and data for building phylogenies of Treponema pallidum subspieces TPE and TPA, divided into lineages SS14 and Nichols. The pipeline utilizes Nextstrain's Augur and Auspice tools for phylogeny visualization.)
## File Structure

- **Scripts:** Contains scripts for read processing and analysis.
- **Pre-processing:** Contains Python scripts for metadata cleanup.
- **Comparison:** Includes phylogeny JSON files for comparison with external data.
- **Configurations:** Contains Auspice configuration files.
- **Snakefile:** Defines Nextstrain pipeline commands.
- **Augur.sh:** Bash script to run Augur commands outside the Nextstrain container.

## Data Processing

- **Read Processing**: 
    - Raw reads in fastq format were processed to obtain mapped sequences in fasta format. 
    - Reference genomes were indexed, adapters were removed, and sequences were mapped using BWA mem. 
    - Mapped genomes were sorted, cleaned, duplicates were marked and removed, and raw variants were called using GATK toolkit.

- **Variant Filtering**: 
    - Variants were filtered based on quality metrics using GATK VariantFiltration, and base recalibration was performed.

- **Consensus Sequence Generation**: 
    - Consensus sequences were generated using bcftools consensus.

## Nextstrain Pipeline

- **Augur Commands**: 
    - The Snakefile contains Augur commands for generating Auspice files. Adjust paths based on the lineage being analyzed.

- **Metadata Cleaning**: 
    - Metadata was cleaned using Python scripts in the pre-processing folder. Adjust date formats and extract geographical coordinates for Nextstrain compatibility.

- **Tree Building**: 
    - Different trees were built for different lineages of T. pallidum using IQtree. TreeTime was used for refinement, ancestral sequence inference, and trait assignment (e.g., location).

- **Limitations**: 
    - Technical difficulties limited the input of a variant file containing shared SNPs, restricting Nextstrain's capabilities. Clades only present nucleotide mutations due to this limitation.

## Usage

1. Ensure dependencies (BWA, Samtools, Cutadapt, Picard toolkit, GATK toolkit, bcftools) are installed.
2. Execute scripts in the provided order for read processing and variant calling: download raw sequences with down.sh, index reference with index.sh, master_script.sh contains the commands used to process the reads and get a recalibrated vcf file containing recalibrated SNPs, call consensus sequences with consensus.sh, merge with merge.sh and rename with rename.py. The final fasta containing all sequences can be used as inpu in Augur pipeline.
4. Run the Snakefile using Nextstrain container to generate JSON files. Don't forget to change the filepaths!
5. View phylogenies from the generated JSON file with Auspice.

## References
- Li, H. (2013). Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997.
- Danecek, P. et al. (2021). Twelve years of SAMtools and BCFtools. Genome Res. 31(1): 23-30.
- Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal. 17(1): 10-12.
- O’Connor, B., and van der Auwera, G. (2020). Genomics in the cloud: Using docker, GATK, and WDL in terra. O’Reilly Media, Incorporated. tex.lccn: 2021275634
- Picard Tools - By Broad Institute. [Online]. Available: https://broadinstitute.github.io/picard/
- Hadfield, J.  et al. (2018). Nextstrain: real-time tracking of pathogen evolution. Bioinformatics. 34(23): 4121–4123.
