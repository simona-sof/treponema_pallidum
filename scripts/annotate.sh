#!/bin/bash

#SBATCH --ntasks=6
#SBATCH --time=5:5:5
#SBATCH --mem-per-cpu=4G

for file in `cat raw_trimmed/accessions.txt`; do
	java -jar snpEff/snpEff.jar -v Treponema_pallidum_subsp_pallidum_ss14 snp/${file}.vcf.gz > annotated/${file}.vcf.gz
done
