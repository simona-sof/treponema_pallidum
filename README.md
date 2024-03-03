# Treponema pallidum Nextstrain build

This is a Nextstrain build for Treponema pallidum accompanied by the data of two subspieces of Treponema Pllidum: TPE and TPA. TPA subspieces is then devided into two lineages: SS14 and Nichols.

Scripts folder contains scripts used to process reads. Read processing in detail is described in the project report. In short, the reference was stored in a separate folder, and indexed. The sequence files were downladed with down.sh script. The master_script.sh contains the commands used to process the reads and get a recalibrated vcf file containing recalibrated SNPs. The consensus sequences then can be found with consensus.sh, the sequences then can be merged (merge.sh) and renamed for Nextstrain analysis (rename.py).

The Snakefile contains Augur commands that can be run on Nextstrain container to generate Auspice files for phylogeny viewing. It is a general command file therefore, the paths should be adjusted depending on which lineage is being analysed. The metadata for each phylogeny, the combined fasta file of all sequences, the raw variants and SNPs after recalibration are each in separate folder for the respective subspieces/lineage. The configurations folder contains the Auspice configurations. Additional bash script, augur.sh, was created to run the Augur commands without the Nextstrain container.

The metadata was cleaned up, the python scripts in pre-processing folder contain the meta pre-processing, clean_meta.py contains script to find longitutes and latitudes, adjust the date format, so it is usable on Nextstrain.

Comparison folder contains the SS14-lineage phylogeny JSON files for 13 sequences processed in this project and phylogeny from Majander et al 2023. The sequences used for the comparison are not masked for recombinant sites.

