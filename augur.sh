#!/bin/bash

# Names of files used in the analysis
seq_file="sequence.fasta"
ref_file="reference/reference.fasta"
meta_file="meta.csv"
exclude_file="config/dropped_strains_TPE.txt"
auspice_config_file="config/auspice_config.json"
colors_file="config/color.tsv"
geo_info_file="config/lat_longs.tsv"

filter() {
    augur filter --sequences "$seq_file" \
        --metadata "$meta_file" \
        --exclude "$exclude_file" \
        --output "results/filtered.fasta"
}

mask() {
    augur mask --sequences "${1}" \
        --mask "$mask_file" \
        --output "results/masked.fasta"
}

tree() {
    augur tree --alignment "${1}" \
        --vcf-reference "$ref_file" \
        --method 'iqtree' \
        --output "results/tree_raw.nwk"
}



refine() {
    augur refine --tree "${1}" \
        --alignment "${2}" \
        --vcf-reference "$ref_file" \
        --metadata "$meta_file" \
        --timetree \
        --root 'min_dev' \
        --coalescent 'opt' \
        --output-tree "results/tree.nwk" \
        --output-node-data "results/branch_lengths.json"
}

ancestral() {
    augur ancestral --tree "${1}" \
        --alignment "${2}" \
        --inference 'joint' \
        --output-node-data "results/nt_muts.json" \
        --output-sequences "results/ancestral.fasta"
}


traits() {
    augur traits --tree "${1}" \
        --metadata "$meta_file" \
        --columns 'location' \
        --output-node-data "results/traits.json"
}

export_data() {
    augur export v2 \
        --tree "${1}" \
        --metadata "$meta_file" \
        --node-data "${2}" "${3}" "${4}" \
        --auspice-config "$auspice_config_file" \
        --colors "$colors_file" \
        --lat-longs "$geo_info_file" \
        --output "auspice/tp.json"
}


# Execute rules. Don't forget to change the directories as needed!
filter
#mask "results/filtered_.fasta"
tree "results/filtered.fasta"
refine "results/tree_raw.nwk" "results/filtered.fasta"
ancestral "results/tree.nwk" "results/filtered.fasta"
clades "results/tree.nwk" "results/nt_muts.json"
traits "results/tree.nwk"
export_data "results/tree.nwk" "results/branch_lengths.json" "results/traits.json" "results/nt_muts.json"

