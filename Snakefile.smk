# Snakefile

# Names of files used in the analysis
seq_file = "sequence.fasta"
ref_file = "reference/reference.fasta"
meta_file = "meta.csv"
exclude_file = "config/dropped_strains.txt"
auspice_config_file = "config/auspice_config.json"
colors_file = "config/color.tsv"
geo_info_file = "config/lat_longs.tsv"
# mask_file = "config/mask.bed"

rule filter:
    input:
        seq_file = seq_file,
        meta_file = meta_file,
        exclude_file = exclude_file
    output:
        filtered = "results/filtered.fasta"
    shell:
        '''
        augur filter \
        --sequences {input.seq_file} \
        --metadata {input.meta_file} \
        --exclude {input.exclude_file} \
        --output {output.filtered}
        '''

# rule mask:
#     input:
#         seq_file = "results/filtered.fasta"
#     output:
#         masked = "results/masked.fasta"
#     shell:
#         '''
#         augur mask 
#         --sequences {input.seq_file} 
#         --mask {mask_file} 
#         --output {output.masked}
#         '''

rule tree:
    input:
        alignment = rules.filter.output.filtered,
        ref = ref_file
    output:
        tree_raw = "results/tree_raw.nwk"
    params:
        method = 'iqtree'
    shell:
        '''
        augur tree 
        --alignment {input.alignment} \
        --vcf-reference {ref_file} \
        --method {params.method} \
        --output {output.tree_raw}
        '''

rule refine:
    input:
        tree = rules.tree.output.tree_raw,
        alignment = rules.filter.output.filtered,
        metadata = meta_file,
        ref = ref_file
    output:
        tree = "results/tree.nwk",
        node_data = "results/branch_lengths.json"
    params:
        root = 'min_dev',
        coal = 'opt'
    shell:
        '''
        augur refine \
        --tree {input.tree} \
        --alignment {input.alignment} \
        --vcf-reference {ref_file} \
        --metadata {meta_file} \
        --timetree \
        --root {params.root} \
        --coalescent {params.coal} \
        --output-tree {output.tree} \
        --output-node-data {output.node_data}
        '''

rule ancestral:
    input:
        tree = rules.refine.output.tree,
        alignment = rules.filter.output.filtered,
        ref = ref_file
    output:
        nt_muts = "results/nt_muts.json",
        ancestral_seq = "results/ancestral.fasta"
    params:
        inference = "joint"
    shell:
        '''
        augur ancestral \
        --tree {input.tree} \
        --alignment {input.alignment} \
        --inference {params.inference} \
        --output-node-data {output.nt_muts} \
        --output-sequences {output.ancestral_seq}
        '''

rule traits:
    input:
        tree = rules.refine.output.tree,
    output:
        traits_data = "results/traits.json"
    params:
        traits = 'location'
    shell:
        '''
        augur traits \
        --tree {input.tree} \
        --metadata {meta_file} \
        --columns {params.traits} \
        --output-node-data {output.traits_data}'''

rule export_data:
    input:
        tree = rules.refine.output.tree,
        metadata = meta_file,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.traits_data,
        nt_muts = rules.ancestral.output.nt_data,
        color_defs = colors_file,
        config = "config/auspice_config.json",
        geo_info = geo_info_file,
    output:
        auspice_json = "auspice/tp.json"
    shell:
        '''
        augur export v2 \
        --tree {input.tree} \
        --metadata {meta_file} \
        --node-data {input.branch_lengths} {input.traits} {input.nt_muts} \
        --auspice-config {input.config} \
        --colors {input.color_defs} \
        --lat-longs {input.geo_info} \
        --output {output.auspice_json} \
        """