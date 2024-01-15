#!/bin/bash

# align genes with omm_macse from a provided newline-separated list of multifasta files
# Required: *prealign_species.fasta (output of idrename_splitpreal.sh) containing sequences from species of the corresponding alignment clade
# example usage: bash run_omm_macse.sh genelist

for preal_fasta in $(cat $1); do

	singularity run ~/bin/omm_macse_v11.05b.sif --out_dir ./"${preal_fasta}"_aligned --out_file_prefix ${preal_fasta} --in_seq_file $preal_fasta --genetic_code_number 1 --java_mem 20000m

done
