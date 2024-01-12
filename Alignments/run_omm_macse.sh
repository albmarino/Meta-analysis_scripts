#!/bin/bash

# align genes with omm_macse from a provided newline-separated list of multifasta files
# example usage: bash run_omm_macse.sh genelist

for preal_fasta in $(cat $1); do

	singularity run ~/bin/omm_macse_v11.05b.sif --out_dir ./"${gene}"_aligned --out_file_prefix ${gene} --in_seq_file $preal_fasta --genetic_code_number 1 --java_mem 20000m

done