#!/bin/bash

# Example script to re-align aln profiles of two clades (Mammalia and Aves)

# Required input: newline-separated suffixes of genes used for branch lenght recomputation (50_sharedmarkers_brlengths, from Table S6), corresponding *final_align_NT.fasta files (output from run_omm_macse.sh) of the two clades
# Example usage: bash align_profiles_example.sh 50_sharedmarkers_brlengths

# Outputs the alignment of each selected gene across the two clades. The same procedure can be further applied to the other taxa until alignment of the markers for the full phylogeny is reached (i.e., tetrapoda + actinopteri, mollusca + insecta, then vertebrata + protostomata)

for gene in $(cat $1); do

java -jar ~/bin/macse_v2.06.jar -prog alignTwoProfiles -p1 mammalia_"${gene}"_final_mask_align_NT.aln -p2 aves_"${gene}"_final_mask_align_NT.aln -out_NT tetrapoda_"${gene}"_NT.fasta -out_AA tetrapoda_"${gene}"_AA.fasta
sed -i 's/!/-/g' tetrapoda_"${gene}"_AA.fasta

done