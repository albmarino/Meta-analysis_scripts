#!/bin/bash

# Script to concatenate aminoacid gene sequences of a clade
# Required input: AA alignment files output by run_omm_macse.sh selected for the phylogeny

# Example usage: bash concatenate.sh actinopteri

grep '>' "${1}"*align_AA.aln | cut -d: -f2 | sort -u | sed 's/>//' > list_esp.txt # get the species list

ls "${1}"*_final_align_AA.aln > list_file.txt # list of alignment files

./Concatenate list_file.txt "${1}"_out_107.fst list_esp.txt AA Fasta