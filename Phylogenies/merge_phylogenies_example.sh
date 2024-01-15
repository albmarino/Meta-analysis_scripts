#!/bin/bash

# Example script to merge two clade phylogenies (Mammalia and Aves)

# Required input: newline-separated suffixes of genes used for the phylogeny (shared 107_sharedmarkers_cladephylogeny, from Table S6), *final_align_NT.fasta files (output from run_omm_macse.sh) of the clades to merge
# Example usage: bash merge_phylogenies_example.sh 107_sharedmarkers_cladephylogeny

# Outputs mammalian and avian phylogenies rooted with avian and mammalian outgroup. After manual merging of the trees, the same procedure can be further applied to the other taxa until the full phylogeny is reached. 

for gene in $(cat $1); do

# align platypus sequence to the avian alignment
seqkit grep –n -p mammalia_"${gene}"_final_align_NT.fasta -s Ornithorhynchus_anatinus | sed "s/-//g"  > Ornithorhynchus_anatinus.fasta # extract gapless outgroup gene
java -jar ~/bin/macse_v2.06.jar -prog enrichAlignment -align aves_"${gene}"_final_mask_align_NT.aln -seq Ornithorhynchus_anatinus.fasta # align it to the avian alignment; output is aves_${gene}_final_mask_align_NT_AA.aln
sed -i 's/!/-/g' aves_"${gene}"_final_mask_align_NT_AA.aln # Replace ! with - (symbol for frameshift mutations in MACSE)
sed -i 's/!/-/g' aves_"${gene}"_final_mask_align_NT.aln

# align kiwi sequence to the mammalian alignment
seqkit grep –n -p aves_"${gene}"_final_align_NT.fasta -s Apteryx_rowi | sed "s/-//g"  > Apteryx_rowi.fasta
java -jar ~/bin/macse_v2.06.jar -prog enrichAlignment -align mammalia_"${gene}"_final_mask_align_NT.aln -seq Apteryx_rowi.fasta
sed -i 's/!/-/g' mammalia_"${gene}"_final_mask_align_NT_AA.aln
sed -i 's/!/-/g' mammalia_"${gene}"_final_mask_align_NT.aln

done

# concatenate
grep '>' mammalia*_final_mask_align_NT_AA.aln | cut -d: -f2 | sort -u | sed 's/>//' > list_esp1.txt # mammalia + outgroup species list
ls mammalia*_final_mask_align_NT_AA.aln > list_file1.txt # list of mammalia + outgroup alignment files
Concatenate list_file1.txt mammalia_outgroup.fst list_esp1.txt AA Fasta

grep '>' aves*_final_mask_align_NT_AA.aln | cut -d: -f2 | sort -u | sed 's/>//' > list_esp2.txt # aves + outgroup species list
ls $aves*_NT_AA.fasta > list_file2.txt # list of aves + outgroup alignment files
Concatenate list_file2.txt aves_outgroup.fst list_esp2.txt AA Fasta

# root trees with their outgroup
Rscript excludeSpeciesTree_edit.R mammalia_outgroup.fst mammalia_out_100.fst.treefile # outputs pruned new_mammalia_out_100.fst.treefile
iqtree -s mammalia_outgroup.fst -m LG+G -g new_mammalia_out_100.fst.treefile -pre mammalia_outgroup # use previously calculated mammalian topology to constrain the rooted tree
Rscript excludeSpeciesTree_edit.R aves_outgroup.fst aves_out_107.fst.treefile
iqtree -s aves_outgroup.fst -m LG+G -g new_aves_out_107.fst.treefile -pre aves_outgroup