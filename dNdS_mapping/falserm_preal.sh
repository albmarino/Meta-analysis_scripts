#!/bin/bash

# script to remove alignments with more than n% of missing information before alignment.
# Required: *prealign_species.fasta (output of idrename_splitpreal.sh) containing sequences from species of the corresponding alignment clade, and <clade>_<genename>_final_mask_align_NT.aln already computed 

# Provide a newline-separated list of pre-alignment multifasta files as 1st argument, minimum threshold to keep a sequence (i.e. 0.9 to exclude sequences with more than 10% of missing info), clade name
# example usage: bash falserm_preal.sh genelist "0.9" actinopteri

for preal_fasta in $(cat $1); do

	gene=$(preal_fasta%_prealign_species.fasta)
	java -jar ~/bin/macse_v2.06.jar -prog trimNonHomologousFragments -seq $preal_fasta -min_internal_homology_to_keep_seq 0.9 -out_trim_info "$clade"_"$gene"_stats.csv
	done
	ls *stats.csv > statlist
	ls $1*aln > alnlist 

	for statfile in $(cat statlist); do
	statgene=${statfile%_stats.csv}
	if grep -q $statgene alnlist; then
	grep true $statfile | cut -d';' -f1 > sp2keep.txt
	seqkit grep -f sp2keep.txt "$statgene"_final_mask_align_NT.aln -o falserm_"$statgene"_final_mask_align_NT.aln
	fi
	done

done

# outputs falserm_<clade>_<gene>_final_mask_align_NT.aln to be used for dN/dS analyses