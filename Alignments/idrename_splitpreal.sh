#!/bin/bash

# renames *singlecopy.fasta file headers and outputs two same multifasta files for each gene: one with a long ID (assembly accession, genename, contig/scaffold name, starting and ending position) (genename_prealign.fasta), and one with just species name (genename_prealign_species.fasta).
# run in the directory with all *singlecopy.fasta files obtained with faatofa_singlecopy_busco_editgenecol.py and acc_species mapfile with tab-separated assembly accessions and species names 

# example run: bash idrename_splitpreal.sh acc_species

## rename *singlecopy.fasta files headers to shorter IDs

for fasta in *.singlecopy.fasta; do

acc=${fasta%.fa.singlecopy.fasta}
grep ">" ${fasta} | sed 's/>//g' > oldid_list # newline-separated long IDs
grep -Po ">[a-zA-Z0-9]+(?=_)(?=|)" ${fasta} | sed 's/>//g' >> all_genenames # make a complete list of all genenames
cat oldid_list | cut -d'|' -f1,2,7,8 | sed "s/^/${acc}|/g"> newid_list # newline-separated new IDs (with assembly accession, genename, contig/scaffold name, starting and ending position)
paste oldid_list newid_list > mapfile

seqkit replace -p '(.+)' -r '{kv}' -k mapfile ${fasta} > ${acc}_singlecopy_renamed.fasta

done

sort -u all_genenames > unique_genenames # filter to a unique list of busco genenames
rm oldid_list newid_list mapfile all_genenames

#########################################################

## For each gene set, create a fasta file containing the corresponding sequence from all the genomes

for gene in $(cat unique_genenames); do # loop through each gene name

echo "" > ${gene}_prealign.fasta

for renamed_fasta in *renamed.fasta; do # retrieve the corresponding sequence from every renamed fasta

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' ${renamed_fasta} > renamed_fasta_oneliner # set the sequence on a single line for grepping
grep -P -A1 "${gene}(?=_)(?=|)" renamed_fasta_oneliner >> ${gene}_prealign.fasta # and append it to the pre-alignment fasta file of this geneset

done

done
rm renamed_fasta_oneliner
rm unique_genenames


#######################################################

## Rename IDs in gene_prealign.fasta with just the species name (gene_prealign_species.fasta)
## acc_species file with old IDs in the first columns and new IDs in the second column (tab-separated) must be already in the pwd

for old_prealign in *_prealign.fasta; do

gene=${old_prealign%.fasta}

for line in ${old_prealign}; do

sed -i 's/|.*//' ${line} # reduce header to just the assembly accession to match the info in the mapfile

done

grep ">" ${old_prealign} | sed 's/>//g' > current_accessions_list # list all the accessions in the present prealign fasta
grep -F -f current_accessions_list $1 > current_mapfile # extract the accession/species couples that we need to map in this file

seqkit replace -p '(.+)' -r '{kv}' -k current_mapfile ${old_prealign} > ${gene}_species.fasta # replace the accesion IDs with species names in <genename>_prealign_species.fasta

done

rm current_accessions_list current_mapfile

#rm *prealign.fasta # keep just the fasta with species ID
