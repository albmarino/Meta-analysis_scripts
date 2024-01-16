#!/bin/bash

# Script to count the number of genes per species over overall and clade total, excluding the species in a list

# Required input: all the falserm_<clade>_<genename>_final_mask_NT.aln files and the sp2rm list output by gc3plots.R in the pwd
# Example run: bash genespersp_filtergenespersp.sh sp2rm_50percgenespersp
# Output: genespersp_50percgenespersp.tsv, can be used again as input with gc3plots.R

pref=${1#sp2rm_}
cat falserm_mollusca*aln|grep ">" | sort -u | sed 's/>//g' | grep -Fv -f $1 > spmol
cat falserm_mammalia*aln|grep ">" | sort -u | sed 's/>//g' | grep -Fv -f $1 > spmam
cat falserm_insecta*aln|grep ">" | sort -u | sed 's/>//g' | grep -Fv -f $1 > spins
cat falserm_aves*aln|grep ">" | sort -u | sed 's/>//g' | grep -Fv -f $1 > spave
cat falserm_actinopteri*aln|grep ">" | sort -u | sed 's/>//g' | grep -Fv -f $1 > spact
cat spmol spmam spins spave spact > sptot

echo -e "Clade\tSpecies\tngenes_sp\tngenes_clade\tngenes_total" > genespersp_"${pref}".tsv

ngenes_mol=$(ls falserm_mollusca*aln| wc -l)
ngenes_mam=$(ls falserm_mammalia*aln| wc -l)
ngenes_ins=$(ls falserm_insecta*aln| wc -l)
ngenes_ave=$(ls falserm_aves*aln| wc -l)
ngenes_act=$(ls falserm_actinopteri*aln| wc -l)
ngenes_tot=$(ls *_final_mask_align_NT.aln | cut -d'_' -f3 | sort -u | wc -l)

# Genes level across total dataset
for sp in $(cat spmol); do

	ngenes_sp=$(grep -E ${sp}$ falserm_mollusca*aln | wc -l)
	echo -e "Mollusca\t$sp\t$ngenes_sp\t$ngenes_mol\t$ngenes_tot" >> genespersp_"${pref}".tsv
done

for sp in $(cat spmam); do

        ngenes_sp=$(grep -E ${sp}$ falserm_mammalia*aln | wc -l)
        echo -e "Mammalia\t$sp\t$ngenes_sp\t$ngenes_mam\t$ngenes_tot" >> genespersp_"${pref}".tsv
done

for sp in $(cat spins); do

        ngenes_sp=$(grep -E ${sp}$ falserm_insecta*aln | wc -l)
        echo -e "Insecta\t$sp\t$ngenes_sp\t$ngenes_ins\t$ngenes_tot" >> genespersp_"${pref}".tsv
done

for sp in $(cat spave); do

        ngenes_sp=$(grep -E ${sp}$ falserm_aves*aln | wc -l)
        echo -e "Aves\t$sp\t$ngenes_sp\t$ngenes_ave\t$ngenes_tot" >> genespersp_"${pref}".tsv
done

for sp in $(cat spact); do

        ngenes_sp=$(grep -E ${sp}$ falserm_actinopteri*aln | wc -l)
        echo -e "Actinopteri\t$sp\t$ngenes_sp\t$ngenes_act\t$ngenes_tot" >> genespersp_"${pref}".tsv
done

rm spmol spmam spins spave spact sptot
