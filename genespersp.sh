#!/bin/bash

# bash genespersp.sh

cat falserm_mollusca*aln|grep ">" | sort -u | sed 's/>//g' > spmol
cat falserm_mammalia*aln|grep ">" | sort -u | sed 's/>//g' > spmam
cat falserm_insecta*aln|grep ">" | sort -u | sed 's/>//g' > spins
cat falserm_aves*aln|grep ">" | sort -u | sed 's/>//g' > spave
cat falserm_actinopteri*aln|grep ">" | sort -u | sed 's/>//g' > spact
cat spmol spmam spins spave spact > sptot

echo -e "Clade\tSpecies\tngenes_sp\tngenes_clade\tngenes_total" > genespersp.tsv

ngenes_mol=$(ls falserm_mollusca*aln| wc -l)
ngenes_mam=$(ls falserm_mammalia*aln| wc -l)
ngenes_ins=$(ls falserm_insecta*aln| wc -l)
ngenes_ave=$(ls falserm_aves*aln| wc -l)
ngenes_act=$(ls falserm_actinopteri*aln| wc -l)
ngenes_tot=$(ls *_final_mask_align_NT.aln | cut -d'_' -f3 | sort -u | wc -l)

# Genes level across total dataset
for sp in $(cat spmol); do

	ngenes_sp=$(grep -E ${sp}$ falserm_mollusca*aln | wc -l)
	echo -e "Mollusca\t$sp\t$ngenes_sp\t$ngenes_mol\t$ngenes_tot" >> genespersp.tsv
done

for sp in $(cat spmam); do

        ngenes_sp=$(grep -E ${sp}$ falserm_mammalia*aln | wc -l)
        echo -e "Mammalia\t$sp\t$ngenes_sp\t$ngenes_mam\t$ngenes_tot" >> genespersp.tsv
done

for sp in $(cat spins); do

        ngenes_sp=$(grep -E ${sp}$ falserm_insecta*aln | wc -l)
        echo -e "Insecta\t$sp\t$ngenes_sp\t$ngenes_ins\t$ngenes_tot" >> genespersp.tsv
done

for sp in $(cat spave); do

        ngenes_sp=$(grep -E ${sp}$ falserm_aves*aln | wc -l)
        echo -e "Aves\t$sp\t$ngenes_sp\t$ngenes_ave\t$ngenes_tot" >> genespersp.tsv
done

for sp in $(cat spact); do

        ngenes_sp=$(grep -E ${sp}$ falserm_actinopteri*aln | wc -l)
        echo -e "Actinopteri\t$sp\t$ngenes_sp\t$ngenes_act\t$ngenes_tot" >> genespersp.tsv
done

rm spmol spmam spins spave spact sptot
