#!/bin/bash

# Script to concatenate GC3-selected genes for Coevol

# Required input: lists of geneset and of species to exclude output from gc3plots.R {or gc3plots_gcrich.R} 
# Example usage: bash concatenate_4coevol.sh gcpoorset_tot_50percgenespersp_95_20_20.tsv sp2rm_50percgenespersp

settype=$(echo $1 | cut -d'_' -f1) # gcpoorset or gcrichset
dataset=$(echo $1 | cut -d'_' -f2) # tot or clade
nmarkers=$(echo $1 | cut -d'_' -f5) # 20 or 50

if [[ "${dataset}" == "tot" ]] ; then # if we are dealing with total dataset...
	echo -e "$nmarkers genes from $settype for $dataset dataset are going to be aligned: this might take some time...\n"
	echo "Removing species with more than 50% of genes missing"
	for gene in $(cat $1); do

		for clade in $(echo -e "actinopteri\naves\ninsecta\nmammalia\nmollusca"); do
			if [ -e falserm_${clade}_${gene}_final_mask_align_NT.aln ]; then
				seqkit grep -n -v -f sp2rm_50percgenespersp falserm_${clade}_${gene}_final_mask_align_NT.aln > spremoved_falserm_${clade}_${gene}_final_mask_align_NT.aln
			else
				echo "$clade_$gene does not exist, skipping it..."
			fi
		done
		# align profiles Mollusca + Insecta -> Protostomata
		echo -e "Starting alignment for $gene\n"
		if [ -e spremoved_falserm_mollusca_${gene}_final_mask_align_NT.aln ] && [ -e spremoved_falserm_insecta_${gene}_final_mask_align_NT.aln ]; then
			java -jar ~/bin/macse_v2.06.jar -prog alignTwoProfiles -p1 spremoved_falserm_insecta_${gene}_final_mask_align_NT.aln -p2 spremoved_falserm_mollusca_${gene}_final_mask_align_NT.aln -out_NT protostomata_${gene}_NT.fasta -out_AA protostomata_${gene}_AA.fasta
		elif [ -e spremoved_falserm_mollusca_${gene}_final_mask_align_NT.aln ] && [ ! -e spremoved_falserm_insecta_${gene}_final_mask_align_NT.aln ]; then
			cp spremoved_falserm_mollusca_${gene}_final_mask_align_NT.aln protostomata_${gene}_NT.fasta
		elif [ ! -e spremoved_falserm_mollusca_${gene}_final_mask_align_NT.aln ] && [ -e spremoved_falserm_insecta_${gene}_final_mask_align_NT.aln ]; then
			cp spremoved_falserm_insecta_${gene}_final_mask_align_NT.aln protostomata_${gene}_NT.fasta
		fi
		sed -i 's/!/-/g' protostomata_${gene}_NT.fasta

		# align profiles Mammalia + Aves -> Tetrapoda
		if [ -e spremoved_falserm_mammalia_${gene}_final_mask_align_NT.aln ] && [ -e spremoved_falserm_aves_${gene}_final_mask_align_NT.aln ]; then
			java -jar ~/bin/macse_v2.06.jar -prog alignTwoProfiles -p1 spremoved_falserm_mammalia_${gene}_final_mask_align_NT.aln -p2 spremoved_falserm_aves_${gene}_final_mask_align_NT.aln -out_NT tetrapoda_${gene}_NT.fasta -out_AA tetrapoda_${gene}_AA.fasta
		elif [ -e spremoved_falserm_mammalia_${gene}_final_mask_align_NT.aln ] && [ ! -e spremoved_falserm_aves_${gene}_final_mask_align_NT.aln ]; then
			cp spremoved_falserm_mammalia_${gene}_final_mask_align_NT.aln tetrapoda_${gene}_NT.fasta
		elif [ ! -e spremoved_falserm_mammalia_${gene}_final_mask_align_NT.aln ] && [ -e spremoved_falserm_aves_${gene}_final_mask_align_NT.aln ]; then
			cp spremoved_falserm_aves_${gene}_final_mask_align_NT.aln tetrapoda_${gene}_NT.fasta
		fi
		sed -i 's/!/-/g' tetrapoda_${gene}_NT.fasta

		# align profiles Tetrapoda+ Actinopteri -> Vertebrata
		if [ -e spremoved_falserm_actinopteri_${gene}_final_mask_align_NT.aln ]; then
			java -jar ~/bin/macse_v2.06.jar -prog alignTwoProfiles -p1 spremoved_falserm_actinopteri_${gene}_final_mask_align_NT.aln -p2 tetrapoda_${gene}_NT.fasta -out_NT vertebrata_${gene}_NT.fasta -out_AA vertebrata_${gene}_AA.fasta
		elif [ ! -e spremoved_falserm_actinopteri_${gene}_final_mask_align_NT.aln ]; then
			cp tetrapoda_${gene}_NT.fasta vertebrata_${gene}_NT.fasta
		fi
		sed -i 's/!/-/g' vertebrata_${gene}_NT.fasta

		# align profiles Vertebrata + Protostomata -> Total
		java -jar ~/bin/macse_v2.06.jar -prog alignTwoProfiles -p1 protostomata_${gene}_NT.fasta -p2 vertebrata_${gene}_NT.fasta -out_NT ${dataset}_${gene}_NT.fasta -out_AA ${dataset}_${gene}_AA.fasta
		sed -i 's/!/-/g' ${dataset}_${gene}_NT.fasta
	done

	grep ">" ${dataset}_*_NT.fasta | cut -d: -f2 | sort -u | sed 's/>//' > list_esp.txt
        ls ${dataset}_*_NT.fasta > list_file.txt

	${HOME}/softwares/Concatenate list_file.txt "concatenate_${settype}_${dataset}.fasta" list_esp.txt nuc Fasta

	rm tetrapoda_*fasta vertebrata_*fasta protostomata_*fasta

else # otherwise if we are dealing with a clade dataset...
	echo -e "$nmarkers genes from $settype for $dataset dataset are going to be concatenated\n"
        echo "Removing species with more than 50% of genes missing"
	for gene in $(cat $1); do
		seqkit grep -n -v -f sp2rm_50percgenespersp falserm_${dataset}_${gene}_final_mask_align_NT.aln > spremoved_falserm_${dataset}_${gene}_final_mask_align_NT.aln
	done

	echo "Concatenating..."
	grep ">" spremoved_falserm_${dataset}*final_mask_align_NT.aln | cut -d: -f2 | sort -u | sed 's/>//' > list_esp.txt
	ls spremoved_falserm_${dataset}*final_mask_align_NT.aln > list_file.txt

	~/softwares/Concatenate list_file.txt concatenate_${settype}_${dataset}.fasta list_esp.txt nuc Fasta
fi

rm spremoved* list_file.txt list_esp.txt
