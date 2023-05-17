#!/bin/bash
# Extract the GC3 average of every busco gene across all the clades - Actinotperi, Aves, Insecta, Mammalia, Mollusca
# Extract the number of species present in each alignment file - to understand if we have a mostly conserved geneset
# Order them by increasing GC to understand if we have a big enough GC-poor shared geneset
# bash gcstat_filtgenespersp.sh sp2rm_70percgenespersp

thresh=${1#sp2rm_}
sed 's/$/\$/g' $1 > regexdollar_${1}
sed 's/$/,/g' $1 > regexcomma_${1}
#ls *_final_mask_align_NT.aln | cut -d'_' -f3 | sort -u > genelist
echo -e "Gene\tGene_nsp_total\tSp_total\tMissingdata_average_total\tGC3_average_total\tGC3_variance_total\tGene_nsp_actinopteri\tSp_actinopteri\tMissingdata_actinopteri\tGC3_average_actinopteri\tGC3_variance_actinopteri\tGene_nsp_aves\tSp_aves\tMissingdata_aves\tGC3_average_aves\tGC3_variance_aves\tGene_nsp_insecta\tSp_insecta\tMissingdata_insecta\tGC3_average_insecta\tGC3_variance_insecta\tGene_nsp_mammalia\tSp_mammalia\tMissingdata_mammalia\tGC3_average_mammalia\tGC3_variance_mammalia\tGene_nsp_mollusca\tSp_mollusca\tMissingdata_mollusca\tGC3_average_mollusca\tGC3_variance_mollusca" > gc3_report_${thresh}.tsv

suff="_final_mask_align_NT.aln"
act_pref="falserm_actinopteri_"
ave_pref="falserm_aves_"
ins_pref="falserm_insecta_"
mam_pref="falserm_mammalia_"
mol_pref="falserm_mollusca_"

# Count total species overall and per clade
total_sp=$(grep '>' *$suff | grep -vE -f regexdollar_$1 | cut -d'>' -f2| sort -u |wc -l)
actinopteri_sp=$(grep '>' $act_pref*$suff | grep -vE -f regexdollar_$1 | cut -d'>' -f2| sort -u |wc -l)
aves_sp=$(grep '>' $ave_pref*$suff | grep -vE -f regexdollar_$1 |cut -d'>' -f2| sort -u |wc -l)
insecta_sp=$(grep '>' $ins_pref*$suff | grep -vE -f regexdollar_$1 |cut -d'>' -f2| sort -u |wc -l)
mammalia_sp=$(grep '>' $mam_pref*$suff | grep -vE -f regexdollar_$1 | cut -d'>' -f2| sort -u |wc -l)
mollusca_sp=$(grep '>' $mol_pref*$suff | grep -vE -f regexdollar_$1 |cut -d'>' -f2| sort -u |wc -l)



for gene in $(cat genelist); do

# compute average GC3, GC variance and average amount of missing data in each clade
computeSeqStat $act_pref$gene$suff fasta GC | grep -vE -f regexcomma_$1 | cut -d',' -f6 | sed '1d'> gc3temp_actinopteri
computeSeqStat ~/alignedbusco_phylogenies/aligned_busco/extractbusco_actinopteri/${gene}_prealign_species_NT.fasta fasta missing_data | grep -vE -f regexcomma_$1 | cut -d',' -f5 > missingdatatemp_actinopteri

computeSeqStat $ave_pref$gene$suff fasta GC | grep -vE -f regexcomma_$1 | cut -d',' -f6 | sed '1d'> gc3temp_aves
computeSeqStat ~/alignedbusco_phylogenies/aligned_busco/extractbusco_aves/${gene}_prealign_species_NT.fasta fasta missing_data | grep -vE -f regexcomma_$1 | cut -d',' -f5 > missingdatatemp_aves

computeSeqStat $ins_pref$gene$suff fasta GC | grep -vE -f regexcomma_$1 | cut -d',' -f6 | sed '1d'> gc3temp_insecta
computeSeqStat ~/alignedbusco_phylogenies/aligned_busco/extractbusco_insecta/${gene}_prealign_species_NT.fasta fasta missing_data | grep -vE -f regexcomma_$1 | cut -d',' -f5 > missingdatatemp_insecta

computeSeqStat $mam_pref$gene$suff fasta GC | grep -vE -f regexcomma_$1 | cut -d',' -f6 | sed '1d'> gc3temp_mammalia
computeSeqStat ~/alignedbusco_phylogenies/aligned_busco/extractbusco_mammalia/${gene}_prealign_species_NT.fasta fasta missing_data | grep -vE -f regexcomma_$1 | cut -d',' -f5 > missingdatatemp_mammalia

computeSeqStat $mol_pref$gene$suff fasta GC | grep -vE -f regexcomma_$1 | cut -d',' -f6 | sed '1d'> gc3temp_mollusca
computeSeqStat ~/alignedbusco_phylogenies/aligned_busco/extractbusco_mollusca/${gene}_prealign_species_NT.fasta fasta missing_data | grep -vE -f regexcomma_$1 | cut -d',' -f5 > missingdatatemp_mollusca

# Variance formula: var = sommatory i:1->n {(X(i) - X(mean))^2}/n
# Compute GC3 stats and missing data for actinopteri
if [ ! -s gc3temp_actinopteri ]; then
	sumgc_actinopteri=0
	summissingdata_actinopteri=0
	nsp_actinopteri=0
	averagegc_actinopteri="NA"
	averagemissingdata_actinopteri="NA"
	vargc_actinopteri="NA"
else
	sumgc_actinopteri=$(awk '{print $1}' gc3temp_actinopteri | paste -sd+ | bc)
	summissingdata_actinopteri=$(awk '{print $1}' missingdatatemp_actinopteri | paste -sd+ | bc)
	nsp_actinopteri=$(grep -vE -f regexdollar_$1 $act_pref$gene$suff | grep -c '>')
	averagegc_actinopteri=$(echo "$sumgc_actinopteri / $nsp_actinopteri" | bc -l)
	averagemissingdata_actinopteri=$(echo "$summissingdata_actinopteri / $nsp_actinopteri" | bc -l)
	sommsq=0
	for gc3 in $(cat gc3temp_actinopteri); do
		currnum=$(echo "($gc3 - $averagegc_actinopteri) ^ 2" | bc -l)
		sommsq=$(echo "$sommsq + $currnum" | bc -l)
	done
	vargc_actinopteri=$(echo "$sommsq / ($nsp_actinopteri-1)" | bc -l)
fi

# for aves
if [ ! -s gc3temp_aves ]; then
        sumgc_aves=0
        summissingdata_aves=0
        nsp_aves=0
        averagegc_aves="NA"
        averagemissingdata_aves="NA"
	vargc_aves="NA"
else
        sumgc_aves=$(awk '{print $1}' gc3temp_aves | paste -sd+ | bc)
        summissingdata_aves=$(awk '{print $1}' missingdatatemp_aves | paste -sd+ | bc)
        nsp_aves=$(grep -vE -f regexdollar_$1 $ave_pref$gene$suff | grep -c '>')
        averagegc_aves=$(echo "$sumgc_aves / $nsp_aves" | bc -l)
        averagemissingdata_aves=$(echo "$summissingdata_aves / $nsp_aves" | bc -l)
        sommsq=0
        for gc3 in $(cat gc3temp_aves); do
                currnum=$(echo "($gc3 - $averagegc_aves) ^ 2" | bc -l)
                sommsq=$(echo "$sommsq + $currnum" | bc -l)
        done
        vargc_aves=$(echo "$sommsq / ($nsp_aves-1)" | bc -l)
fi

# for insecta
if [ ! -s gc3temp_insecta ]; then
        sumgc_insecta=0
        summissingdata_insecta=0
        nsp_insecta=0
        averagegc_insecta="NA"
        averagemissingdata_insecta="NA"
        vargc_insecta="NA"
else
        sumgc_insecta=$(awk '{print $1}' gc3temp_insecta | paste -sd+ | bc)
        summissingdata_insecta=$(awk '{print $1}' missingdatatemp_insecta | paste -sd+ | bc)
        nsp_insecta=$(grep -vE -f regexdollar_$1 $ins_pref$gene$suff | grep -c '>')
        averagegc_insecta=$(echo "$sumgc_insecta / $nsp_insecta" | bc -l)
        averagemissingdata_insecta=$(echo "$summissingdata_insecta / $nsp_insecta" | bc -l)
        sommsq=0
        for gc3 in $(cat gc3temp_insecta); do
                currnum=$(echo "($gc3 - $averagegc_insecta) ^ 2" | bc -l)
                sommsq=$(echo "$sommsq + $currnum" | bc -l)
        done
        vargc_insecta=$(echo "$sommsq / ($nsp_insecta-1)" | bc -l)
fi

# for mammalia
if [ ! -s gc3temp_mammalia ]; then
	sumgc_mammalia=0
	summissingdata_mammalia=0
	nsp_mammalia=0
	averagegc_mammalia="NA"
	averagemissingdata_mammalia="NA"
	vargc_mammalia="NA"
else
	sumgc_mammalia=$(awk '{print $1}' gc3temp_mammalia | paste -sd+ | bc)
	summissingdata_mammalia=$(awk '{print $1}' missingdatatemp_mammalia | paste -sd+ | bc)
	nsp_mammalia=$(grep -vE -f regexdollar_$1 $mam_pref$gene$suff | grep -c '>')
	averagegc_mammalia=$(echo "$sumgc_mammalia / $nsp_mammalia" | bc -l)
	averagemissingdata_mammalia=$(echo "$summissingdata_mammalia / $nsp_mammalia" | bc -l)
	sommsq=0
	for gc3 in $(cat gc3temp_mammalia); do
		currnum=$(echo "($gc3 - $averagegc_mammalia) ^ 2" | bc -l)
		sommsq=$(echo "$sommsq + $currnum" | bc -l)
	done
	vargc_mammalia=$(echo "$sommsq / ($nsp_mammalia-1)" | bc -l)
fi

# for mollusca

if [ ! -s gc3temp_mollusca ]; then
	sumgc_mollusca=0
	summissingdata_mollusca=0
	nsp_mollusca=0
	averagegc_mollusca="NA"
	averagemissingdata_mollusca="NA"
	vargc_mollusca="NA"
else
	sumgc_mollusca=$(awk '{print $1}' gc3temp_mollusca | paste -sd+ | bc)
	summissingdata_mollusca=$(awk '{print $1}' missingdatatemp_mollusca | paste -sd+ | bc)
	nsp_mollusca=$(grep -vE -f regexdollar_$1 $mol_pref$gene$suff | grep -c '>')
	averagegc_mollusca=$(echo "$sumgc_mollusca / $nsp_mollusca" | bc -l)
	averagemissingdata_mollusca=$(echo "$summissingdata_mollusca / $nsp_mollusca" | bc -l)
	sommsq=0
	for gc3 in $(cat gc3temp_mollusca); do
		currnum=$(echo "($gc3 - $averagegc_mollusca) ^ 2" | bc -l)
		sommsq=$(echo "$sommsq + $currnum" | bc -l)
	done
	vargc_mollusca=$(echo "$sommsq / ($nsp_mollusca-1)" | bc -l)
fi

# Average GC3, variance and missing data overall

sumgc_total=$( echo "$sumgc_actinopteri + $sumgc_aves + $sumgc_insecta + $sumgc_mammalia + $sumgc_mollusca" | bc -l)

summissingdata_total=$(echo "$summissingdata_actinopteri + $summissingdata_aves + $summissingdata_insecta + $summissingdata_mammalia + $summissingdata_mollusca" | bc -l)

nsp_total=$(echo "$nsp_actinopteri + $nsp_aves + $nsp_insecta + $nsp_mammalia + $nsp_mollusca" | bc -l)

averagegc_total=$(echo "$sumgc_total" / "$nsp_total" | bc -l)

averagemissingdata_total=$(echo "$summissingdata_total" / "$nsp_total" | bc -l)

cat gc3temp_* > all_gc3temp
sommsq=0
for gc3 in $(cat all_gc3temp); do
	currnum=$(echo "($gc3 - $averagegc_total) ^ 2" | bc -l)
	sommsq=$(echo "$sommsq + $currnum" | bc -l)
done

vargc_total=$(echo "$sommsq / ($nsp_total-1)" | bc -l)

echo -e "$gene\t$nsp_total\t$total_sp\t$averagemissingdata_total\t$averagegc_total\t$vargc_total\t$nsp_actinopteri\t$actinopteri_sp\t$averagemissingdata_actinopteri\t$averagegc_actinopteri\t$vargc_actinopteri\t$nsp_aves\t$aves_sp\t$averagemissingdata_aves\t$averagegc_aves\t$vargc_aves\t$nsp_insecta\t$insecta_sp\t$averagemissingdata_insecta\t$averagegc_insecta\t$vargc_insecta\t$nsp_mammalia\t$mammalia_sp\t$averagemissingdata_mammalia\t$averagegc_mammalia\t$vargc_mammalia\t$nsp_mollusca\t$mollusca_sp\t$averagemissingdata_mollusca\t$averagegc_mollusca\t$vargc_mollusca" >> gc3_report_${thresh}.tsv

done

rm all_gc3temp gc3temp_* missingdatatemp_* regex*
