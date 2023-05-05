#!/bin/bash
# Sum up the overall TE content for each species (% genome coverage and bp, overall and by TE subclass) and the recent TE content (% of read counts over the total number - defined as % genome - and bp, overall and by TE subclass) as defined by $1
# Usage: ./dnaPT_summary.sh 5 recentTEs.tsv overallTEs.tsv
# Arguments: maximum divergence to consider TEs as recent, output filename for recent TEs info, output filename for overall TEs info
# launch from the working directory containing only the output directories of the dnaPT pipeline

# Prepare the univocal list of subclasses and proceed only if we have all the subclasses taken into account, otherwise see message

awk '{print$6}' */final_dnapipete_output/reads_per_component_and_annotation | grep -v "^$"|cut -d'/' -f1|sort -u | grep -v "RNA\|Satellite\|Simple_repeat\|Low_complexity\|ARTEFACT" > subclasslist
echo -e "Overall_TEs\nOther" >> subclasslist
sort -o subclasslist subclasslist
nrow=$(($(cat subclasslist|wc -l)))

if (( nrow < 10 )); then

	echo "Subclasses sampled in this round do not include all the known ones: check subclasslist"

elif (( nrow > 10 )); then

	echo "Subclasses sampled in this round are more than the known ones: check subclasslist"

elif (( nrow = 10 )); then

	echo "Subclasses sampled in this round as many as expected: continuing..."
	startpath=$PWD

	if [ ! -f $2 ]; then
		echo -n "" > $2 # if not existing, create new empty output file for recent TE for this session
		sed 's/$/_recent/g' subclasslist | tr "\n" "\t" >> $2
        	sed -i 's/\t$/\n/g;s/^/Species\t/g' $2
	else
		echo "$2 already exists and won't be overwritten, exiting..."
		exit
	fi

	if [ ! -f $3 ]; then
		echo -n "" > $3 # if not existing, create new file for overall TE content
		echo -e "Species\tOverall_repeat_bp\tOverall_repeat_gsperc\tdnaPT_DNA_bp\tdnaPT_DNA_gsperc\tdnaPT_RC_bp\tdnaPT_RC_gsperc\tdnaPT_LTR_bp\tdnaPT_LTR_gsperc\tdnaPT_PLE_bp\tdnaPT_PLE_gsperc\tdnaPT_LINE_bp\tdnaPT_LINE_gsperc\tdnaPT_SINE_bp\tdnaPT_SINE_gsperc\tdnaPT_Retroposon_bp\tdnaPT_Retroposon_gsperc\tdnaPT_Unknown_bp\tdnaPT_Unknown_gsperc\tdnaPT_Other_bp\tdnaPT_Other_gsperc" > $3
	else
		echo "$3 already exists and won't be overwritten, exiting..."
		exit
	fi

# break if either $2 or $3 files already exist, otherwise keep on

	for dir in "$startpath"/*/; do
		dir=${dir%\/}
		species=${dir##*/}
		finaldir="$startpath"/$species/final_dnapipete_output
		annotdir="$startpath"/$species/final_dnapipete_output/Annotation
		cp subclasslist "$annotdir"/
# Only use Counts.txt to get the total bp (will be 0.25x the given genome size)

		sed -i "/^$/d" "$finaldir"/Counts.txt
		Total_bp=$(grep "Total" "$finaldir"/Counts.txt | cut -f2)

#and use instead reads_per_component_per_annotation to get the proportions of each subclass
# Overall composition of repeat subclasses in base pairs
		awk '{FS=" "} {gsub(/\/.*/, "", $6)}1 {print$6}' "$finaldir"/reads_per_component_and_annotation | awk '{if ($1=="" || $1=="Unknown") {print "Unknown"}
		else if ($1=="DNA") {print $1}
		else if ($1=="LINE") {print$1}
		else if ($1=="LTR") {print$1}
		else if ($1=="PLE") {print$1}
		else if ($1=="RC") {print$1}
		else if ($1=="Retroposon") {print$1}
		else if ($1=="SINE") {print$1}
		else if ($1=="Low_complexity" || $1=="rRNA" || $1=="Satellite" || $1=="Simple_repeat" || $1=="snRNA" ||
		$1=="tRNA" || $1=="srpRNA" || $1=="ARTEFACT") {print "Other_repeats"}}' > "$finaldir"/subclass_col

		paste -d ' ' "$finaldir"/reads_per_component_and_annotation "$finaldir"/subclass_col > "$finaldir"/reads_per_component_and_annotation_subclass

		LTR_bp=$(awk '{FS=" "} $NF=="LTR" {print$2}' "$finaldir"/reads_per_component_and_annotation_subclass |paste -sd+ | bc)
		if [ -z "$LTR_bp" ]; then LTR_bp=$((0)); fi
		LINE_bp=$(awk '{FS=" "} $NF=="LINE" {print$2}' "$finaldir"/reads_per_component_and_annotation_subclass |paste -sd+ | bc)
		if [ -z "$LINE_bp" ]; then LINE_bp=$((0)); fi
		SINE_bp=$(awk '{FS=" "} $NF=="SINE" {print$2}' "$finaldir"/reads_per_component_and_annotation_subclass |paste -sd+ | bc)
		if [ -z "$SINE_bp" ]; then SINE_bp=$((0)); fi
		DNA_bp=$(awk '{FS=" "} $NF=="DNA" {print$2}' "$finaldir"/reads_per_component_and_annotation_subclass |paste -sd+ | bc)
		if [ -z "$DNA_bp" ]; then DNA_bp=$((0)); fi
		PLE_bp=$(awk '{FS=" "} $NF=="PLE" {print$2}' "$finaldir"/reads_per_component_and_annotation_subclass |paste -sd+ | bc)
		if [ -z "$PLE_bp" ]; then PLE_bp=$((0)); fi
		RC_bp=$(awk '{FS=" "} $NF=="RC" {print$2}' "$finaldir"/reads_per_component_and_annotation_subclass |paste -sd+ | bc)
		if [ -z "$RC_bp" ]; then RC_bp=$((0)); fi
		Retroposon_bp=$(awk '{FS=" "} $NF=="Retroposon" {print$2}' "$finaldir"/reads_per_component_and_annotation_subclass |paste -sd+ | bc)
		if [ -z "$Retroposon_bp" ]; then Retroposon_bp=$((0)); fi
		Other_bp=$(awk '{FS=" "} $NF=="Other_repeats" {print$2}' "$finaldir"/reads_per_component_and_annotation_subclass |paste -sd+ | bc)
		if [ -z "$Other_bp" ]; then Other_bp=$((0)); fi
		Unknown_bp=$(awk '{FS=" "} $NF=="Unknown" {print$2}' "$finaldir"/reads_per_component_and_annotation_subclass |paste -sd+ | bc)
		if [ -z "$Unknown_bp" ]; then Unknown_bp=$((0)); fi
		Overall_repeat_bp=$(awk '{FS=" "} {print$2}' "$finaldir"/reads_per_component_and_annotation_subclass | paste -sd+ - | bc)

# Overall composition of repeat subclasses as % of genome size (subclass bp/total bp * 100)

		LTR_perc=$(echo "scale=4 ; $LTR_bp / $Total_bp * 100" | bc)
		LINE_perc=$(echo "scale=4 ; $LINE_bp / $Total_bp * 100" | bc)
		SINE_perc=$(echo "scale=4 ; $SINE_bp / $Total_bp * 100" | bc)
		DNA_perc=$(echo "scale=4 ; $DNA_bp / $Total_bp * 100" | bc)
		PLE_perc=$(echo "scale=4 ; $PLE_bp / $Total_bp * 100" | bc)
		RC_perc=$(echo "scale=4 ; $RC_bp / $Total_bp * 100" | bc)
		Retroposon_perc=$(echo "scale=4 ; $Retroposon_bp / $Total_bp * 100" | bc)
		Unknown_perc=$(echo "scale=4 ; $Unknown_bp / $Total_bp * 100" | bc)
		Other_perc=$(echo "scale=4 ; $Other_bp / $Total_bp * 100" | bc)
		Overall_repeat_perc=$(echo "scale=4 ; $Overall_repeat_bp / $Total_bp * 100" | bc)

		echo -e "$species\t$Overall_repeat_bp\t$Overall_repeat_perc\t$DNA_bp\t$DNA_perc\t$RC_bp\t$RC_perc\t$LTR_bp\t$LTR_perc\t$PLE_bp\t$PLE_perc\t$LINE_bp\t$LINE_perc\t$SINE_bp\t$SINE_perc\t$Retroposon_bp\t$Retroposon_perc\t$Unknown_bp\t$Unknown_perc\t$Other_bp\t$Other_perc" >> $3

		rm "$finaldir"/subclass_col "$finaldir"/reads_per_component_and_annotation_subclass
		cd "$annotdir"

		Rscript ~/bin/dnapt_recentTEs.R $1 $species "$startpath" $2

		cd "$startpath"

	done

fi
