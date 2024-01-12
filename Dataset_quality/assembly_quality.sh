#!/bin/bash

# Script to run BUSCO and Quast on a directory containing genomes

# run: buscorun.sh assemblydir/

busco -i $1 -o buscoout -l metazoa_odb10 --update-data -m genome --cpu 30

for genome in $(ls $1); do

	python3 quast.py $genome -e -s -o "${genome}"_quast_report

done