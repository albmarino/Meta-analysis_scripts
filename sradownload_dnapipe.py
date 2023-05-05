#!/usr/bin/env python3 
# -*- coding: utf-8 -*
# python3 download_sra.py dnaPT_annotation.tsv

import sys,os,re, time, csv
import subprocess


### variable
table=sys.argv[1]
dict_sra={} ### dict_sra[experiment]=[species,Unified_Cvalues_methods,]
initial_folder="/media/bigvol/alba/downloads/run1/"
dnapt_folder="/media/bigvol/alba/bin/pipeline_dnapipe/"
output_folder="/media/bigvol/alba/TEannotation/dnapipeTE_2rounds_outputs/"


### parse table
with open(table,newline='') as csvfile:
	reader = csv.DictReader(csvfile, delimiter='\t')
	for row in reader:
		sra_acc=row['Experiment']
		species=row['Species'].replace(" ", "_")
		size=row['Unified_Cvalues_methods']
		if sra_acc != "":
			if sra_acc not in dict_sra.keys():
				dict_sra[sra_acc]=[species,size]
			else:
				print(f"oups we have two {sra_acc}")
				break


### Download data drom sra and put them in the good folder
for sra_id in dict_sra.keys():
	species=dict_sra[sra_id][0]
	size=dict_sra[sra_id][1]
	sr_folder=f"{initial_folder}{species}"
	sr_folder2=f"{initial_folder}{species}/short_reads"
	if not os.path.isdir(sr_folder):
		os.mkdir(sr_folder)
	if not os.path.isdir(sr_folder2):
		os.mkdir(sr_folder2)
	print("Currently downloading: " + sra_id)
	cmd_prefetch = f"cd {sr_folder2} && prefetch {sra_id} --max-size 100GB -f yes -p"
	subprocess.run(cmd_prefetch, shell=True)
	print ("Generating fastq for: " + sra_id)
	# slow fastq convertion
	#cmd_fastq= f"cd {sr_folder2}/ && fastq-dump --outdir {sr_folder2} --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip {sr_folder2}/{sra_id}/{sra_id}.sra"
	#subprocess.run(cmd_fastq, shell=True)
	#print(f"done for {sra_id}")

	# WAY faster fastq convertion - could include --print-read-nr but the option is not recognized and fasterq-dump crushes
	cmd_fasterq=f"cd {sr_folder2}/ && fasterq-dump --outdir {sr_folder2} --mem 1G --split-3 --threads 8 --skip-technical {sr_folder2}/{sra_id}/{sra_id}.sra"
	subprocess.run(cmd_fasterq, shell=True)
	print(f"piggingz {sra_id}")
	subprocess.run(f"pigz {sr_folder2}/*fastq", shell=True)
	print(f"fastq created for {sra_id}")


### Run dnaPT pipeline on the downloaded fastq

	print(f"Running 2-rounds dnaPT on {sra_id} --> {species}")
	if os.path.exists(f"{sr_folder2}/{sra_id}.fastq.gz"): # if the fastq is not split into 1 and 2 rename it accordingly
		cmd_rename = (f"mv {sr_folder2}/{sra_id}.fastq.gz {sr_folder2}/{sra_id}_1.fastq.gz")
		subprocess.run(cmd_rename, shell=True)
	cmd_rundnapt = f'cd {dnapt_folder} && snakemake all --use-conda -j 8 -C genome_size={size} sampling_size="0.25" out_dir={output_folder}{species} short_reads={sr_folder2}/{sra_id}"_1.fastq.gz" species={species} acc_num={sra_id}'
	subprocess.run(cmd_rundnapt, shell=True)

### If finished successfully, remove trimmed reads and move to backup folder

	if os.path.exists(f"{output_folder}{species}/final_dnapipete_output/dnaPipeTE_landscapes_subclass.pdf"):
		subprocess.run(f"rm -r {sr_folder}", shell=True)
		subprocess.run(f"rm -r {output_folder}{species}/trim_reads", shell=True)
		subprocess.run(f"mv {output_folder}{species} {output_folder}run1_2backup/", shell=True)
