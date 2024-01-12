#! /usr/bin/env python

# this script takes busco output directories from one genome and extract nucleotide fasta files of single copy orhtologs
# example run:   python3 faatofa_singlecopy_busco_editgenecol.py /your/full/path/to/dir/with/buscodirs

import sys
import os
from Bio import SeqIO

parentpath = sys.argv[1] # argument must be the parent directory containing all the busco directories WITH ITS FULL PATH
list_buscodir = os.listdir(parentpath) # make sure the input parent path only contains the directories with busco results

for buscodir in list_buscodir:

	with open("%s.summary.out" % buscodir, "a") as outsummary: # output file for reporting steps
		pathto_buscodir = os.path.join(parentpath, buscodir) # get the full path to the busco output directory (e.g. `/media/bigvol/alba/busco_output/GCA_000003815.2.fa`)
		currdir = os.path.basename(pathto_buscodir)
		outsummary.write("Processing " + currdir + "\n")

		path_tofaa = "run_metazoa_odb10/busco_sequences/single_copy_busco_sequences/"
		faafiles_fullpath = os.path.join(pathto_buscodir, path_tofaa) # get the full path to the directory containing all the faa files of single busco genes (e.g. `/media/bigvol/alba/busco_output/GCA_000003815.2.fa/run_metazoa_odb10/busco_sequences/single_copy_busco_sequences/`)
		path_tometaeuk_initial = "run_metazoa_odb10/metaeuk_output/initial_results/"
		metaeuk_initial_fullpath = os.path.join(pathto_buscodir, path_tometaeuk_initial, buscodir) # get the full path basis for .gff and .codon.fas files in the initial_results folder
		path_tometaeuk_rerun = "run_metazoa_odb10/metaeuk_output/rerun_results/"
		metaeuk_rerun_fullpath = os.path.join(pathto_buscodir, path_tometaeuk_rerun, buscodir) # get the full path basis for .gff and .codon.fas files in the rerun_results folder
		list_faafiles = os.listdir(faafiles_fullpath) # make a list containing all the .faa filenames of the current buscodir

		if list_faafiles:
			outsummary.write("List of faa filenames created for " + currdir + "\n")

			with open("%s.singlecopy.fasta" % buscodir, "w") as outfasta:

				outsummary.write("Processing faa files for " + currdir + "\n")
				records = [] # create empty list for storing matching records of single copy sequences
				nomatch_genenames_list = [] # create empty list for storing non matching records of single copy sequences
				for faafile in list_faafiles: # for every .faa filename in the current buscodir
					faafile_fullpath = os.path.join(faafiles_fullpath, faafile) # get its full path

					with open(faafile_fullpath, "r") as in_faa: # and open it

						genename=os.path.basename(faafile_fullpath) # get the busco gene name

						start_end_strlist = in_faa.readline().partition(":")[2].rstrip().split("-") # read the header, extract the position information after the column, delete the newline, and split the strings separated by a dash in a two-objects list
						start_end_intlist = [int(pos) for pos in start_end_strlist] # convert the list of string of characters to a list of integers
						start_gff,end_gff = start_end_intlist[0] + 1, start_end_intlist[1] + 1 # start and end position of the gene in the gff (incremented by 1)
						start_nt, end_nt = start_end_intlist[0], start_end_intlist[1] # start and end position in the nucleotide sequence
						start_gff,end_gff = str(start_gff), str(end_gff) # reconvert the list objects to strings for matching them to strings from the gff file
						start_nt,end_nt = str(start_nt), str(end_nt) # reconvert the list of objects to strings for matching them to strings from the codon.fas file

						idseq = "" # create empty string for storing the ID of the single copy gene

						with open("%s.gff" % metaeuk_initial_fullpath, "r") as gff_initial, open("%s.codon.fas" % metaeuk_initial_fullpath, "r") as codon_initial:
							for line in gff_initial: # loop through all the lines (annotations) of the gff in the initial_results folder
								if line.split("\t")[2] == "gene" and line.split("\t")[3] == start_gff and line.split("\t")[4] == end_gff and genename.replace(".faa","") in line.split("\t")[-1]: # if the desired features and positions of the current busco gene are matched,
									idseq = line.split("\t")[-1].partition("=")[2].partition(";")[0] # save the sequence ID
									for record in SeqIO.parse(codon_initial, "fasta"): # loop through every record of the fasta file containing the sequences of the initial_results,
										if idseq in record.id and start_nt in record.id and end_nt in record.id: # and if the header contains the saved ID
											records.append(record) # append the record to the records list
											outsummary.write(genename + " sequence retrieved in initial_results" + "\n")

						if not idseq: # if the ID was not retrieved in the initial_results folder
							with open("%s.gff" % metaeuk_rerun_fullpath, "r") as gff_rerun, open("%s.codon.fas" % metaeuk_rerun_fullpath, "r") as codon_rerun:
								for line in gff_rerun: # loop through all the lines (annotations) of the gff in the rerun_results folder
									if line.split("\t")[2] == "gene" and line.split("\t")[3] == start_gff and line.split("\t")[4] == end_gff and genename.replace(".faa","") in line.split("\t")[-1]: # if the desired features and positions of the current busco gene are matched,
										idseq = line.split("\t")[-1].partition("=")[2].partition(";")[0] # save the sequence ID
										for record in SeqIO.parse(codon_rerun, "fasta"): # loop through every record of the fasta file containing the sequences of the rerun_results,
											if idseq in record.id and start_nt in record.id and end_nt in record.id: # and if the header contains the saved ID
												records.append(record) # append the record to the records list
												outsummary.write(genename + " sequence retrieved in rerun_results" + "\n")

						if not idseq: # if after looking in rerun_results still no match is found,
							outsummary.write(genename + " sequence not retrieved" + "\n")
							nomatch_genenames_list.append(genename)
							with open("%s.singlecopy_missing_list" % buscodir, "a") as nomatch_summary:
								nomatch_summary.write(genename + "\n")

				SeqIO.write(records, outfasta, "fasta") # write all the retrieved records of single copy genes in the current output file
				outsummary.write("All the retrieved single copy genes for " + buscodir + " are in " + outfasta.name + "\n")
				outsummary.write(buscodir + " had " + str(len(nomatch_genenames_list)) + " missing single copy sequences;" + " if present, list of missing genes reported in " + "%s.singlecopy_missing_list" % buscodir + "\n")

		outsummary.write("\nFinished!")
