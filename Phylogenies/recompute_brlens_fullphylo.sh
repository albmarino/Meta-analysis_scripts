#!/bin/bash

# Example script to recompute branch lengths of the full topology with the aligned profiles

# Required input: full tree topology (produced as in merge_phylogenies_example.sh) as 1st argument, full profile alignments of the 50 markers (fulltree* prefix)
# Example run: bash recompute_brlens_fullphylo.sh fulltree.treefile
# Outputs the full phylogeny with branch lenghts based on the 50 markers

grep '>' fulltree*AA.fasta | cut -d: -f2 | sort -u | sed 's/>//' > list_esp.txt
ls fulltree*AA.fasta > list_file.txt

Concatenate list_file.txt fulltree_aligned.fasta list_esp.txt AA Fasta

# Estimate branch lengths on the full phylogeny topology
iqtree -s fulltree_aligned.fasta -st AA -m LG+G  -te $1 -pre FULLTREE_brlengths