#!/bin/bash

# Command to run coevol 

# Required input:
# -d = concatenated genes (concatenate_4coevol.sh output)
# -t = clade phylogeny
# -c = characters matrix (make_coevol_table.R output)
# run prefix

~/bin/coevol -d concatenate_gcpoorset_actinopteri.fasta -t actinopteri_brlen_rooted.treefile -fixtimes -c Actinopteri_coevol.txt -dsom actinopteri_gcpoor_dsom


