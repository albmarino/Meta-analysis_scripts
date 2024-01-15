#!/bin/bash

# Command to run iqtree

# Required input: concatenated AA sequences generated with concatenate.sh
# Example usage: bash iqtree.sh actinopteri_out_107.fst

iqtree -s ${1} -st AA -m JTT+F+R10 -nt AUTO -ntmax 20 -bb 1000