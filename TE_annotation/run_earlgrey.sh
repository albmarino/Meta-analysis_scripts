#!/bin/bash

# script for TE annotation with EarlGrey
# run example: bash run_earlgrey.sh GCA_001542645.1.fa Anopheles_gambiae

earlGrey -g $1 -s $2 -r metazoa -o ./"${2}"_earlgrey -t 8