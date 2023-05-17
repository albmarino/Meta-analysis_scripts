## Extracting single copy BUSCO genes
The nucleotide sequences of single copy genes were extracted from BUSCO outputs with `faatofa_singlecopy_busco_editgenecol.py`.

## Gene alignments and phylogeny
The 954 BUSCO genes were grouped according to the taxa Aves, Actinopteri, Insecta, Mammalia, Mollusca. Alignments and phylogenies were first computed clade-wise with the OMM_MACSE pipeline (version 11.05) employing MACSE 2.06 (Ranwez et al., 2018; Scornavacca et al., 2019) and IQ-TREE 1.6.12 (Nguyen et al., 2015), and then progressively merged into a whole phylogeny.

## Aligning BUSCO genes
954 gene alignments were computed for each clade as below:
```
singularity run ~/bin/omm_macse_v11.05b.sif --out_dir ./actinopteri_${gene}_aligned --out_file_prefix actinopteri_${gene} --in_seq_file ${gene}_prealign_species.fasta --genetic_code_number 1 --java_mem 20000m
```
