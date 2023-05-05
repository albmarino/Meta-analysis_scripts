#!/usr/bin/env Rscript

# To run: Rscript excludeSpeciesTree.R falserm_mollusca_561943at33208_final_mask_align_NT.aln mollusca_out_107.fst.treefile

args = commandArgs(trailingOnly=TRUE)
library(ape)

mytree<-read.tree(args[1])
#alns<-list.files(path=".",pattern="^falserm_*NT.aln$", all.files=FALSE, full.names=FALSE)


al<-read.dna(args[1], format="fasta", as.character=T, as.matrix=F)
sp2rm<-	mytree$tip.label[! mytree$tip.label %in% names(al)] # store species to remove by excluding the tree labels that are not in the alignment

print(sp2rm)

prunedtree<-drop.tip(mytree, sp2rm) # drop those species from the tree and write down the "pruned" tree

write.tree(prunedtree, file=paste0(args[1], ".treefile")) # outputs mollusca_561943at33208_final_mask_align_NT.aln.treefile

